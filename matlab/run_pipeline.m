function outputs = run_pipeline(config_path)
%RUN_PIPELINE Run steps 1-3 from a JSON config file and write every output.
%   outputs = run_pipeline('examples/san_emidio_2016/config.json')
%
%   Config keys (paths are relative to the config file's folder):
%     input_catalog       required  CSV with strike, dip, rake
%     output_dir          required  folder for the outputs
%     friction            optional  default 0.6
%     n_runs              optional  default 1000
%     max_iterations      optional  default 20
%     random_seed         optional  default none (results vary run to run)
%     significance_level  optional  default 0.05
%     make_plots          optional  default true
started = tic;
cfg = read_config(config_path);
if ~isfolder(cfg.output_dir), mkdir(cfg.output_dir); end
out = @(name) fullfile(cfg.output_dir, name);

fprintf('Step 1/3: computing nodal planes\n');
M = compute_nodal_planes(cfg.input_catalog);
writetable(M, out('focal_mechanisms.csv'));

fprintf('Step 2/3: %d iterative inversions (friction %.2f)\n', cfg.n_runs, cfg.friction);
runs = run_stress_inversions(M, 'Friction', cfg.friction, 'NumRuns', cfg.n_runs, ...
    'MaxIterations', cfg.max_iterations, 'RandomSeed', cfg.random_seed);
writetable(runs, out('inversion_runs.csv'));

fprintf('Step 3/3: summarizing and validating\n');
summary = summarize_inversions(M, runs, 'Friction', cfg.friction, 'Alpha', cfg.significance_level);
writetable(summary.stress_solutions, out('stress_solutions.csv'));
writetable(summary.solution_planes, out('solution_planes.csv'));
writetable(summary.nodal_plane_frequency, out('nodal_plane_frequency.csv'));

outputs = struct('focal_mechanisms', M, 'inversion_runs', runs, 'summary', summary);
if cfg.make_plots
    plot_results(summary, cfg.output_dir);
end

meta = struct('language', 'MATLAB', 'version', version, ...
    'run_at', datestr(now, 'yyyy-mm-ddTHH:MM:SS'), 'runtime_seconds', round_to(toc(started), 2), ...
    'config', cfg, 'n_events', height(M), 'n_stress_solutions', height(summary.stress_solutions));
fid = fopen(out('run_metadata.json'), 'w');
fprintf(fid, '%s', jsonencode(meta));
fclose(fid);

top = summary.stress_solutions(1, :);
fprintf('Done in %.1fs. Top solution reached by %.0f%% of runs: sigma1 %.0f/%.0f, R = %.2f\n', ...
    toc(started), 100 * top.fraction_of_runs, top.sigma1_trend, top.sigma1_plunge, top.shape_ratio);
end

function cfg = read_config(config_path)
if ~isfile(config_path)
    error('focal_stress:fileNotFound', 'Config file not found: %s', config_path);
end
raw = jsondecode(fileread(config_path));
for key = {'input_catalog', 'output_dir'}
    if ~isfield(raw, key{1})
        error('focal_stress:badConfig', '%s: missing required key ''%s''', config_path, key{1});
    end
end
defaults = struct('friction', 0.6, 'n_runs', 1000, 'max_iterations', 20, ...
    'random_seed', [], 'significance_level', 0.05, 'make_plots', true);
cfg = defaults;
for f = fieldnames(raw)'
    cfg.(f{1}) = raw.(f{1});
end
base = fileparts(config_path);
cfg.input_catalog = resolve(base, cfg.input_catalog);
cfg.output_dir = resolve(base, cfg.output_dir);
end

function p = resolve(base, p)
is_absolute = startsWith(p, '/') || startsWith(p, '\') || ~isempty(regexp(p, '^[A-Za-z]:', 'once'));
if ~is_absolute
    p = fullfile(base, p);
end
end
