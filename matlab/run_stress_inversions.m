function runs = run_stress_inversions(M, varargin)
%RUN_STRESS_INVERSIONS Step 2: repeated iterative stress inversions (Vavrycuk, 2014).
%   runs = run_stress_inversions(M)
%   runs = run_stress_inversions(M, 'Friction', 0.6, 'NumRuns', 1000, ...
%                                   'MaxIterations', 20, 'RandomSeed', 42)
%
%   M    : output of compute_nodal_planes (needs strike/dip/rake _1 and _2).
%   runs : one row per run: run_id, plane_choice ('1|2|...'), iterations,
%          status, sigma1..3 trend/plunge, shape_ratio, misfit mean/std, friction.
%
%   Each run picks plane 1 or 2 at random for every event, inverts for stress,
%   switches every event to its plane closest to failure, and repeats until
%   the choice stops changing ('converged'), starts repeating ('cycling'),
%   or MaxIterations is reached ('max_iterations').
p = inputParser;
p.addParameter('Friction', 0.6, @(x) isnumeric(x) && x > 0 && x <= 2);
p.addParameter('NumRuns', 1000, @(x) isnumeric(x) && x >= 1);
p.addParameter('MaxIterations', 20, @(x) isnumeric(x) && x >= 1);
p.addParameter('RandomSeed', [], @(x) isempty(x) || isnumeric(x));
p.parse(varargin{:});
opt = p.Results;

require_columns(M, {'event_id', 'strike_1', 'dip_1', 'rake_1', 'strike_2', 'dip_2', 'rake_2'}, ...
    'focal mechanisms (run compute_nodal_planes first)');
n_events = height(M);
if n_events < 3
    error('focal_stress:tooFewEvents', 'Stress inversion needs at least 3 events, got %d', n_events);
elseif n_events < 10
    warning('focal_stress:fewEvents', 'Only %d events; results are more reliable with 10 or more', n_events);
end
if ~isempty(opt.RandomSeed)
    rng(opt.RandomSeed);
end

planes = plane_set(M);
n = opt.NumRuns;
plane_choice = cell(n, 1); status = cell(n, 1);
iterations = zeros(n, 1); axes = zeros(n, 6); shape_ratio = zeros(n, 1);
misfit_mean = zeros(n, 1); misfit_std = zeros(n, 1);
for run = 1:n
    start = 1 + (rand(n_events, 1) < 0.5);
    [choice, result, iterations(run), status{run}] = ...
        iterate_once(planes, start, opt.Friction, opt.MaxIterations);
    plane_choice{run} = encode_choice(choice);
    axes(run, :) = reshape(result.axes_trend_plunge', 1, []);
    shape_ratio(run) = result.shape_ratio;
    misfit_mean(run) = mean(result.misfit_deg);
    misfit_std(run) = std(result.misfit_deg);
end

runs = table((1:n)', plane_choice, iterations, status, ...
    round_to(axes(:, 1), 2), round_to(axes(:, 2), 2), round_to(axes(:, 3), 2), ...
    round_to(axes(:, 4), 2), round_to(axes(:, 5), 2), round_to(axes(:, 6), 2), ...
    round_to(shape_ratio, 4), round_to(misfit_mean, 2), round_to(misfit_std, 2), ...
    repmat(opt.Friction, n, 1), 'VariableNames', {'run_id', 'plane_choice', 'iterations', ...
    'status', 'sigma1_trend', 'sigma1_plunge', 'sigma2_trend', 'sigma2_plunge', ...
    'sigma3_trend', 'sigma3_plunge', 'shape_ratio', 'misfit_mean_deg', 'misfit_std_deg', 'friction'});

fprintf('%d runs: %d converged, %d cycling, %d hit max_iterations; %d distinct plane choices\n', ...
    n, sum(strcmp(status, 'converged')), sum(strcmp(status, 'cycling')), ...
    sum(strcmp(status, 'max_iterations')), numel(unique(plane_choice)));
end
