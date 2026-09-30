function summary = summarize_inversions(M, runs, varargin)
%SUMMARIZE_INVERSIONS Step 3: group runs into stress solutions and validate them.
%   summary = summarize_inversions(M, runs)
%   summary = summarize_inversions(M, runs, 'Friction', 0.6, 'Alpha', 0.05)
%
%   M    : output of compute_nodal_planes
%   runs : output of run_stress_inversions
%   summary.stress_solutions      one row per distinct set of chosen planes,
%                                 ranked by how many runs reached it
%   summary.solution_planes       the chosen plane of every event, per solution
%   summary.nodal_plane_frequency how often each plane of each event was chosen
%
%   If M has uncertainty_deg, each solution gets a variance test of
%   misfit / uncertainty against a variance of 1 (same as vartest).
p = inputParser;
p.addParameter('Friction', [], @(x) isempty(x) || isnumeric(x));
p.addParameter('Alpha', 0.05, @isnumeric);
p.parse(varargin{:});
friction = p.Results.Friction;
if isempty(friction), friction = runs.friction(1); end
alpha = p.Results.Alpha;

M.event_id = as_text_ids(M.event_id);
planes = plane_set(M);
has_unc = ismember('uncertainty_deg', M.Properties.VariableNames);
n_runs = height(runs);
n_events = height(M);
choices_text = cellstr(runs.plane_choice);
converged = strcmp(cellstr(runs.status), 'converged');

[keys, ~, group] = unique(choices_text);
group = group(:);
counts = accumarray(group, 1);
n_conv = accumarray(group, double(converged));
first_seen = accumarray(group, (1:n_runs)', [], @min);
% Most runs first; ties broken by which solution appeared first.
[~, order] = sortrows([-counts, first_seen]);
n_sol = numel(keys);

S = zeros(n_sol, 13); tests = nan(n_sol, 3);
plane_rows = cell(n_sol, 1);
for s = 1:n_sol
    g = order(s);
    choice = decode_choice(keys{g});
    [normals, slips] = select_planes(planes, choice);
    result = invert_stress(normals, slips);
    inst = fault_instability(normals, result, friction);
    S(s, :) = [s, counts(g), counts(g) / n_runs, n_conv(g), ...
        reshape(result.axes_trend_plunge', 1, []), result.shape_ratio, ...
        mean(result.misfit_deg), std(result.misfit_deg)];

    strike = M.strike_1; dip = M.dip_1; rake = M.rake_1;
    two = choice == 2;
    strike(two) = M.strike_2(two); dip(two) = M.dip_2(two); rake(two) = M.rake_2(two);
    P = table(repmat(s, n_events, 1), M.event_id, choice, strike, dip, rake, ...
        round_to(result.misfit_deg, 4), round_to(inst, 4), 'VariableNames', ...
        {'solution_id', 'event_id', 'chosen_plane', 'strike', 'dip', 'rake', 'misfit_deg', 'instability'});
    if has_unc
        normalized = result.misfit_deg ./ M.uncertainty_deg;
        [stat, pval] = variance_test(normalized, 1);
        tests(s, :) = [stat, pval, pval >= alpha];
        P.normalized_misfit = round_to(normalized, 4);
    end
    plane_rows{s} = P;
end

names = {'solution_id', 'n_runs', 'fraction_of_runs', 'n_runs_converged', ...
    'sigma1_trend', 'sigma1_plunge', 'sigma2_trend', 'sigma2_plunge', ...
    'sigma3_trend', 'sigma3_plunge', 'shape_ratio', 'misfit_mean_deg', 'misfit_std_deg'};
sol = array2table(round_to(S, 4), 'VariableNames', names);
if has_unc
    sol.variance_test_statistic = round_to(tests(:, 1), 4);
    sol.variance_test_p_value = round_to(tests(:, 2), 4);
    sol.fits_within_uncertainty = logical(tests(:, 3));
end
sol.plane_choice = keys(order);
summary.stress_solutions = sol;
summary.solution_planes = vertcat(plane_rows{:});

% How often each nodal plane was chosen
C = cell2mat(cellfun(@(k) decode_choice(k)', choices_text, 'UniformOutput', false));
times1 = sum(C == 1, 1)'; times2 = sum(C == 2, 1)';
idx = reshape([1:n_events; 1:n_events], [], 1);
plane = repmat([1; 2], n_events, 1);
times = reshape([times1'; times2'], [], 1);
pref = reshape([(times1 >= times2)'; (times1 < times2)'], [], 1);
f_strike = reshape([M.strike_1'; M.strike_2'], [], 1);
f_dip = reshape([M.dip_1'; M.dip_2'], [], 1);
f_rake = reshape([M.rake_1'; M.rake_2'], [], 1);
summary.nodal_plane_frequency = table(M.event_id(idx), plane, f_strike, f_dip, f_rake, ...
    times, round_to(times / n_runs, 4), pref, 'VariableNames', {'event_id', 'plane', ...
    'strike', 'dip', 'rake', 'times_chosen', 'fraction_chosen', 'is_preferred'});
end
