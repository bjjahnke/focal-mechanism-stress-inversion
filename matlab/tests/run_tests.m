function run_tests()
%RUN_TESTS Checks for the MATLAB pipeline. Run from the repository root:
%   addpath('matlab'); addpath('matlab/tests'); run_tests
%
% Plain asserts, so it runs in any MATLAB release (and GNU Octave).
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
example = fullfile(root, 'examples', 'san_emidio_2016', 'catalog.csv');
synthetic = fullfile(root, 'examples', 'synthetic');
tests = {
    @() test_auxiliary_plane_is_consistent()
    @() test_known_strike_slip_mechanism()
    @() test_catalog_rejects_bad_dip()
    @() test_catalog_adds_event_ids()
    @() test_matches_original_inversion(example)
    @() test_recovers_synthetic_stress(synthetic)
    @() test_summary_tables(example)
};
names = {'auxiliary plane is consistent', 'known strike-slip mechanism', ...
    'catalog rejects bad dip', 'catalog adds event ids', 'matches original inversion.m', ...
    'recovers synthetic stress', 'summary tables'};
failed = 0;
for i = 1:numel(tests)
    try
        tests{i}();
        fprintf('  PASS  %s\n', names{i});
    catch err
        failed = failed + 1;
        fprintf('  FAIL  %s: %s\n', names{i}, err.message);
    end
end
fprintf('%d passed, %d failed\n', numel(tests) - failed, failed);
if failed > 0
    error('focal_stress:testsFailed', '%d test(s) failed', failed);
end
end

function test_auxiliary_plane_is_consistent()
rand('seed', 1); %#ok<RAND>
n = 100;
T = table(360 * rand(n, 1), 1 + 88 * rand(n, 1), -179 + 358 * rand(n, 1), ...
    'VariableNames', {'strike', 'dip', 'rake'});
M = compute_nodal_planes(T);
% Plane 2 of plane 2 must be plane 1.
T2 = table(M.strike_2, M.dip_2, M.rake_2, 'VariableNames', {'strike', 'dip', 'rake'});
back = compute_nodal_planes(T2);
% Compare as vectors: plane 2 is rounded to 0.01 deg, and strike of a nearly
% flat plane is very sensitive to that rounding.
[n0, u0] = plane_vectors_public(T.strike, T.dip, T.rake);
[n1, u1] = plane_vectors_public(back.strike_2, back.dip_2, back.rake_2);
assert(max(acosd(min(sum(n0 .* n1, 2), 1))) < 0.1, 'plane does not round-trip');
assert(max(acosd(min(sum(u0 .* u1, 2), 1))) < 0.1, 'slip does not round-trip');
% Both planes must give the same P axis.
assert(max(angle_diff(back.p_trend, M.p_trend)) < 0.05, 'P axis differs between planes');
end

function test_known_strike_slip_mechanism()
M = compute_nodal_planes(table(0, 90, 180, 'VariableNames', {'strike', 'dip', 'rake'}));
assert(abs(M.dip_2 - 90) < 1e-6);
assert(abs(mod(M.strike_2, 180) - 90) < 1e-6);
assert(abs(M.p_plunge) < 1e-6);
end

function test_catalog_rejects_bad_dip()
T = table([1; 2; 3], [10; 95; 30], [0; 0; 0], 'VariableNames', {'strike', 'dip', 'rake'});
try
    load_catalog(T);
    error('no error raised');
catch err
    assert(~isempty(strfind(err.message, 'dip')), err.message);
end
end

function test_catalog_adds_event_ids()
T = load_catalog(table([1; 2; 3], [10; 20; 30], [0; 0; 0], 'VariableNames', {'strike', 'dip', 'rake'}));
assert(isequal(T.event_id, {'EV001'; 'EV002'; 'EV003'}));
end

function test_matches_original_inversion(example)
% Reference values from the original inversion.m (Michael 1984 method), run on
% the example catalog with plane 1 for odd events and plane 2 for even events.
M = compute_nodal_planes(example);
even = mod((1:height(M))', 2) == 0;
s = M.strike_1; d = M.dip_1; r = M.rake_1;
s(even) = M.strike_2(even); d(even) = M.dip_2(even); r(even) = M.rake_2(even);
[normals, slips] = plane_vectors_public(s, d, r);
result = invert_stress(normals, slips);
assert(abs(result.shape_ratio - 0.287057) < 1e-5, 'shape ratio differs');
expected = [33.998 69.411; 164.512 13.715; 258.276 15.054];
assert(max(abs(result.axes_trend_plunge(:) - expected(:))) < 2e-3, 'stress axes differ');
assert(abs(mean(result.misfit_deg) - 47.6214) < 1e-3, 'misfit differs');
end

function test_recovers_synthetic_stress(folder)
truth = jsondecode(fileread(fullfile(folder, 'true_answer.json')));
M = compute_nodal_planes(fullfile(folder, 'catalog.csv'));
runs = run_stress_inversions(M, 'NumRuns', 50, 'RandomSeed', 3);
S = summarize_inversions(M, runs);
best = S.stress_solutions(1, :);
sigma1 = trend_plunge_vector(best.sigma1_trend, best.sigma1_plunge);
true1 = trend_plunge_vector(truth.sigma1_trend_plunge(1), truth.sigma1_trend_plunge(2));
assert(acosd(min(abs(sigma1 * true1'), 1)) < 5, 'sigma1 not recovered');
P = S.solution_planes(S.solution_planes.solution_id == 1, :);
true_planes = cellfun(@(id) truth.true_plane_by_event.(id), P.event_id);
assert(mean(P.chosen_plane == true_planes) > 0.9, 'true fault planes not recovered');
end

function test_summary_tables(example)
M = compute_nodal_planes(example);
runs = run_stress_inversions(M, 'NumRuns', 20, 'RandomSeed', 1);
S = summarize_inversions(M, runs);
assert(sum(S.stress_solutions.n_runs) == 20);
assert(height(S.solution_planes) == height(M) * height(S.stress_solutions));
F = S.nodal_plane_frequency;
assert(height(F) == 2 * height(M));
assert(all(abs(F.fraction_chosen(1:2:end) + F.fraction_chosen(2:2:end) - 1) < 1e-9));
assert(ismember('fits_within_uncertainty', S.stress_solutions.Properties.VariableNames));
end

% ---- helpers (the pipeline's own helpers are private to matlab/) ----

function d = angle_diff(a, b)
d = abs(mod(a - b + 180, 360) - 180);
end

function v = trend_plunge_vector(trend, plunge)
v = [cosd(trend) * cosd(plunge), sind(trend) * cosd(plunge), sind(plunge)];
end

function [normal, slip] = plane_vectors_public(strike, dip, rake)
normal = [-sind(dip).*sind(strike), sind(dip).*cosd(strike), -cosd(dip)];
slip = [cosd(rake).*cosd(strike) + cosd(dip).*sind(rake).*sind(strike), ...
        cosd(rake).*sind(strike) - cosd(dip).*sind(rake).*cosd(strike), -sind(rake).*sind(dip)];
end
