function M = compute_nodal_planes(catalog)
%COMPUTE_NODAL_PLANES Step 1: add the second nodal plane and the P, T, N axes.
%   M = compute_nodal_planes(catalog)
%
%   catalog : table or CSV path with strike, dip, rake of one nodal plane.
%   M       : table with event_id, strike_1..rake_1, strike_2..rake_2,
%             p/t/n_trend and p/t/n_plunge, uncertainty_deg (if given),
%             then all other input columns.
%
%   The second plane's normal is the first plane's slip direction, and its
%   slip direction is the first plane's normal.
T = load_catalog(catalog);
[normal1, slip1] = plane_vectors(T.strike, T.dip, T.rake);
[strike2, dip2, rake2] = vectors_to_plane(slip1, normal1);

[p_trend, p_plunge] = to_trend_plunge((normal1 - slip1) / sqrt(2));
[t_trend, t_plunge] = to_trend_plunge((normal1 + slip1) / sqrt(2));
[n_trend, n_plunge] = to_trend_plunge(cross(normal1, slip1, 2));

r = @(x) round_to(x, 2);
new = table(T.event_id, double(T.strike), double(T.dip), double(T.rake), ...
    r(strike2), r(dip2), r(rake2), r(p_trend), r(p_plunge), r(t_trend), r(t_plunge), ...
    r(n_trend), r(n_plunge), 'VariableNames', {'event_id', ...
    'strike_1', 'dip_1', 'rake_1', 'strike_2', 'dip_2', 'rake_2', ...
    'p_trend', 'p_plunge', 't_trend', 't_plunge', 'n_trend', 'n_plunge'});

names = T.Properties.VariableNames;
rest = names(~ismember(names, {'event_id', 'strike', 'dip', 'rake', 'uncertainty_deg'}));
if ismember('uncertainty_deg', names)
    rest = [{'uncertainty_deg'}, rest];
end
M = new;
for k = 1:numel(rest)
    M.(rest{k}) = T.(rest{k});
end
end
