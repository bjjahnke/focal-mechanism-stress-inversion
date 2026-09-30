function files = plot_results(summary, output_dir, max_panels)
%PLOT_RESULTS Lower-hemisphere stereonets (equal-angle) of the step 3 results.
%   files = plot_results(summary, output_dir) writes
%   plots/stress_solutions.png and plots/nodal_plane_frequency.png.
if nargin < 3, max_panels = 4; end
plot_dir = fullfile(output_dir, 'plots');
if ~isfolder(plot_dir), mkdir(plot_dir); end
colors = [0.75 0.22 0.17; 0.18 0.55 0.34; 0.17 0.44 0.73];

sol = summary.stress_solutions;
n = min(max_panels, height(sol));
fig = figure('Visible', 'off', 'Position', [100 100 420 * n 480]);
for s = 1:n
    subplot(1, n, s); draw_stereonet();
    P = summary.solution_planes(summary.solution_planes.solution_id == s, :);
    [normals, ~] = plane_vectors(P.strike, P.dip, P.rake);
    [tr, pl] = to_trend_plunge(normals);
    [x, y] = stereo_project(tr, pl);
    plot(x, y, 'o', 'Color', [0.35 0.35 0.35], 'MarkerFaceColor', [0.35 0.35 0.35], 'MarkerSize', 4);
    h = zeros(1, 3);
    for k = 1:3
        [x, y] = stereo_project(sol.(sprintf('sigma%d_trend', k))(s), sol.(sprintf('sigma%d_plunge', k))(s));
        h(k) = plot(x, y, 's', 'MarkerSize', 11, 'MarkerFaceColor', colors(k, :), 'MarkerEdgeColor', 'k');
    end
    if s == 1
        legend(h, {'\sigma_1', '\sigma_2', '\sigma_3'}, 'Location', 'southoutside', 'Orientation', 'horizontal');
    end
    title({sprintf('Solution %d: %.0f%% of runs', s, 100 * sol.fraction_of_runs(s)), ...
        sprintf('R = %.2f, misfit = %.0f^{\\circ} \\pm %.0f^{\\circ}', sol.shape_ratio(s), ...
        sol.misfit_mean_deg(s), sol.misfit_std_deg(s))}, 'FontWeight', 'normal');
end
files{1} = fullfile(plot_dir, 'stress_solutions.png');
print(fig, files{1}, '-dpng', '-r150'); close(fig);

F = summary.nodal_plane_frequency;
fig = figure('Visible', 'off', 'Position', [100 100 560 500]);
draw_stereonet();
[normals, ~] = plane_vectors(F.strike, F.dip, F.rake);
[tr, pl] = to_trend_plunge(normals);
[x, y] = stereo_project(tr, pl);
scatter(x, y, 36, F.fraction_chosen, 'filled', 'MarkerEdgeColor', 'k');
colormap(flipud(gray(256))); caxis([0 1]);
cb = colorbar; ylabel(cb, 'Fraction of runs the plane was chosen');
title('Poles of both nodal planes for every event', 'FontWeight', 'normal');
files{2} = fullfile(plot_dir, 'nodal_plane_frequency.png');
print(fig, files{2}, '-dpng', '-r150'); close(fig);
end

function draw_stereonet()
hold on; axis equal; axis off;
t = linspace(0, 2*pi, 361);
fill(cos(t), sin(t), 'w', 'EdgeColor', 'k');
for p = 10:10:80
    r = tand((90 - p) / 2);
    plot(r * cos(t), r * sin(t), '-', 'Color', [0.88 0.88 0.88]);
end
for a = 0:30:150
    plot([-sind(a) sind(a)], [-cosd(a) cosd(a)], '-', 'Color', [0.88 0.88 0.88]);
end
text(0, 1.1, 'N', 'HorizontalAlignment', 'center');
text(1.1, 0, 'E', 'HorizontalAlignment', 'center');
text(0, -1.12, 'S', 'HorizontalAlignment', 'center');
text(-1.12, 0, 'W', 'HorizontalAlignment', 'center');
xlim([-1.2 1.2]); ylim([-1.2 1.2]);
end

function [x, y] = stereo_project(trend, plunge)
r = tand((90 - plunge) / 2);
x = r .* sind(trend);
y = r .* cosd(trend);
end
