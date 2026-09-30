function result = invert_stress(normals, slips)
%INVERT_STRESS Linear stress inversion of fault-slip data (Michael, 1984).
%   result = invert_stress(normals, slips)
%
%   normals, slips : n-by-3 unit vectors (North-East-Down), one fault per row.
%   result fields  :
%     tensor             3x3 deviatoric stress, compression negative
%     principal_values   [sigma1; sigma2; sigma3], sigma1 most compressive
%     principal_axes     3x3, rows are unit vectors of sigma1, sigma2, sigma3
%     axes_trend_plunge  3x2 [trend plunge] of sigma1, sigma2, sigma3
%     shape_ratio        R = (sigma1 - sigma2) / (sigma1 - sigma3)
%     misfit_deg         n-by-1 angle between slip and predicted shear traction
%     variance           sum of squared residuals / (3n - 5)
n = size(normals, 1);
if n < 2
    error('focal_stress:tooFewFaults', 'Stress inversion needs at least 2 faults');
end
G = zeros(3 * n, 5);
for i = 1:n
    n1 = normals(i, 1); n2 = normals(i, 2); n3 = normals(i, 3);
    G(3*i-2:3*i, :) = [ ...
        n1-n1^3+n1*n3^2,  n2-2*n2*n1^2,  n3-2*n3*n1^2,  -n1*(n2^2-n3^2),   -2*n1*n2*n3; ...
        -n2*(n1^2-n3^2),  n1-2*n1*n2^2,  -2*n1*n2*n3,   n2-n2^3+n2*n3^2,   n3-2*n3*n2^2; ...
        -n3*n1^2-n3+n3^3, -2*n1*n2*n3,   n1-2*n1*n3^2,  -n2^2*n3-n3+n3^3,  n2-2*n2*n3^2];
end
d = reshape(slips', [], 1);
m = G \ d;

tensor = [m(1) m(2) m(3); m(2) m(4) m(5); m(3) m(5) -m(1)-m(4)];
[V, D] = eig((tensor + tensor') / 2);
[values, order] = sort(diag(D));          % ascending: sigma1 first
V = V(:, order);

predicted = reshape(G * m, 3, [])';
cos_fit = sum(predicted .* slips, 2) ./ sqrt(sum(predicted.^2, 2));
[trend, plunge] = to_trend_plunge(V');

result.tensor = tensor;
result.principal_values = values;
result.principal_axes = V';
result.axes_trend_plunge = [trend, plunge];
result.shape_ratio = (values(1) - values(2)) / (values(1) - values(3));
result.misfit_deg = acosd(max(min(cos_fit, 1), -1));
result.variance = sum((G * m - d).^2) / max(3 * n - 5, 1);
end
