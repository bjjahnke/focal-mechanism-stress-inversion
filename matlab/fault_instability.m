function I = fault_instability(normals, result, friction)
%FAULT_INSTABILITY Fault instability (Vavrycuk, 2014) of planes under a stress result.
%   I = 1 for a plane optimally oriented for slip, 0 for the most stable plane.
R = result.shape_ratio;
c = normals * result.principal_axes';
c1 = c(:, 1); c2 = c(:, 2); c3 = c(:, 3);
sigma_n = c1.^2 + (1 - 2*R) * c2.^2 - c3.^2;
tau = sqrt(max(c1.^2 + (1 - 2*R)^2 * c2.^2 + c3.^2 - sigma_n.^2, 0));
I = (tau - friction * (sigma_n - 1)) / (friction + sqrt(1 + friction^2));
end
