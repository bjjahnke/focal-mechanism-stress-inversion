function [strike, dip, rake] = vectors_to_plane(normal, slip)
%VECTORS_TO_PLANE Strike, dip, rake (degrees) from normal and slip vectors (NED).
n = normal ./ sqrt(sum(normal.^2, 2));
u = slip ./ sqrt(sum(slip.^2, 2));

% The normal must point up; flipping it also flips the slip vector.
down = n(:, 3) > 0;
n(down, :) = -n(down, :);
u(down, :) = -u(down, :);

dip = acosd(max(min(-n(:, 3), 1), -1));
strike = mod(atan2(-n(:, 1), n(:, 2)) * 180 / pi, 360);
sin_dip = sind(dip);
along_strike = u(:, 1) .* cosd(strike) + u(:, 2) .* sind(strike);
rake = atan2(-u(:, 3) ./ sin_dip, along_strike) * 180 / pi;

% Horizontal plane: strike is undefined, so fix it at 0 and measure rake from North.
horizontal = sin_dip < 1e-10;
strike(horizontal) = 0;
rake(horizontal) = atan2(-u(horizontal, 2), u(horizontal, 1)) * 180 / pi;
end
