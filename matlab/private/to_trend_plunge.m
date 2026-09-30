function [trend, plunge] = to_trend_plunge(v)
%TO_TREND_PLUNGE Trend and plunge (degrees, lower hemisphere) of NED vectors (one per row).
v = v ./ sqrt(sum(v.^2, 2));
up = v(:, 3) < 0;
v(up, :) = -v(up, :);
plunge = asind(max(min(v(:, 3), 1), -1));
trend = mod(atan2(v(:, 2), v(:, 1)) * 180 / pi, 360);
end
