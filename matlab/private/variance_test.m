function [statistic, p_value] = variance_test(x, expected_variance)
%VARIANCE_TEST Two-sided chi-square test that x has the expected variance.
%   Same result as vartest(x, expected_variance) but needs no toolbox.
if nargin < 2, expected_variance = 1; end
x = x(:);
dof = numel(x) - 1;
statistic = dof * var(x) / expected_variance;
cdf = gammainc(statistic / 2, dof / 2);   % chi-square CDF
p_value = min(1, 2 * min(cdf, 1 - cdf));
end
