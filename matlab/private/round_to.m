function y = round_to(x, decimals)
%ROUND_TO Round to a number of decimal places (works in MATLAB and Octave).
y = round(x * 10^decimals) / 10^decimals;
end
