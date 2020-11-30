function [errbarx, errbary] = MakeErrorBars(x, value, sem);
% function [errbarx, errbary] = MakeErrorBars(x, value, sem);
%
% Make error bars for display using PLOT.
% All arguments are ROWS.

len = length(x);
nans = NaN * ones(1, len);
errbarx = reshape([x; x; nans], 3 * len, 1);
errbary = reshape([value - sem; value + sem; nans], 3 * len, 1);
