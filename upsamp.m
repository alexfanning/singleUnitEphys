function sampRateChange = upsamp(var2scale,var2scale2)
%
%   Upsamples data using the built-in matlab function interp1
%
%   Alex Fanning, December 2019
% *************************************************************************

scale = (length(var2scale2)/length(var2scale));
x = [1 scale*(2:length(var2scale))];
xq = 1:length(var2scale2);
sampRateChange = interp1(x, var2scale, xq);