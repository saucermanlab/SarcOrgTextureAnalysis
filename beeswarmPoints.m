function [x,y] = beeswarmPoints(yvals,xspacing,yspacing,ylimit)
%% FUNCTION DESCRIPTION
% This function generates scatter points in a beeswarm 'center' layout,
% allowing one to simply call scatter(x,y) after calling this function.
%
% INPUTS
% yvals : array with data points
% xspacing : horizontal spacing between scatter points, highly dependent on
%   figure size and marker size
% yspacing : vertical spacing between scatter points, highly dependent on
%   figure size and marker size
% ylimit : figure axis y-limits, must be set before calling this function,
%   in the form of [ymin ymax]
%
% OUTPUTS
% x : array of x values for scatter points
% y : array of y values for scatter points

%%
x = zeros(1,numel(yvals));

ybins = ylimit(1):yspacing:ylimit(2);
nvals = histcounts(yvals,ylimit(1):yspacing:ylimit(2));
y = repelem(ybins(1:end-1)+yspacing/2,nvals);
ind = 1;
for iY = 1:length(ybins)-1
    if nvals(iY) > 1
        x(ind:(ind+nvals(iY)-1)) = (((-nvals(iY)/2+1)*xspacing):xspacing:((nvals(iY)/2)*xspacing))-xspacing/2;
    elseif nvals == 1
        x(ind) = 0;
    end
    ind = ind + nvals(iY);
end