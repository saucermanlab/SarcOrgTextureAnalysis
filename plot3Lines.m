function plot3Lines(yvals,xpos,xwidth,linewidth)
%% FUNCTION DESCRIPTION
% This function plots the median, 25th and 75th percentiles on top of the
% beeswarm plots.
%
% INPUTS
% yvals : array with data points
% xpos : centered x position where the data is plotted
% xwidth : width of median line
% linewidth : thickness of the lines

%%
prctilevals = prctile(yvals,[25 50 75]);
plot([xpos-xwidth xpos+xwidth],prctilevals(2).*[1 1],'Color',[0 0 0],'LineWidth',linewidth);
plot([xpos-(xwidth/2) xpos+(xwidth/2)],prctilevals(1).*[1 1],'Color',[0 0 0],'LineWidth',linewidth);
plot([xpos-(xwidth/2) xpos+(xwidth/2)],prctilevals(3).*[1 1],'Color',[0 0 0],'LineWidth',linewidth);
plot(xpos.*[1 1],[prctilevals(1) prctilevals(3)],'Color',[0 0 0],'LineWidth',linewidth);
end
