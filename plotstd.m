function plotstd(x,mean,std,width,linewidth,bottom)
%% FUNCTION DESCRIPTION
% This function plots the standard deviation lines used in Figure 5.
%
% INPUTS
% x : horizontal location centered x position where the bar is plotted
% mean : mean of data points
% std : standard deviation of data points
% width : horizontal width of the error bars
% linewidth : thickness of the lines
% bottom : 1 if plotting both the upper and lower error bars, 0 if only
%   upper

%%
if bottom
    plot([x-width x+width],[1 1].*(mean-std),'k','LineWidth',linewidth);
    plot([1 1].*x,[mean-std mean+std],'k','LineWidth',linewidth);
else
    plot([1 1].*x,[mean mean+std],'k','LineWidth',linewidth);
end
plot([x-width x+width],[1 1].*(mean+std),'k','LineWidth',linewidth);


end