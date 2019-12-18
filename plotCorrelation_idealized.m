function f = plotCorrelation_idealized(dataStruct,annotate,ang)
%% FUNCTION DESCRIPTION
% This function plots the correlation vs. pixel offset distance, as shown
% in Figure 3.
%
% INPUTS
% dataStruct : metric structure for a single idealized image with
%   correlation data
% d1_conversion : pixels to microns conversion for a single dimension
% ang : angle to plot in black.

% OUTPUTS
% f : function handle

%%
data = dataStruct.CorrelationData;
f = figure('units','pixels','position',[20 -10 450 300]);
hold on;
cmap = hsv(size(data,1));
x = 0:1:(size(data,2)-1);
for i = 1:size(data,1)
    plot(x,data(i,:),'Color',[cmap(i,:) 1]);
end
if ang
    if dataStruct.SarcomereOrientation > 0
        plot(x,data(dataStruct.SarcomereOrientationIndex,:),'Color',[0 0 0],'LineWidth',2);
    else
        plot(x,data(dataStruct.SarcomereOrientationIndex,:),'Color',[0 0 0],'LineWidth',2);
    end
end
if annotate
    [~,minLocs] = findpeaks(-data(dataStruct.SarcomereOrientation+1,:));
    [~,maxLocs] = findpeaks(data(dataStruct.SarcomereOrientation+1,:));
    
    plot([maxLocs(1).*d1_conversion.*1.15./10 minLocs(1).*d1_conversion./10],data(dataStruct.SarcomereOrientation+1,minLocs(1)).*[1 1],'k:','LineWidth',3);
    plot([maxLocs(1).*d1_conversion.*1.15./10 maxLocs(1).*d1_conversion./10],data(dataStruct.SarcomereOrientation+1,maxLocs(1)).*[1 1],'k:','LineWidth',3);
    
    plot([maxLocs(1).*d1_conversion.*1.15./10 maxLocs(1).*d1_conversion.*1.25./10],data(dataStruct.SarcomereOrientation+1,minLocs(1)).*[1 1],'k','LineWidth',3);
    plot([maxLocs(1).*d1_conversion.*1.15./10 maxLocs(1).*d1_conversion.*1.25./10],data(dataStruct.SarcomereOrientation+1,maxLocs(1)).*[1 1],'k','LineWidth',3);
    plot(maxLocs(1).*d1_conversion./10.*[1 1],[data(dataStruct.SarcomereOrientationIndex,maxLocs(1)).*1.05 data(dataStruct.SarcomereOrientationIndex,maxLocs(1)).*1.2],'k','LineWidth',3);
end
set(gca,'FontName','Arial');
set(gca,'FontSize',12);
set(gca,'Linewidth',1);
xlabel('Offset distance (pixels)');
ylabel('Haralick correlation');
set(gca,'TickDir','Out');
set(gca,'XColor','k','YColor','k');
set(gca,'Color','None');
xlim([0 (size(data,2)-1)]);
ylim([-1 1]);
set(gca,'TickLength',[0.015 0.025]);
if ang
    if dataStruct.SarcomereOrientation > 0
        cmap(dataStruct.SarcomereOrientationIndex-1,:) = [0 0 0];
        cmap(dataStruct.SarcomereOrientationIndex,:) = [0 0 0];
    else
        cmap(dataStruct.SarcomereOrientationIndex+1,:) = [0 0 0];
    end
    if dataStruct.SarcomereOrientationIndex <= size(cmap,1)
        cmap(dataStruct.SarcomereOrientationIndex+1,:) = [0 0 0];
    end
end
colormap(cmap);
b = colorbar('Ticks',[0 0.25 0.5 0.75 1],'TickLabels',{'0°','45°','90°','135°','180°'});
b.LineWidth = 1;
b.TickLength = 0.015;
b.Color = 'k';
b.FontSize = 12;
% b.Position = [b.Position(1:2) 0.02 b.Position(4)];
end