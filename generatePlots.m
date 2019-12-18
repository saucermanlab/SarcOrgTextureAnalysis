%% generatePlots
% Generates the plots used in the figures.
close all;
eCM_color = [241,163,64]./255;
rCM_color = [153,142,195]./255;

xspacing = 0.05;
xwidth = 0.2;
nBins = 64;

figureSize = [300 400];
gcaPosition = [0.2 0.1 0.75 0.8];
scattergcaPosition = [0.175 0.2 0.75 0.7];
sarcomereOrganizationThreshold = 0.1;
elongationThreshold = 1;
markerSize = 12;
markerSizeScatter = 24;
alphaValue = 1;
fontSize = 12;
fitLineColor = [0.1 0.1 0.1 1];
belowThreshScatterColor = [0.7 0.7 0.7];
xticklabels = {'CM','iCLM'};
scatterFigureSize = [300 300];
scatterFontSize = 12;
xlimits = [0.4 2.6];
ticklengths = [0.02 0.025];
lw = 1;

%% Cell Area
ylimits = [0 16000];
[eX,eY] = beeswarmPoints([eCM.all.Area_um2],xspacing,(ylimits(2)-ylimits(1))/nBins,ylimits);
[rX,rY] = beeswarmPoints([rCM.all.Area_um2],xspacing,(ylimits(2)-ylimits(1))/nBins,ylimits);
figure('units','pixels','position',[50 50 figureSize]); hold on;
scatter(eX+1,eY,markerSize,'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',eCM_color,'MarkerFaceAlpha',alphaValue);
plot3Lines([eCM.all.Area_um2],1,xwidth,lw);
scatter(rX+2,rY,markerSize,'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',rCM_color,'MarkerFaceAlpha',alphaValue);
plot3Lines([rCM.all.Area_um2],2,xwidth,lw);
set(gca,'FontSize',fontSize,'LineWidth',lw,'Color','None','TickDir','out','TickLength',ticklengths);
ax = gca;
ax.XAxis.Color = 'k';
ax.YAxis.Color = 'k';
xlim(xlimits);
ylim(ylimits);
ylabel('Cell area (µm²)');
set(gca,'XTick',[1 2],'XTickLabel',xticklabels);
set(gca,'Position',gcaPosition);
set(gca,'YTick',0:2000:16000);
ax = gca;
ax.YAxis.Exponent = 3;
%% Cell Elongation
ylimits = [1 6];
[eX,eY] = beeswarmPoints([eCM.all.Elongation],xspacing,(ylimits(2)-ylimits(1))/nBins,ylimits);
[rX,rY] = beeswarmPoints([rCM.all.Elongation],xspacing,(ylimits(2)-ylimits(1))/nBins,ylimits);
figure('units','pixels','position',[350 50 figureSize]); hold on;
scatter(eX+1,eY,markerSize,'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',eCM_color,'MarkerFaceAlpha',alphaValue);
plot3Lines([eCM.all.Elongation],1,xwidth,lw);
scatter(rX+2,rY,markerSize,'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',rCM_color,'MarkerFaceAlpha',alphaValue);
plot3Lines([rCM.all.Elongation],2,xwidth,lw);
set(gca,'FontSize',fontSize,'LineWidth',lw,'Color','None','TickDir','out','TickLength',ticklengths);
ax = gca;
ax.XAxis.Color = 'k';
ax.YAxis.Color = 'k';
xlim(xlimits);
ylim(ylimits);
ylabel('Cell elongation');
set(gca,'XTick',[1 2],'XTickLabel',xticklabels);
set(gca,'Position',gcaPosition);
%% Cell Eccentricity
ylimits = [0 1];
[eX,eY] = beeswarmPoints([eCM.all.Eccentricity],xspacing,(ylimits(2)-ylimits(1))/nBins,ylimits);
[rX,rY] = beeswarmPoints([rCM.all.Eccentricity],xspacing,(ylimits(2)-ylimits(1))/nBins,ylimits);
figure('units','pixels','position',[650 50 figureSize]); hold on;
scatter(eX+1,eY,markerSize,'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',eCM_color,'MarkerFaceAlpha',alphaValue);
plot3Lines([eCM.all.Eccentricity],1,xwidth,lw);
scatter(rX+2,rY,markerSize,'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',rCM_color,'MarkerFaceAlpha',alphaValue);
plot3Lines([rCM.all.Eccentricity],2,xwidth,lw);
set(gca,'FontSize',fontSize,'LineWidth',lw,'Color','None','TickDir','out','TickLength',ticklengths);
ax = gca;
ax.XAxis.Color = 'k';
ax.YAxis.Color = 'k';
xlim(xlimits);
ylim(ylimits);
ylabel('Cell eccentricity');
set(gca,'XTick',[1 2],'XTickLabel',xticklabels);
set(gca,'Position',gcaPosition);
%% Cell Circularity
ylimits = [0 1];
[eX,eY] = beeswarmPoints([eCM.all.Circularity],xspacing,(ylimits(2)-ylimits(1))/nBins,ylimits);
[rX,rY] = beeswarmPoints([rCM.all.Circularity],xspacing,(ylimits(2)-ylimits(1))/nBins,ylimits);
figure('units','pixels','position',[950 50 figureSize]); hold on;
scatter(eX+1,eY,markerSize,'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',eCM_color,'MarkerFaceAlpha',alphaValue);
plot3Lines([eCM.all.Circularity],1,xwidth,lw);
scatter(rX+2,rY,markerSize,'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',rCM_color,'MarkerFaceAlpha',alphaValue);
plot3Lines([rCM.all.Circularity],2,xwidth,lw);
set(gca,'FontSize',fontSize,'LineWidth',lw,'Color','None','TickDir','out','TickLength',ticklengths);
ax = gca;
ax.XAxis.Color = 'k';
ax.YAxis.Color = 'k';
xlim(xlimits);
ylim(ylimits);
ylabel('Cell circularity');
set(gca,'XTick',[1 2],'XTickLabel',xticklabels);
set(gca,'Position',gcaPosition);

%% Sarcomere Organization
ylimits = [0 0.6];
xlimits = [0.25 2.75];
[eX,eY] = beeswarmPoints([eCM.all.SarcomereOrganizationScore],xspacing,(ylimits(2)-ylimits(1))/nBins,ylimits);
[rX,rY] = beeswarmPoints([rCM.all_aact.SarcomereOrganizationScore],xspacing,(ylimits(2)-ylimits(1))/nBins,ylimits);
figure('units','pixels','position',[50 500 figureSize]); hold on;
plot(xlimits,sarcomereOrganizationThreshold.*[1 1],'k--','LineWidth',lw);
scatter(eX+1,eY,markerSize,'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',eCM_color,'MarkerFaceAlpha',alphaValue);
plot3Lines([eCM.all.SarcomereOrganizationScore],1,xwidth,lw);
scatter(rX+2,rY,markerSize,'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',rCM_color,'MarkerFaceAlpha',alphaValue);
plot3Lines([rCM.all_aact.SarcomereOrganizationScore],2,xwidth,lw);

% plot([0.9 1.1],mean([eCM.all.SarcomereOrganizationScore])*[1 1],'k');
% plot([1.9 2.1],mean([rCM.all_aact.SarcomereOrganizationScore])*[1 1],'k');

set(gca,'FontSize',fontSize,'LineWidth',lw,'Color','None','TickDir','out','TickLength',[0.02 0.025]);
xlim(xlimits);
ylim(ylimits);
set(gca,'YTick',0:0.1:0.6);
ylabel('Sarcomere organization');
set(gca,'XTick',[1 2],'XTickLabel',xticklabels);
set(gca,'Position',gcaPosition);
ax = gca;
ax.XAxis.Color = 'k';
ax.YAxis.Color = 'k';
%% Sarcomere Length
eCM_sarcOrg = [eCM.all.SarcomereOrganizationScore];
rCM_sarcOrg = [rCM.all_aact.SarcomereOrganizationScore];

eCM_sarcLength = [eCM.all.SarcomereLength_um];
rCM_sarcLength = [rCM.all_aact.SarcomereLength_um];

eCM_sarcLength = eCM_sarcLength(eCM_sarcOrg >= sarcomereOrganizationThreshold);
rCM_sarcLength = rCM_sarcLength(rCM_sarcOrg >= sarcomereOrganizationThreshold);

ylimits = [0 3];
[eX,eY] = beeswarmPoints(eCM_sarcLength,xspacing,(ylimits(2)-ylimits(1))/nBins,ylimits);
[rX,rY] = beeswarmPoints(rCM_sarcLength,xspacing,(ylimits(2)-ylimits(1))/nBins,ylimits);
figure('units','pixels','position',[350 500 figureSize]); hold on;
scatter(eX+1,eY,markerSize,'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',eCM_color,'MarkerFaceAlpha',alphaValue);
plot3Lines(eCM_sarcLength,1,xwidth,lw);
scatter(rX+2,rY,markerSize,'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',rCM_color,'MarkerFaceAlpha',alphaValue);
plot3Lines(rCM_sarcLength,2,xwidth,lw);
set(gca,'FontSize',fontSize,'LineWidth',lw,'Color','None','TickDir','out','TickLength',[0.02 0.025]);
xlim([0.25 2.75]);
ylim(ylimits);
ylabel('Sarcomere length (µm)');
set(gca,'YTick',0:0.5:3);
set(gca,'XTick',[1 2],'XTickLabel',xticklabels);
set(gca,'Position',gcaPosition);
ax = gca;
ax.XAxis.Color = 'k';
ax.YAxis.Color = 'k';
%% Cell-Sarcomere Misalignment
eCM_sarcOrg = [eCM.all.SarcomereOrganizationScore];
rCM_sarcOrg = [rCM.all_aact.SarcomereOrganizationScore];

eCM_elong = [eCM.all.Elongation];
rCM_elong = [rCM.all_aact.Elongation];

eCM_sarcMisalignment = [eCM.all.CellSarcomereAlignment];
rCM_sarcMisalignment = [rCM.all_aact.CellSarcomereAlignment];

eCM_sarcMisalignment = eCM_sarcMisalignment(eCM_sarcOrg > sarcomereOrganizationThreshold & eCM_elong > elongationThreshold);
rCM_sarcMisalignment = rCM_sarcMisalignment(rCM_sarcOrg > sarcomereOrganizationThreshold & rCM_elong > elongationThreshold);

ylimits = [0 90];
[eX,eY] = beeswarmPoints(eCM_sarcMisalignment,xspacing,(ylimits(2)-ylimits(1))/nBins,ylimits);
[rX,rY] = beeswarmPoints(rCM_sarcMisalignment,xspacing,(ylimits(2)-ylimits(1))/nBins,ylimits);
figure('units','pixels','position',[650 500 figureSize]); hold on;
scatter(eX+1,eY,markerSize,'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',eCM_color,'MarkerFaceAlpha',alphaValue);
plot3Lines(eCM_sarcMisalignment,1,xwidth,lw);
scatter(rX+2,rY,markerSize,'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',rCM_color,'MarkerFaceAlpha',alphaValue);
plot3Lines(rCM_sarcMisalignment,2,xwidth,lw);
set(gca,'FontSize',fontSize,'LineWidth',lw,'Color','None','TickDir','out','TickLength',[0.02 0.025]);
xlim([0.25 2.75]);
ylim(ylimits);
ylabel('Cell-sarcomere misalignment');
set(gca,'XTick',[1 2],'XTickLabel',xticklabels);
set(gca,'YTick',0:10:90,'YTickLabel',{'0°','10°','20°','30°','40°','50°','60°','70°','80°','90°'});
set(gca,'Position',gcaPosition);
ax = gca;
ax.XAxis.Color = 'k';
ax.YAxis.Color = 'k';

%% Cell Elongation vs Cell Area
xlimits = [0 16000];
ylimits = [1 6];

% Endogenous CM
figure('units','pixels','position',[50 50 scatterFigureSize]); hold on;
scatter([eCM.all.Area_um2],[eCM.all.Elongation],markerSizeScatter,'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',eCM_color,'MarkerFaceAlpha',alphaValue);
set(gca,'FontSize',scatterFontSize,'LineWidth',1,'Color','None','TickDir','out','TickLength',[0.02 0.025]);
xlim(xlimits);
ylim(ylimits);
xlabel('Cell area (µm²)');
ylabel('Cell elongation');
set(gca,'Position',scattergcaPosition);
ax = gca;
ax.XAxis.Color = 'k';
ax.YAxis.Color = 'k';
ax.XAxis.Exponent = 3;
set(gca,'XTick',0:4000:16000);

X = [ones(length([eCM.all.Area_um2]),1),[eCM.all.Area_um2]'];
k = X\[eCM.all.Elongation]';
xfit = [min([eCM.all.Area_um2]) max([eCM.all.Area_um2])];
yfit = k(1)+k(2).*xfit;
plot(xfit,yfit,'LineWidth',2,'Color',fitLineColor);

P = fitlm([eCM.all.Area_um2],[eCM.all.Elongation],'linear');

text(xlimits(2),ylimits(2),['R²=' sprintf('%.2f',P.Rsquared.Ordinary)],'FontSize',12,'FontName','Arial','HorizontalAlignment','Right','VerticalAlignment','Top');

% Reprogrammed CM
figure('units','pixels','position',[50 50 scatterFigureSize]); hold on;
scatter([rCM.all.Area_um2],[rCM.all.Elongation],markerSizeScatter,'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',rCM_color,'MarkerFaceAlpha',alphaValue);
set(gca,'FontSize',scatterFontSize,'LineWidth',1,'Color','None','TickDir','out','TickLength',[0.02 0.025]);
xlim(xlimits);
ylim(ylimits);
xlabel('Cell area (µm²)');
ylabel('Cell elongation');
set(gca,'Position',scattergcaPosition);
ax = gca;
ax.XAxis.Color = 'k';
ax.YAxis.Color = 'k';
ax.XAxis.Exponent = 3;
set(gca,'XTick',0:4000:16000);

X = [ones(length([rCM.all.Area_um2]),1),[rCM.all.Area_um2]'];
k = X\[rCM.all.Elongation]';
xfit = [min([rCM.all.Area_um2]) max([rCM.all.Area_um2])];
yfit = k(1)+k(2).*xfit;
plot(xfit,yfit,'LineWidth',2,'Color',fitLineColor);

P = fitlm([rCM.all.Area_um2],[rCM.all.Elongation],'linear');

text(xlimits(2),ylimits(2),['R²=' sprintf('%.2f',P.Rsquared.Ordinary)],'FontSize',12,'FontName','Arial','HorizontalAlignment','Right','VerticalAlignment','Top');
%% Cell Area vs Sarcomere Organization
xlimits = [0 0.6];
ylimits = [0 16000];

% Endogenous CM
figure('units','pixels','position',[50 50 scatterFigureSize]); hold on;
scatter([eCM.all.SarcomereOrganizationScore],[eCM.all.Area_um2],markerSizeScatter,'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',eCM_color,'MarkerFaceAlpha',alphaValue);
set(gca,'FontSize',scatterFontSize,'LineWidth',lw,'Color','None','TickDir','out','TickLength',[0.02 0.025]);
xlim(xlimits);
ylim(ylimits);
xlabel('Sarcomere organization');
ylabel('Cell area (µm²)');
set(gca,'Position',scattergcaPosition);
ax = gca;
ax.XAxis.Color = 'k';
ax.YAxis.Color = 'k';
ax.YAxis.Exponent = 3;
set(gca,'YTick',0:4000:16000);

X = [ones(length([eCM.all.SarcomereOrganizationScore]),1),[eCM.all.SarcomereOrganizationScore]'];
k = X\[eCM.all.Area_um2]';
xfit = [min([eCM.all.SarcomereOrganizationScore]) max([eCM.all.SarcomereOrganizationScore])];
yfit = k(1)+k(2).*xfit;
plot(xfit,yfit,'LineWidth',2,'Color',fitLineColor);

P = fitlm([eCM.all.SarcomereOrganizationScore],[eCM.all.Area_um2],'linear');

text(xlimits(2),ylimits(2),['R²=' sprintf('%.2f',P.Rsquared.Ordinary)],'FontSize',scatterFontSize,'FontName','Arial','HorizontalAlignment','Right','VerticalAlignment','Top');

% Reprogrammed CM
figure('units','pixels','position',[50 50 scatterFigureSize]); hold on;
scatter([rCM.all.SarcomereOrganizationScore],[rCM.all.Area_um2],markerSizeScatter,'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',rCM_color,'MarkerFaceAlpha',alphaValue);
set(gca,'FontSize',scatterFontSize,'LineWidth',lw,'Color','None','TickDir','out','TickLength',[0.02 0.025]);
xlim(xlimits);
ylim(ylimits);
xlabel('Sarcomere organization');
ylabel('Cell area (µm²)');
set(gca,'Position',scattergcaPosition);
ax = gca;
ax.XAxis.Color = 'k';
ax.YAxis.Color = 'k';
ax.YAxis.Exponent = 3;
set(gca,'YTick',0:4000:16000);

X = [ones(length([rCM.all.SarcomereOrganizationScore]),1),[rCM.all.SarcomereOrganizationScore]'];
k = X\[rCM.all.Area_um2]';
xfit = [min([rCM.all.SarcomereOrganizationScore]) max([rCM.all.SarcomereOrganizationScore])];
yfit = k(1)+k(2).*xfit;
plot(xfit,yfit,'LineWidth',2,'Color',fitLineColor);

P = fitlm([rCM.all.SarcomereOrganizationScore],[rCM.all.Area_um2],'linear');

text(xlimits(2),ylimits(2),['R²=' sprintf('%.2f',P.Rsquared.Ordinary)],'FontSize',scatterFontSize,'FontName','Arial','HorizontalAlignment','Right','VerticalAlignment','Top');
%% Cell Elongation vs Sarcomere Organization
xlimits = [0 0.6];
ylimits = [1 6];

% Endogenous CM
figure('units','pixels','position',[50 50 scatterFigureSize]); hold on;
scatter([eCM.all.SarcomereOrganizationScore],[eCM.all.Elongation],markerSizeScatter,'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',eCM_color,'MarkerFaceAlpha',alphaValue);
set(gca,'FontSize',scatterFontSize,'LineWidth',lw,'Color','None','TickDir','out','TickLength',[0.02 0.025]);
xlim(xlimits);
ylim(ylimits);
ax = gca;
ax.XAxis.Color = 'k';
ax.YAxis.Color = 'k';
xlabel('Sarcomere organization');
ylabel('Cell elongation');
set(gca,'Position',scattergcaPosition);

X = [ones(length([eCM.all.SarcomereOrganizationScore]),1),[eCM.all.SarcomereOrganizationScore]'];
k = X\[eCM.all.Elongation]';
xfit = [min([eCM.all.SarcomereOrganizationScore]) max([eCM.all.SarcomereOrganizationScore])];
yfit = k(1)+k(2).*xfit;
plot(xfit,yfit,'LineWidth',2,'Color',fitLineColor);

P = fitlm([eCM.all.SarcomereOrganizationScore],[eCM.all.Elongation],'linear');

text(xlimits(2),ylimits(2),['R²=' sprintf('%.2f',P.Rsquared.Ordinary)],'FontSize',scatterFontSize,'FontName','Arial','HorizontalAlignment','Right','VerticalAlignment','Top');

% Reprogrammed CM
figure('units','pixels','position',[50 50 scatterFigureSize]); hold on;
scatter([rCM.all.SarcomereOrganizationScore],[rCM.all.Elongation],markerSizeScatter,'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',rCM_color,'MarkerFaceAlpha',alphaValue);
set(gca,'FontSize',scatterFontSize,'LineWidth',lw,'Color','None','TickDir','out','TickLength',[0.02 0.025]);
xlim(xlimits);
ylim(ylimits);
xlabel('Sarcomere organization');
ylabel('Cell elongation');
set(gca,'Position',scattergcaPosition);
ax = gca;
ax.XAxis.Color = 'k';
ax.YAxis.Color = 'k';
X = [ones(length([rCM.all.SarcomereOrganizationScore]),1),[rCM.all.SarcomereOrganizationScore]'];
k = X\[rCM.all.Elongation]';
xfit = [min([rCM.all.SarcomereOrganizationScore]) max([rCM.all.SarcomereOrganizationScore])];
yfit = k(1)+k(2).*xfit;
plot(xfit,yfit,'LineWidth',2,'Color',fitLineColor);

P = fitlm([rCM.all.SarcomereOrganizationScore],[rCM.all.Elongation],'linear');

text(xlimits(2),ylimits(2),['R²=' sprintf('%.2f',P.Rsquared.Ordinary)],'FontSize',scatterFontSize,'FontName','Arial','HorizontalAlignment','Right','VerticalAlignment','Top');
%% Sarcomere Length vs Sarcomere Organization
xlimits = [0 0.6];
ylimits = [0 6];

eCM_sarcOrg = [eCM.all.SarcomereOrganizationScore];
rCM_sarcOrg = [rCM.all_aact.SarcomereOrganizationScore];

eCM_sarcLength = [eCM.all.SarcomereLength_um];
rCM_sarcLength = [rCM.all_aact.SarcomereLength_um];

eCM_sarcLength_org = eCM_sarcLength(eCM_sarcOrg >= sarcomereOrganizationThreshold);
rCM_sarcLength_org = rCM_sarcLength(rCM_sarcOrg >= sarcomereOrganizationThreshold);

eCM_sarcLength_disorg = eCM_sarcLength(eCM_sarcOrg < sarcomereOrganizationThreshold);
rCM_sarcLength_disorg = rCM_sarcLength(rCM_sarcOrg < sarcomereOrganizationThreshold);

eCM_sarcOrg_org = eCM_sarcOrg(eCM_sarcOrg >= sarcomereOrganizationThreshold);
rCM_sarcOrg_org = rCM_sarcOrg(rCM_sarcOrg >= sarcomereOrganizationThreshold);

eCM_sarcOrg_disorg = eCM_sarcOrg(eCM_sarcOrg < sarcomereOrganizationThreshold);
rCM_sarcOrg_disorg = rCM_sarcOrg(rCM_sarcOrg < sarcomereOrganizationThreshold);

% Endogenous CM
figure('units','pixels','position',[50 50 scatterFigureSize]); hold on;
scatter(eCM_sarcOrg_disorg,eCM_sarcLength_disorg,markerSizeScatter*2/3,'MarkerEdgeColor',eCM_color,'MarkerFaceColor','none','MarkerFaceAlpha',alphaValue,'LineWidth',1);
scatter(eCM_sarcOrg_org,eCM_sarcLength_org,markerSizeScatter,'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',eCM_color,'MarkerFaceAlpha',alphaValue);
area([-10 sarcomereOrganizationThreshold],[10 10],-10,'FaceColor',[0.9 0.9 0.9],'FaceAlpha',1,'LineStyle','--','LineWidth',1);
set(gca,'FontSize',scatterFontSize,'LineWidth',lw,'Color','None','TickDir','out','TickLength',[0.02 0.025]);
xlim(xlimits);
ylim(ylimits);
xlabel('Sarcomere organization');
ylabel('Sarcomere length (µm)');
set(gca,'Position',scattergcaPosition);
ax = gca;
ax.XAxis.Color = 'k';
ax.YAxis.Color = 'k';
X = [ones(length(eCM_sarcOrg_org),1),eCM_sarcOrg_org'];
k = X\eCM_sarcLength_org';
xfit = [min(eCM_sarcOrg_org) max(eCM_sarcOrg_org)];
yfit = k(1)+k(2).*xfit;
plot(xfit,yfit,'LineWidth',2,'Color',fitLineColor);

P = fitlm(eCM_sarcOrg_org,eCM_sarcLength_org,'linear');

text(xlimits(2),ylimits(2),['R²=' sprintf('%.2f',P.Rsquared.Ordinary)],'FontSize',scatterFontSize,'FontName','Arial','HorizontalAlignment','Right','VerticalAlignment','Top');

% Reprogrammed CM
figure('units','pixels','position',[50 50 scatterFigureSize]); hold on;
scatter(rCM_sarcOrg_disorg,rCM_sarcLength_disorg,markerSizeScatter*2/3,'MarkerEdgeColor',rCM_color,'MarkerFaceColor','none','MarkerFaceAlpha',alphaValue,'LineWidth',1);
scatter(rCM_sarcOrg_org,rCM_sarcLength_org,markerSizeScatter,'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',rCM_color,'MarkerFaceAlpha',alphaValue);
area([-10 sarcomereOrganizationThreshold],[10 10],-10,'FaceColor',[0.9 0.9 0.9],'FaceAlpha',1,'LineStyle','--','LineWidth',1);
set(gca,'FontSize',scatterFontSize,'LineWidth',lw,'Color','None','TickDir','out','TickLength',[0.02 0.025]);
xlim(xlimits);
ylim(ylimits);
xlabel('Sarcomere organization');
ylabel('Sarcomere length (µm)');
set(gca,'Position',scattergcaPosition);
ax = gca;
ax.XAxis.Color = 'k';
ax.YAxis.Color = 'k';
X = [ones(length(rCM_sarcOrg_org),1),rCM_sarcOrg_org'];
k = X\rCM_sarcLength_org';
xfit = [min(rCM_sarcOrg_org) max(rCM_sarcOrg_org)];
yfit = k(1)+k(2).*xfit;
plot(xfit,yfit,'LineWidth',2,'Color',fitLineColor);

P = fitlm(rCM_sarcOrg_org,rCM_sarcLength_org,'linear');

text(xlimits(2),ylimits(2),['R²=' sprintf('%.2f',P.Rsquared.Ordinary)],'FontSize',scatterFontSize,'FontName','Arial','HorizontalAlignment','Right','VerticalAlignment','Top');
%% Elongation vs Cell-Sarcomere Misalignment
xlimits = [0 90];
ylimits = [1 6];

eCM_sarcOrg = [eCM.all.SarcomereOrganizationScore];
rCM_sarcOrg = [rCM.all_aact.SarcomereOrganizationScore];

eCM_Elongation = [eCM.all.Elongation];
rCM_Elongation = [rCM.all_aact.Elongation];

eCM_CellSarcomereAlignment = [eCM.all.CellSarcomereAlignment];
rCM_CellSarcomereAlignment = [rCM.all_aact.CellSarcomereAlignment];

eCM_Elongation_org = eCM_Elongation(eCM_sarcOrg >= sarcomereOrganizationThreshold);
rCM_Elongation_org = rCM_Elongation(rCM_sarcOrg >= sarcomereOrganizationThreshold);

eCM_Elongation_disorg = eCM_Elongation(eCM_sarcOrg < sarcomereOrganizationThreshold);
rCM_Elongation_disorg = rCM_Elongation(rCM_sarcOrg < sarcomereOrganizationThreshold);

eCM_CellSarcomereAlignment_org = eCM_CellSarcomereAlignment(eCM_sarcOrg >= sarcomereOrganizationThreshold);
rCM_CellSarcomereAlignment_org = rCM_CellSarcomereAlignment(rCM_sarcOrg >= sarcomereOrganizationThreshold);

eCM_CellSarcomereAlignment_disorg = eCM_CellSarcomereAlignment(eCM_sarcOrg < sarcomereOrganizationThreshold);
rCM_CellSarcomereAlignment_disorg = rCM_CellSarcomereAlignment(rCM_sarcOrg < sarcomereOrganizationThreshold);

% Endogenous CM
figure('units','pixels','position',[50 50 scatterFigureSize]); hold on;
scatter(eCM_CellSarcomereAlignment_disorg,eCM_Elongation_disorg,markerSizeScatter*2/3,'MarkerEdgeColor',eCM_color,'MarkerFaceColor','none','MarkerFaceAlpha',alphaValue,'LineWidth',1);
scatter(eCM_CellSarcomereAlignment_org,eCM_Elongation_org,markerSizeScatter,'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',eCM_color,'MarkerFaceAlpha',alphaValue);

set(gca,'FontSize',scatterFontSize,'LineWidth',lw,'Color','None','TickDir','out','TickLength',[0.02 0.025]);
xlim(xlimits);
ylim(ylimits);
xlabel('Cell-sarcomere misalignment');
ylabel('Cell elongation');
set(gca,'Position',scattergcaPosition);
set(gca,'XTick',0:30:90,'XTickLabel',{'0°','30°','60°','90°'});
ax = gca;
ax.XAxis.Color = 'k';
ax.YAxis.Color = 'k';
X = [ones(length(eCM_CellSarcomereAlignment_org),1),eCM_CellSarcomereAlignment_org'];
k = X\eCM_Elongation_org';
xfit = [min(eCM_CellSarcomereAlignment_org) max(eCM_CellSarcomereAlignment_org)];
yfit = k(1)+k(2).*xfit;
plot(xfit,yfit,'LineWidth',2,'Color',fitLineColor);

P = fitlm(eCM_CellSarcomereAlignment_org,eCM_Elongation_org,'linear');

text(xlimits(2),ylimits(2),['R²=' sprintf('%.2f',P.Rsquared.Ordinary)],'FontSize',scatterFontSize,'FontName','Arial','HorizontalAlignment','Right','VerticalAlignment','Top');

% Reprogrammed CM
figure('units','pixels','position',[50 50 scatterFigureSize]); hold on;
scatter(rCM_CellSarcomereAlignment_disorg,rCM_Elongation_disorg,markerSizeScatter*2/3,'MarkerEdgeColor',rCM_color,'MarkerFaceColor','none','MarkerFaceAlpha',alphaValue,'LineWidth',1);
scatter(rCM_CellSarcomereAlignment_org,rCM_Elongation_org,markerSizeScatter,'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',rCM_color,'MarkerFaceAlpha',alphaValue);

set(gca,'FontSize',scatterFontSize,'LineWidth',lw,'Color','None','TickDir','out','TickLength',[0.02 0.025]);
xlim(xlimits);
ylim(ylimits);
xlabel('Cell-sarcomere misalignment');
ylabel('Cell elongation');
set(gca,'Position',scattergcaPosition);
set(gca,'XTick',0:30:90,'XTickLabel',{'0°','30°','60°','90°'});
ax = gca;
ax.XAxis.Color = 'k';
ax.YAxis.Color = 'k';
X = [ones(length(rCM_CellSarcomereAlignment_org),1),rCM_CellSarcomereAlignment_org'];
k = X\rCM_Elongation_org';
xfit = [min(rCM_CellSarcomereAlignment_org) max(rCM_CellSarcomereAlignment_org)];
yfit = k(1)+k(2).*xfit;
plot(xfit,yfit,'LineWidth',2,'Color',fitLineColor);

P = fitlm(rCM_CellSarcomereAlignment_org,rCM_Elongation_org,'linear');

text(xlimits(2),ylimits(2),['R²=' sprintf('%.2f',P.Rsquared.Ordinary)],'FontSize',scatterFontSize,'FontName','Arial','HorizontalAlignment','Right','VerticalAlignment','Top');