%% generatePlots_subtype
% Generates the plots used in the supplemental figures (where data is shown by subtype). 
close all;
eCM_color = [241,163,64]./255;
rCM_color = [153,142,195]./255;

edgecolor = [0 0 0];
xspacing = 0.085;
xwidth = 0.2;
nBins = 64;

figureSize = [50 100 566 400];
figureSize2 = [50 100 465 400];
gcaPosition = [0.15 0.1 0.8 0.825];
scattergcaPosition = [0.175 0.2 0.75 0.7];
sarcomereOrganizationThreshold = 0.1;
elongationThreshold = 1;
markerSize = 12;
markerSizeScatter = 96;
alphaValue = 1;
fontSize = 13;
fitLineColor = [0.1 0.1 0.1 1];
belowThreshScatterColor = [0.7 0.7 0.7];
xticklabels = {};
scatterFontSize = 32;

%% Cell Area
ylimits = [0 16000];
[e110X,e110Y] = beeswarmPoints([eCM.aacthcn4.Area_um2],xspacing,(ylimits(2)-ylimits(1))/nBins,ylimits);
[e101X,e101Y] = beeswarmPoints([eCM.aactnppa.Area_um2],xspacing,(ylimits(2)-ylimits(1))/nBins,ylimits);
[r100X,r100Y] = beeswarmPoints([rCM.aact.Area_um2],xspacing,(ylimits(2)-ylimits(1))/nBins,ylimits);
[r110X,r110Y] = beeswarmPoints([rCM.aacthcn4.Area_um2],xspacing,(ylimits(2)-ylimits(1))/nBins,ylimits);
[r101X,r101Y] = beeswarmPoints([rCM.aactnppa.Area_um2],xspacing,(ylimits(2)-ylimits(1))/nBins,ylimits);
[r010X,r010Y] = beeswarmPoints([rCM.hcn4.Area_um2],xspacing,(ylimits(2)-ylimits(1))/nBins,ylimits);
% r011X = []; r011Y = [];
[r111X,r111Y] = beeswarmPoints([rCM.aacthcn4nppa.Area_um2],xspacing,(ylimits(2)-ylimits(1))/nBins,ylimits);
figure('units','pixels','position',figureSize); hold on;
scatter(e110X+1,e110Y,markerSize,'MarkerEdgeColor',edgecolor,'MarkerFaceColor',eCM_color,'MarkerFaceAlpha',alphaValue,'LineWidth',0.25);
scatter(e101X+2,e101Y,markerSize,'MarkerEdgeColor',edgecolor,'MarkerFaceColor',eCM_color,'MarkerFaceAlpha',alphaValue,'LineWidth',0.25);
scatter(r100X+3,r100Y,markerSize,'MarkerEdgeColor',edgecolor,'MarkerFaceColor',rCM_color,'MarkerFaceAlpha',alphaValue,'LineWidth',0.25);
scatter(r110X+4,r110Y,markerSize,'MarkerEdgeColor',edgecolor,'MarkerFaceColor',rCM_color,'MarkerFaceAlpha',alphaValue,'LineWidth',0.25);
scatter(r101X+5,r101Y,markerSize,'MarkerEdgeColor',edgecolor,'MarkerFaceColor',rCM_color,'MarkerFaceAlpha',alphaValue,'LineWidth',0.25);
scatter(r010X+6,r010Y,markerSize,'MarkerEdgeColor',edgecolor,'MarkerFaceColor',rCM_color,'MarkerFaceAlpha',alphaValue,'LineWidth',0.25);
% scatter(r011X+7,r011Y,markerSize,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',rCM_color,'MarkerFaceAlpha',alphaValue,'LineWidth',0.25);
scatter(r111X+7,r111Y,markerSize,'MarkerEdgeColor',edgecolor,'MarkerFaceColor',rCM_color,'MarkerFaceAlpha',alphaValue,'LineWidth',0.25);

set(gca,'FontSize',fontSize,'LineWidth',1,'Color','None','TickDir','out','TickLength',[0.0125 0.025]);
ax = gca;
ax.XAxis.Color = 'k';
ax.YAxis.Color = 'k';
xlim([0.25 7.75]);
ylim(ylimits);
ylabel('Cell area (µm²)');
set(gca,'XTick',1:8,'XTickLabel',xticklabels);
set(gca,'Position',gcaPosition);
set(gca,'YTick',0:2000:16000);
ax = gca;
ax.YAxis.Exponent = 3;

%% Cell Elongation
ylimits = [1 6];
[e110X,e110Y] = beeswarmPoints([eCM.aacthcn4.Elongation],xspacing,(ylimits(2)-ylimits(1))/nBins,ylimits);
[e101X,e101Y] = beeswarmPoints([eCM.aactnppa.Elongation],xspacing,(ylimits(2)-ylimits(1))/nBins,ylimits);
[r100X,r100Y] = beeswarmPoints([rCM.aact.Elongation],xspacing,(ylimits(2)-ylimits(1))/nBins,ylimits);
[r110X,r110Y] = beeswarmPoints([rCM.aacthcn4.Elongation],xspacing,(ylimits(2)-ylimits(1))/nBins,ylimits);
[r101X,r101Y] = beeswarmPoints([rCM.aactnppa.Elongation],xspacing,(ylimits(2)-ylimits(1))/nBins,ylimits);
[r010X,r010Y] = beeswarmPoints([rCM.hcn4.Elongation],xspacing,(ylimits(2)-ylimits(1))/nBins,ylimits);
r011X = []; r011Y = [];
[r111X,r111Y] = beeswarmPoints([rCM.aacthcn4nppa.Elongation],xspacing,(ylimits(2)-ylimits(1))/nBins,ylimits);
figure('units','pixels','position',figureSize); hold on;
scatter(e110X+1,e110Y,markerSize,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',eCM_color,'MarkerFaceAlpha',alphaValue,'LineWidth',0.25);
scatter(e101X+2,e101Y,markerSize,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',eCM_color,'MarkerFaceAlpha',alphaValue,'LineWidth',0.25);
scatter(r100X+3,r100Y,markerSize,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',rCM_color,'MarkerFaceAlpha',alphaValue,'LineWidth',0.25);
scatter(r110X+4,r110Y,markerSize,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',rCM_color,'MarkerFaceAlpha',alphaValue,'LineWidth',0.25);
scatter(r101X+5,r101Y,markerSize,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',rCM_color,'MarkerFaceAlpha',alphaValue,'LineWidth',0.25);
scatter(r010X+6,r010Y,markerSize,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',rCM_color,'MarkerFaceAlpha',alphaValue,'LineWidth',0.25);
scatter(r011X+7,r011Y,markerSize,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',rCM_color,'MarkerFaceAlpha',alphaValue,'LineWidth',0.25);
scatter(r111X+8,r111Y,markerSize,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',rCM_color,'MarkerFaceAlpha',alphaValue,'LineWidth',0.25);

set(gca,'FontSize',fontSize,'LineWidth',1,'Color','None','TickDir','out','TickLength',[0.0125 0.025]);
ax = gca;
ax.XAxis.Color = 'k';
ax.YAxis.Color = 'k';
xlim([0.25 8.75]);
ylim(ylimits);
ylabel('Cell elongation');
set(gca,'XTick',1:8,'XTickLabel',xticklabels);
set(gca,'Position',gcaPosition);

%% Cell Eccentricity
ylimits = [0 1];
[e110X,e110Y] = beeswarmPoints([eCM.aacthcn4.Eccentricity],xspacing,(ylimits(2)-ylimits(1))/nBins,ylimits);
[e101X,e101Y] = beeswarmPoints([eCM.aactnppa.Eccentricity],xspacing,(ylimits(2)-ylimits(1))/nBins,ylimits);
[r100X,r100Y] = beeswarmPoints([rCM.aact.Eccentricity],xspacing,(ylimits(2)-ylimits(1))/nBins,ylimits);
[r110X,r110Y] = beeswarmPoints([rCM.aacthcn4.Eccentricity],xspacing,(ylimits(2)-ylimits(1))/nBins,ylimits);
[r101X,r101Y] = beeswarmPoints([rCM.aactnppa.Eccentricity],xspacing,(ylimits(2)-ylimits(1))/nBins,ylimits);
[r010X,r010Y] = beeswarmPoints([rCM.hcn4.Eccentricity],xspacing,(ylimits(2)-ylimits(1))/nBins,ylimits);
r011X = []; r011Y = [];
[r111X,r111Y] = beeswarmPoints([rCM.aacthcn4nppa.Eccentricity],xspacing,(ylimits(2)-ylimits(1))/nBins,ylimits);
figure('units','pixels','position',figureSize); hold on;
scatter(e110X+1,e110Y,markerSize,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',eCM_color,'MarkerFaceAlpha',alphaValue,'LineWidth',0.25);
scatter(e101X+2,e101Y,markerSize,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',eCM_color,'MarkerFaceAlpha',alphaValue,'LineWidth',0.25);
scatter(r100X+3,r100Y,markerSize,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',rCM_color,'MarkerFaceAlpha',alphaValue,'LineWidth',0.25);
scatter(r110X+4,r110Y,markerSize,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',rCM_color,'MarkerFaceAlpha',alphaValue,'LineWidth',0.25);
scatter(r101X+5,r101Y,markerSize,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',rCM_color,'MarkerFaceAlpha',alphaValue,'LineWidth',0.25);
scatter(r010X+6,r010Y,markerSize,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',rCM_color,'MarkerFaceAlpha',alphaValue,'LineWidth',0.25);
scatter(r011X+7,r011Y,markerSize,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',rCM_color,'MarkerFaceAlpha',alphaValue,'LineWidth',0.25);
scatter(r111X+8,r111Y,markerSize,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',rCM_color,'MarkerFaceAlpha',alphaValue,'LineWidth',0.25);

set(gca,'FontSize',fontSize,'LineWidth',1,'Color','None','TickDir','out','TickLength',[0.0125 0.025]);
ax = gca;
ax.XAxis.Color = 'k';
ax.YAxis.Color = 'k';
xlim([0.25 8.75]);
ylim(ylimits);
ylabel('Cell eccentricity');
set(gca,'XTick',1:8,'XTickLabel',xticklabels);
set(gca,'Position',gcaPosition);

%% Cell Circularity
ylimits = [0 1];
[e110X,e110Y] = beeswarmPoints([eCM.aacthcn4.Circularity],xspacing,(ylimits(2)-ylimits(1))/nBins,ylimits);
[e101X,e101Y] = beeswarmPoints([eCM.aactnppa.Circularity],xspacing,(ylimits(2)-ylimits(1))/nBins,ylimits);
[r100X,r100Y] = beeswarmPoints([rCM.aact.Circularity],xspacing,(ylimits(2)-ylimits(1))/nBins,ylimits);
[r110X,r110Y] = beeswarmPoints([rCM.aacthcn4.Circularity],xspacing,(ylimits(2)-ylimits(1))/nBins,ylimits);
[r101X,r101Y] = beeswarmPoints([rCM.aactnppa.Circularity],xspacing,(ylimits(2)-ylimits(1))/nBins,ylimits);
[r010X,r010Y] = beeswarmPoints([rCM.hcn4.Circularity],xspacing,(ylimits(2)-ylimits(1))/nBins,ylimits);
r011X = []; r011Y = [];
[r111X,r111Y] = beeswarmPoints([rCM.aacthcn4nppa.Circularity],xspacing,(ylimits(2)-ylimits(1))/nBins,ylimits);
figure('units','pixels','position',figureSize); hold on;
scatter(e110X+1,e110Y,markerSize,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',eCM_color,'MarkerFaceAlpha',alphaValue,'LineWidth',0.25);
scatter(e101X+2,e101Y,markerSize,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',eCM_color,'MarkerFaceAlpha',alphaValue,'LineWidth',0.25);
scatter(r100X+3,r100Y,markerSize,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',rCM_color,'MarkerFaceAlpha',alphaValue,'LineWidth',0.25);
scatter(r110X+4,r110Y,markerSize,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',rCM_color,'MarkerFaceAlpha',alphaValue,'LineWidth',0.25);
scatter(r101X+5,r101Y,markerSize,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',rCM_color,'MarkerFaceAlpha',alphaValue,'LineWidth',0.25);
scatter(r010X+6,r010Y,markerSize,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',rCM_color,'MarkerFaceAlpha',alphaValue,'LineWidth',0.25);
scatter(r011X+7,r011Y,markerSize,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',rCM_color,'MarkerFaceAlpha',alphaValue,'LineWidth',0.25);
scatter(r111X+8,r111Y,markerSize,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',rCM_color,'MarkerFaceAlpha',alphaValue,'LineWidth',0.25);

set(gca,'FontSize',fontSize,'LineWidth',1,'Color','None','TickDir','out','TickLength',[0.0125 0.025]);
ax = gca;
ax.XAxis.Color = 'k';
ax.YAxis.Color = 'k';
xlim([0.25 8.75]);
ylim(ylimits);
ylabel('Cell circularity');
set(gca,'XTick',1:8,'XTickLabel',xticklabels);
set(gca,'Position',gcaPosition);




%% Sarcomere Organization
ylimits = [0 0.6];
xlimits = [0.25 6.75];
[e110X,e110Y] = beeswarmPoints([eCM.aacthcn4.SarcomereOrganizationScore],xspacing,(ylimits(2)-ylimits(1))/nBins,ylimits);
[e101X,e101Y] = beeswarmPoints([eCM.aactnppa.SarcomereOrganizationScore],xspacing,(ylimits(2)-ylimits(1))/nBins,ylimits);
[r100X,r100Y] = beeswarmPoints([rCM.aact.SarcomereOrganizationScore],xspacing,(ylimits(2)-ylimits(1))/nBins,ylimits);
[r110X,r110Y] = beeswarmPoints([rCM.aacthcn4.SarcomereOrganizationScore],xspacing,(ylimits(2)-ylimits(1))/nBins,ylimits);
[r101X,r101Y] = beeswarmPoints([rCM.aactnppa.SarcomereOrganizationScore],xspacing,(ylimits(2)-ylimits(1))/nBins,ylimits);
[r111X,r111Y] = beeswarmPoints([rCM.aacthcn4nppa.SarcomereOrganizationScore],xspacing,(ylimits(2)-ylimits(1))/nBins,ylimits);
figure('units','pixels','position',figureSize2); hold on;
scatter(e110X+1,e110Y,markerSize,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',eCM_color,'MarkerFaceAlpha',alphaValue,'LineWidth',0.25);
scatter(e101X+2,e101Y,markerSize,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',eCM_color,'MarkerFaceAlpha',alphaValue,'LineWidth',0.25);
scatter(r100X+3,r100Y,markerSize,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',rCM_color,'MarkerFaceAlpha',alphaValue,'LineWidth',0.25);
scatter(r110X+4,r110Y,markerSize,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',rCM_color,'MarkerFaceAlpha',alphaValue,'LineWidth',0.25);
scatter(r101X+5,r101Y,markerSize,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',rCM_color,'MarkerFaceAlpha',alphaValue,'LineWidth',0.25);
scatter(r111X+6,r111Y,markerSize,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',rCM_color,'MarkerFaceAlpha',alphaValue,'LineWidth',0.25);
plot(xlimits,sarcomereOrganizationThreshold*[1 1],'--k');
set(gca,'FontSize',fontSize,'LineWidth',1,'Color','None','TickDir','out','TickLength',[0.0125 0.025]);
xlim(xlimits);
ylim(ylimits);
set(gca,'YTick',0:0.1:0.6);
ylabel('Sarcomere organization');
set(gca,'XTick',1:6,'XTickLabel',xticklabels);
set(gca,'Position',gcaPosition);
ax = gca;
ax.XAxis.Color = 'k';
ax.YAxis.Color = 'k';
%% Sarcomere Length
eCM_sarcOrg = [eCM.all.SarcomereOrganizationScore];
rCM_sarcOrg = [rCM.all_aact.SarcomereOrganizationScore];

e110_sarcOrg = [eCM.aacthcn4.SarcomereOrganizationScore];
e101_sarcOrg = [eCM.aactnppa.SarcomereOrganizationScore];
r100_sarcOrg = [rCM.aact.SarcomereOrganizationScore];
r110_sarcOrg = [rCM.aacthcn4.SarcomereOrganizationScore];
r101_sarcOrg = [rCM.aactnppa.SarcomereOrganizationScore];
r111_sarcOrg = [rCM.aacthcn4nppa.SarcomereOrganizationScore];

eCM_sarcLength = [eCM.all.SarcomereLength_um];
rCM_sarcLength = [rCM.all_aact.SarcomereLength_um];

e110_sarcLength = [eCM.aacthcn4.SarcomereLength_um];
e101_sarcLength = [eCM.aactnppa.SarcomereLength_um];
r100_sarcLength = [rCM.aact.SarcomereLength_um];
r110_sarcLength = [rCM.aacthcn4.SarcomereLength_um];
r101_sarcLength = [rCM.aactnppa.SarcomereLength_um];
r111_sarcLength = [rCM.aacthcn4nppa.SarcomereLength_um];

eCM_sarcLength = eCM_sarcLength(eCM_sarcOrg >= sarcomereOrganizationThreshold);
rCM_sarcLength = rCM_sarcLength(rCM_sarcOrg >= sarcomereOrganizationThreshold);

e110_sarcLength = e110_sarcLength(e110_sarcOrg >= sarcomereOrganizationThreshold);
e101_sarcLength = e101_sarcLength(e101_sarcOrg >= sarcomereOrganizationThreshold);
r100_sarcLength = r100_sarcLength(r100_sarcOrg >= sarcomereOrganizationThreshold);
r110_sarcLength = r110_sarcLength(r110_sarcOrg >= sarcomereOrganizationThreshold);
r101_sarcLength = r101_sarcLength(r101_sarcOrg >= sarcomereOrganizationThreshold);
r111_sarcLength = r111_sarcLength(r111_sarcOrg >= sarcomereOrganizationThreshold);

ylimits = [0 3];
a = 8/6.5;
[e110X,e110Y] = beeswarmPoints(e110_sarcLength,xspacing*a,(ylimits(2)-ylimits(1))/nBins,ylimits);
[e101X,e101Y] = beeswarmPoints(e101_sarcLength,xspacing*a,(ylimits(2)-ylimits(1))/nBins,ylimits);
[r100X,r100Y] = beeswarmPoints(r100_sarcLength,xspacing*a,(ylimits(2)-ylimits(1))/nBins,ylimits);
[r110X,r110Y] = beeswarmPoints(r110_sarcLength,xspacing*a,(ylimits(2)-ylimits(1))/nBins,ylimits);
[r101X,r101Y] = beeswarmPoints(r101_sarcLength,xspacing*a,(ylimits(2)-ylimits(1))/nBins,ylimits);
[r111X,r111Y] = beeswarmPoints(r111_sarcLength,xspacing*a,(ylimits(2)-ylimits(1))/nBins,ylimits);
figure('units','pixels','position',figureSize2); hold on;
scatter(e110X+1,e110Y,markerSize,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',eCM_color,'MarkerFaceAlpha',alphaValue,'LineWidth',0.25);
scatter(e101X+2.75,e101Y,markerSize,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',eCM_color,'MarkerFaceAlpha',alphaValue,'LineWidth',0.25);
scatter(r100X+4.5,r100Y,markerSize,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',rCM_color,'MarkerFaceAlpha',alphaValue,'LineWidth',0.25);
scatter(r110X+5.5,r110Y,markerSize,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',rCM_color,'MarkerFaceAlpha',alphaValue,'LineWidth',0.25);
scatter(r101X+6.5,r101Y,markerSize,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',rCM_color,'MarkerFaceAlpha',alphaValue,'LineWidth',0.25);
scatter(r111X+7.5,r111Y,markerSize,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',rCM_color,'MarkerFaceAlpha',alphaValue,'LineWidth',0.25);
set(gca,'FontSize',fontSize,'LineWidth',1,'Color','None','TickDir','out','TickLength',[0.0125 0.025]);
xlim([0.25 8.25]);
ylim(ylimits);
ylabel('Sarcomere length (µm)');
set(gca,'YTick',0:0.5:3);
set(gca,'XTick',[1 2.75 4.5 5.5 6.5 7.5],'XTickLabel',xticklabels);
set(gca,'Position',gcaPosition);
ax = gca;
ax.XAxis.Color = 'k';
ax.YAxis.Color = 'k';
%% Cell-Sarcomere Misalignment
eCM_sarcOrg = [eCM.all.SarcomereOrganizationScore];
rCM_sarcOrg = [rCM.all_aact.SarcomereOrganizationScore];

e110_sarcOrg = [eCM.aacthcn4.SarcomereOrganizationScore];
e101_sarcOrg = [eCM.aactnppa.SarcomereOrganizationScore];
r100_sarcOrg = [rCM.aact.SarcomereOrganizationScore];
r110_sarcOrg = [rCM.aacthcn4.SarcomereOrganizationScore];
r101_sarcOrg = [rCM.aactnppa.SarcomereOrganizationScore];
r111_sarcOrg = [rCM.aacthcn4nppa.SarcomereOrganizationScore];

eCM_elong = [eCM.all.Elongation];
rCM_elong = [rCM.all_aact.Elongation];

e110_elong = [eCM.aacthcn4.Elongation];
e101_elong = [eCM.aactnppa.Elongation];
r100_elong = [rCM.aact.Elongation];
r110_elong = [rCM.aacthcn4.Elongation];
r101_elong = [rCM.aactnppa.Elongation];
r111_elong = [rCM.aacthcn4nppa.Elongation];

eCM_sarcMisalignment = [eCM.all.CellSarcomereAlignment];
rCM_sarcMisalignment = [rCM.all_aact.CellSarcomereAlignment];

e110_sarcMisalignment = [eCM.aacthcn4.CellSarcomereAlignment];
e101_sarcMisalignment = [eCM.aactnppa.CellSarcomereAlignment];
r100_sarcMisalignment = [rCM.aact.CellSarcomereAlignment];
r110_sarcMisalignment = [rCM.aacthcn4.CellSarcomereAlignment];
r101_sarcMisalignment = [rCM.aactnppa.CellSarcomereAlignment];
r111_sarcMisalignment = [rCM.aacthcn4nppa.CellSarcomereAlignment];

eCM_sarcMisalignment = eCM_sarcMisalignment(eCM_sarcOrg > sarcomereOrganizationThreshold & eCM_elong > elongationThreshold);
rCM_sarcMisalignment = rCM_sarcMisalignment(rCM_sarcOrg > sarcomereOrganizationThreshold & rCM_elong > elongationThreshold);

e110_sarcMisalignment = e110_sarcMisalignment(e110_sarcOrg > sarcomereOrganizationThreshold & e110_elong > elongationThreshold);
e101_sarcMisalignment = e101_sarcMisalignment(e101_sarcOrg > sarcomereOrganizationThreshold & e101_elong > elongationThreshold);
r100_sarcMisalignment = r100_sarcMisalignment(r100_sarcOrg > sarcomereOrganizationThreshold & r100_elong > elongationThreshold);
r110_sarcMisalignment = r110_sarcMisalignment(r110_sarcOrg > sarcomereOrganizationThreshold & r110_elong > elongationThreshold);
r101_sarcMisalignment = r101_sarcMisalignment(r101_sarcOrg > sarcomereOrganizationThreshold & r101_elong > elongationThreshold);
r111_sarcMisalignment = r111_sarcMisalignment(r111_sarcOrg > sarcomereOrganizationThreshold & r111_elong > elongationThreshold);

ylimits = [0 90];
[eX,eY] = beeswarmPoints(eCM_sarcMisalignment,xspacing,(ylimits(2)-ylimits(1))/nBins,ylimits);
[rX,rY] = beeswarmPoints(rCM_sarcMisalignment,xspacing,(ylimits(2)-ylimits(1))/nBins,ylimits);

a = 1;
[e110X,e110Y] = beeswarmPoints(e110_sarcMisalignment,xspacing*a,(ylimits(2)-ylimits(1))/nBins,ylimits);
[e101X,e101Y] = beeswarmPoints(e101_sarcMisalignment,xspacing*a,(ylimits(2)-ylimits(1))/nBins,ylimits);
[r100X,r100Y] = beeswarmPoints(r100_sarcMisalignment,xspacing*a,(ylimits(2)-ylimits(1))/nBins,ylimits);
[r110X,r110Y] = beeswarmPoints(r110_sarcMisalignment,xspacing*a,(ylimits(2)-ylimits(1))/nBins,ylimits);
[r101X,r101Y] = beeswarmPoints(r101_sarcMisalignment,xspacing*a,(ylimits(2)-ylimits(1))/nBins,ylimits);
[r111X,r111Y] = beeswarmPoints(r111_sarcMisalignment,xspacing*a,(ylimits(2)-ylimits(1))/nBins,ylimits);

figure('units','pixels','position',figureSize2); hold on;
scatter(e110X+1,e110Y,markerSize,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',eCM_color,'MarkerFaceAlpha',alphaValue,'LineWidth',0.25);
scatter(e101X+2,e101Y,markerSize,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',eCM_color,'MarkerFaceAlpha',alphaValue,'LineWidth',0.25);
scatter(r100X+3,r100Y,markerSize,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',rCM_color,'MarkerFaceAlpha',alphaValue,'LineWidth',0.25);
scatter(r110X+4,r110Y,markerSize,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',rCM_color,'MarkerFaceAlpha',alphaValue,'LineWidth',0.25);
scatter(r101X+5,r101Y,markerSize,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',rCM_color,'MarkerFaceAlpha',alphaValue,'LineWidth',0.25);
scatter(r111X+6,r111Y,markerSize,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',rCM_color,'MarkerFaceAlpha',alphaValue,'LineWidth',0.25);
set(gca,'FontSize',fontSize,'LineWidth',1,'Color','None','TickDir','out','TickLength',[0.0125 0.025]);
xlim([0.25 6.75]);
ylim(ylimits);
ylabel('Cell-sarcomere misalignment');
set(gca,'XTick',1:6,'XTickLabel',xticklabels);
set(gca,'YTick',0:10:90,'YTickLabel',{'0°','10°','20°','30°','40°','50°','60°','70°','80°','90°'});
set(gca,'Position',gcaPosition);
ax = gca;
ax.XAxis.Color = 'k';
ax.YAxis.Color = 'k';
%% Cell Elongation vs Cell Area
xlimits = [0 16000];
ylimits = [1 6];

% Endogenous CM
figure('units','pixels','position',figureSize); hold on;
scatter([eCM.all.Area_um2],[eCM.all.Elongation],markerSizeScatter,'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',eCM_color,'MarkerFaceAlpha',alphaValue);
set(gca,'FontSize',scatterFontSize,'LineWidth',4,'Color','None','TickDir','out','TickLength',[0.02 0.025]);
xlim(xlimits);
ylim(ylimits);
xlabel('Cell Area (µm²)');
ylabel('Cell Elongation');
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
plot(xfit,yfit,'LineWidth',4,'Color',fitLineColor);

P = fitlm([eCM.all.Area_um2],[eCM.all.Elongation],'linear');

text(xlimits(2),ylimits(2),['R²=' sprintf('%.2f',P.Rsquared.Ordinary)],'FontSize',32,'FontName','Arial','HorizontalAlignment','Right','VerticalAlignment','Top');

% Reprogrammed CM
figure('units','pixels','position',figureSize); hold on;
scatter([rCM.all.Area_um2],[rCM.all.Elongation],markerSizeScatter,'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',rCM_color,'MarkerFaceAlpha',alphaValue);
set(gca,'FontSize',scatterFontSize,'LineWidth',4,'Color','None','TickDir','out','TickLength',[0.02 0.025]);
xlim(xlimits);
ylim(ylimits);
xlabel('Cell Area (µm²)');
ylabel('Cell Elongation');
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
plot(xfit,yfit,'LineWidth',4,'Color',fitLineColor);

P = fitlm([rCM.all.Area_um2],[rCM.all.Elongation],'linear');

text(xlimits(2),ylimits(2),['R²=' sprintf('%.2f',P.Rsquared.Ordinary)],'FontSize',32,'FontName','Arial','HorizontalAlignment','Right','VerticalAlignment','Top');

%% Cell Area vs Sarcomere Organization
xlimits = [0 0.6];
ylimits = [0 16000];

% Endogenous CM
figure('units','pixels','position',figureSize); hold on;
scatter([eCM.all.SarcomereOrganizationScore],[eCM.all.Area_um2],markerSizeScatter,'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',eCM_color,'MarkerFaceAlpha',alphaValue);
set(gca,'FontSize',scatterFontSize,'LineWidth',4,'Color','None','TickDir','out','TickLength',[0.02 0.025]);
xlim(xlimits);
ylim(ylimits);
xlabel('Sarcomere Organization');
ylabel('Cell Area (µm²)');
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
plot(xfit,yfit,'LineWidth',4,'Color',fitLineColor);

P = fitlm([eCM.all.SarcomereOrganizationScore],[eCM.all.Area_um2],'linear');

text(xlimits(2),ylimits(2),['R²=' sprintf('%.2f',P.Rsquared.Ordinary)],'FontSize',32,'FontName','Arial','HorizontalAlignment','Right','VerticalAlignment','Top');

% Reprogrammed CM
figure('units','pixels','position',figureSize); hold on;
scatter([rCM.all.SarcomereOrganizationScore],[rCM.all.Area_um2],markerSizeScatter,'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',rCM_color,'MarkerFaceAlpha',alphaValue);
set(gca,'FontSize',scatterFontSize,'LineWidth',4,'Color','None','TickDir','out','TickLength',[0.02 0.025]);
xlim(xlimits);
ylim(ylimits);
xlabel('Sarcomere Organization');
ylabel('Cell Area (µm²)');
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
plot(xfit,yfit,'LineWidth',4,'Color',fitLineColor);

P = fitlm([rCM.all.SarcomereOrganizationScore],[rCM.all.Area_um2],'linear');

text(xlimits(2),ylimits(2),['R²=' sprintf('%.2f',P.Rsquared.Ordinary)],'FontSize',32,'FontName','Arial','HorizontalAlignment','Right','VerticalAlignment','Top');

%% Cell Elongation vs Sarcomere Organization
xlimits = [0 0.6];
ylimits = [1 6];

% Endogenous CM
figure('units','pixels','position',figureSize); hold on;
scatter([eCM.all.SarcomereOrganizationScore],[eCM.all.Elongation],markerSizeScatter,'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',eCM_color,'MarkerFaceAlpha',alphaValue);
set(gca,'FontSize',scatterFontSize,'LineWidth',4,'Color','None','TickDir','out','TickLength',[0.02 0.025]);
xlim(xlimits);
ylim(ylimits);
ax = gca;
ax.XAxis.Color = 'k';
ax.YAxis.Color = 'k';
xlabel('Sarcomere Organization');
ylabel('Cell Elongation');
set(gca,'Position',scattergcaPosition);

X = [ones(length([eCM.all.SarcomereOrganizationScore]),1),[eCM.all.SarcomereOrganizationScore]'];
k = X\[eCM.all.Elongation]';
xfit = [min([eCM.all.SarcomereOrganizationScore]) max([eCM.all.SarcomereOrganizationScore])];
yfit = k(1)+k(2).*xfit;
plot(xfit,yfit,'LineWidth',4,'Color',fitLineColor);

P = fitlm([eCM.all.SarcomereOrganizationScore],[eCM.all.Elongation],'linear');

text(xlimits(2),ylimits(2),['R²=' sprintf('%.2f',P.Rsquared.Ordinary)],'FontSize',32,'FontName','Arial','HorizontalAlignment','Right','VerticalAlignment','Top');

% Reprogrammed CM
figure('units','pixels','position',figureSize); hold on;
scatter([rCM.all.SarcomereOrganizationScore],[rCM.all.Elongation],markerSizeScatter,'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',rCM_color,'MarkerFaceAlpha',alphaValue);
set(gca,'FontSize',scatterFontSize,'LineWidth',4,'Color','None','TickDir','out','TickLength',[0.02 0.025]);
xlim(xlimits);
ylim(ylimits);
xlabel('Sarcomere Organization');
ylabel('Cell Elongation');
set(gca,'Position',scattergcaPosition);
ax = gca;
ax.XAxis.Color = 'k';
ax.YAxis.Color = 'k';
X = [ones(length([rCM.all.SarcomereOrganizationScore]),1),[rCM.all.SarcomereOrganizationScore]'];
k = X\[rCM.all.Elongation]';
xfit = [min([rCM.all.SarcomereOrganizationScore]) max([rCM.all.SarcomereOrganizationScore])];
yfit = k(1)+k(2).*xfit;
plot(xfit,yfit,'LineWidth',4,'Color',fitLineColor);

P = fitlm([rCM.all.SarcomereOrganizationScore],[rCM.all.Elongation],'linear');

text(xlimits(2),ylimits(2),['R²=' sprintf('%.2f',P.Rsquared.Ordinary)],'FontSize',32,'FontName','Arial','HorizontalAlignment','Right','VerticalAlignment','Top');

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
figure('units','pixels','position',figureSize); hold on;
scatter(eCM_sarcOrg_disorg,eCM_sarcLength_disorg,markerSizeScatter*2/3,'MarkerEdgeColor',eCM_color,'MarkerFaceColor','none','MarkerFaceAlpha',alphaValue,'LineWidth',1);
scatter(eCM_sarcOrg_org,eCM_sarcLength_org,markerSizeScatter,'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',eCM_color,'MarkerFaceAlpha',alphaValue);
area([-10 sarcomereOrganizationThreshold],[10 10],-10,'FaceColor',[0.9 0.9 0.9],'FaceAlpha',1,'LineStyle','--','LineWidth',3);
set(gca,'FontSize',scatterFontSize,'LineWidth',4,'Color','None','TickDir','out','TickLength',[0.02 0.025]);
xlim(xlimits);
ylim(ylimits);
xlabel('Sarcomere Organization');
ylabel('Sarcomere Length (µm)');
set(gca,'Position',scattergcaPosition);
ax = gca;
ax.XAxis.Color = 'k';
ax.YAxis.Color = 'k';
X = [ones(length(eCM_sarcOrg_org),1),eCM_sarcOrg_org'];
k = X\eCM_sarcLength_org';
xfit = [min(eCM_sarcOrg_org) max(eCM_sarcOrg_org)];
yfit = k(1)+k(2).*xfit;
plot(xfit,yfit,'LineWidth',4,'Color',fitLineColor);

P = fitlm(eCM_sarcOrg_org,eCM_sarcLength_org,'linear');

text(xlimits(2),ylimits(2),['R²=' sprintf('%.2f',P.Rsquared.Ordinary)],'FontSize',32,'FontName','Arial','HorizontalAlignment','Right','VerticalAlignment','Top');

% Reprogrammed CM
figure('units','pixels','position',figureSize); hold on;
scatter(rCM_sarcOrg_disorg,rCM_sarcLength_disorg,markerSizeScatter*2/3,'MarkerEdgeColor',rCM_color,'MarkerFaceColor','none','MarkerFaceAlpha',alphaValue);
scatter(rCM_sarcOrg_org,rCM_sarcLength_org,markerSizeScatter,'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',rCM_color,'MarkerFaceAlpha',alphaValue);
area([-10 sarcomereOrganizationThreshold],[10 10],-10,'FaceColor',[0.9 0.9 0.9],'FaceAlpha',1,'LineStyle','--','LineWidth',3);
set(gca,'FontSize',scatterFontSize,'LineWidth',4,'Color','None','TickDir','out','TickLength',[0.02 0.025]);
xlim(xlimits);
ylim(ylimits);
xlabel('Sarcomere Organization');
ylabel('Sarcomere Length (µm)');
set(gca,'Position',scattergcaPosition);
ax = gca;
ax.XAxis.Color = 'k';
ax.YAxis.Color = 'k';
X = [ones(length(rCM_sarcOrg_org),1),rCM_sarcOrg_org'];
k = X\rCM_sarcLength_org';
xfit = [min(rCM_sarcOrg_org) max(rCM_sarcOrg_org)];
yfit = k(1)+k(2).*xfit;
plot(xfit,yfit,'LineWidth',4,'Color',fitLineColor);

P = fitlm(rCM_sarcOrg_org,rCM_sarcLength_org,'linear');

text(xlimits(2),ylimits(2),['R²=' sprintf('%.2f',P.Rsquared.Ordinary)],'FontSize',32,'FontName','Arial','HorizontalAlignment','Right','VerticalAlignment','Top');

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
figure('units','pixels','position',figureSize); hold on;
scatter(eCM_CellSarcomereAlignment_disorg,eCM_Elongation_disorg,markerSizeScatter,'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',belowThreshScatterColor,'MarkerFaceAlpha',alphaValue);
scatter(eCM_CellSarcomereAlignment_org,eCM_Elongation_org,markerSizeScatter,'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',eCM_color,'MarkerFaceAlpha',alphaValue);

set(gca,'FontSize',scatterFontSize,'LineWidth',4,'Color','None','TickDir','out','TickLength',[0.02 0.025]);
xlim(xlimits);
ylim(ylimits);
xlabel('Cell-Sarcomere Misalignment');
ylabel('Cell Elongation');
set(gca,'Position',scattergcaPosition);
set(gca,'XTick',0:30:90,'XTickLabel',{'0°','30°','60°','90°'});
ax = gca;
ax.XAxis.Color = 'k';
ax.YAxis.Color = 'k';
X = [ones(length(eCM_CellSarcomereAlignment_org),1),eCM_CellSarcomereAlignment_org'];
k = X\eCM_Elongation_org';
xfit = [min(eCM_CellSarcomereAlignment_org) max(eCM_CellSarcomereAlignment_org)];
yfit = k(1)+k(2).*xfit;
plot(xfit,yfit,'LineWidth',4,'Color',fitLineColor);

P = fitlm(eCM_CellSarcomereAlignment_org,eCM_Elongation_org,'linear');

text(xlimits(2),ylimits(2),['R²=' sprintf('%.2f',P.Rsquared.Ordinary)],'FontSize',32,'FontName','Arial','HorizontalAlignment','Right','VerticalAlignment','Top');

% Reprogrammed CM
figure('units','pixels','position',figureSize); hold on;
scatter(rCM_CellSarcomereAlignment_disorg,rCM_Elongation_disorg,markerSizeScatter,'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',belowThreshScatterColor,'MarkerFaceAlpha',alphaValue);
scatter(rCM_CellSarcomereAlignment_org,rCM_Elongation_org,markerSizeScatter,'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',rCM_color,'MarkerFaceAlpha',alphaValue);

set(gca,'FontSize',scatterFontSize,'LineWidth',4,'Color','None','TickDir','out','TickLength',[0.02 0.025]);
xlim(xlimits);
ylim(ylimits);
xlabel('Cell-Sarcomere Misalignment');
ylabel('Cell Elongation');
set(gca,'Position',scattergcaPosition);
set(gca,'XTick',0:30:90,'XTickLabel',{'0°','30°','60°','90°'});
ax = gca;
ax.XAxis.Color = 'k';
ax.YAxis.Color = 'k';
X = [ones(length(rCM_CellSarcomereAlignment_org),1),rCM_CellSarcomereAlignment_org'];
k = X\rCM_Elongation_org';
xfit = [min(rCM_CellSarcomereAlignment_org) max(rCM_CellSarcomereAlignment_org)];
yfit = k(1)+k(2).*xfit;
plot(xfit,yfit,'LineWidth',4,'Color',fitLineColor);

P = fitlm(rCM_CellSarcomereAlignment_org,rCM_Elongation_org,'linear');

text(xlimits(2),ylimits(2),['R²=' sprintf('%.2f',P.Rsquared.Ordinary)],'FontSize',32,'FontName','Arial','HorizontalAlignment','Right','VerticalAlignment','Top');