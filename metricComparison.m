%% metricComparison
% This script calculates all potential metrics (Haralick features, Gabor
% features, and Fourier transform features) for organized and disorganized
% cells, and plots them.
readFolder = '.\Images\Metric Comparison\';

nOrganized = 32;
nDisorganized = 26;

[organizedHaralickCorrelation,organizedHaralickContrast,...
    organizedHaralickEnergy,organizedHaralickHomogeneity,...
    organizedHaralickVariance,organizedGaborScore,organizedFourierScore] = deal(zeros(1,nOrganized));

[organizedRadialProfile,organizedFourierTransform,...
    organizedFourierFitExp1,organizedFourierFitExp2,...
    organizedFourierFitGauss,organizedFourierFitParams] = deal(cell(1,nOrganized));
warning('off','all');
for iImage = 1:nOrganized
    disp(['Organized Image #' num2str(iImage)]);
    iI = imread([readFolder 'Organized ' sprintf('%02d',iImage) '.tif']);
    iM = imread([readFolder 'Organized ' sprintf('%02d',iImage) ' Mask.tif']);
    
    iI = im2double(iI);
    iR = iI(:,:,1);
    iM = iM(:,:,1);
    
    iMetrics = morph_texture_function_allmetrics(iR,iM,0:5:175,0:1:40,1,1,readFolder);
    
    organizedHaralickCorrelation(iImage) = iMetrics.SarcomereOrganizationScore;
    organizedHaralickContrast(iImage) = iMetrics.SarcomereOrganizationScoreContrast;
    organizedHaralickEnergy(iImage) = iMetrics.SarcomereOrganizationScoreEnergy;
    organizedHaralickHomogeneity(iImage) = iMetrics.SarcomereOrganizationScoreHomogeneity;
    organizedHaralickVariance(iImage) = iMetrics.SarcomereOrganizationScoreVariance;
    organizedGaborScore(iImage) = iMetrics.GaborScore;
    organizedFourierScore(iImage) = iMetrics.FourierScore;
    organizedRadialProfile{iImage} = iMetrics.RadialProfile;
    organizedFourierTransform{iImage} = iMetrics.FourierTransform;
    organizedFourierFitExp1{iImage} = iMetrics.FitExp1;
    organizedFourierFitExp2{iImage} = iMetrics.FitExp2;
    organizedFourierFitGauss{iImage} = iMetrics.FitGauss;
    organizedFourierFitParams{iImage} = iMetrics.FitParameters;
end


[disorganizedHaralickCorrelation,disorganizedHaralickContrast,...
    disorganizedHaralickEnergy,disorganizedHaralickHomogeneity,...
    disorganizedHaralickVariance,disorganizedGaborScore,disorganizedFourierScore] = deal(zeros(1,nDisorganized));


[disorganizedRadialProfile,disorganizedFourierTransform,...
    disorganizedFourierFitExp1,disorganizedFourierFitExp2,...
    disorganizedFourierFitGauss,disorganizedFourierFitParams] = deal(cell(1,nDisorganized));
for iImage = 1:nDisorganized
    disp(['Disorganized Image #' num2str(iImage)]);
    iI = imread([readFolder 'Disorganized ' sprintf('%02d',iImage) '.tif']);
    iM = imread([readFolder 'Disorganized ' sprintf('%02d',iImage) ' Mask.tif']);
    
    iI = im2double(iI);
    iR = iI(:,:,1);
    iM = iM(:,:,1);
    
    iMetrics = morph_texture_function_allmetrics(iR,iM,0:5:175,0:1:40,1,1,readFolder);
    
    disorganizedHaralickCorrelation(iImage) = iMetrics.SarcomereOrganizationScore;
    disorganizedHaralickContrast(iImage) = iMetrics.SarcomereOrganizationScoreContrast;
    disorganizedHaralickEnergy(iImage) = iMetrics.SarcomereOrganizationScoreEnergy;
    disorganizedHaralickHomogeneity(iImage) = iMetrics.SarcomereOrganizationScoreHomogeneity;
    disorganizedHaralickVariance(iImage) = iMetrics.SarcomereOrganizationScoreVariance;
    disorganizedGaborScore(iImage) = iMetrics.GaborScore;
    disorganizedFourierScore(iImage) = iMetrics.FourierScore;
    disorganizedRadialProfile{iImage} = iMetrics.RadialProfile;
    disorganizedFourierTransform{iImage} = iMetrics.FourierTransform;
    disorganizedFourierFitExp1{iImage} = iMetrics.FitExp1;
    disorganizedFourierFitExp2{iImage} = iMetrics.FitExp2;
    disorganizedFourierFitGauss{iImage} = iMetrics.FitGauss;
    disorganizedFourierFitParams{iImage} = iMetrics.FitParameters;
end

clear i*
%%
save('data_metriccomparison.mat');
%%
[~,pCorrelation] = ttest2(disorganizedHaralickCorrelation,organizedHaralickCorrelation);
[~,pContrast] = ttest2(disorganizedHaralickContrast,organizedHaralickContrast);
[~,pEnergy] = ttest2(disorganizedHaralickEnergy,organizedHaralickEnergy);
[~,pHomogeneity] = ttest2(disorganizedHaralickHomogeneity,organizedHaralickHomogeneity);
[~,pVariance] = ttest2(disorganizedHaralickVariance,organizedHaralickVariance);
[~,pGabor] = ttest2(disorganizedGaborScore,organizedGaborScore);
[~,pFourier] = ttest2(disorganizedFourierScore,organizedFourierScore);

[xCorrelation,yCorrelation,~,AUC_Correlation] = perfcurve([zeros(1,nDisorganized) ones(1,nOrganized)],[disorganizedHaralickCorrelation organizedHaralickCorrelation],1);
[xContrast,yContrast,~,AUC_Contrast] = perfcurve([zeros(1,nDisorganized) ones(1,nOrganized)],[disorganizedHaralickContrast organizedHaralickContrast],1);
[xEnergy,yEnergy,~,AUC_Energy] = perfcurve([zeros(1,nDisorganized) ones(1,nOrganized)],[disorganizedHaralickEnergy organizedHaralickEnergy],1);
[xHomogeneity,yHomogeneity,~,AUC_Homogeneity] = perfcurve([zeros(1,nDisorganized) ones(1,nOrganized)],[disorganizedHaralickHomogeneity organizedHaralickHomogeneity],1);
[xVariance,yVariance,~,AUC_Variance] = perfcurve([zeros(1,nDisorganized) ones(1,nOrganized)],[disorganizedHaralickVariance organizedHaralickVariance],1);
[xGabor,yGabor,~,AUC_Gabor] = perfcurve([zeros(1,nDisorganized) ones(1,nOrganized)],[disorganizedGaborScore organizedGaborScore],1);
[xFourier,yFourier,~,AUC_Fourier] = perfcurve([zeros(1,nDisorganized) ones(1,nOrganized)],[disorganizedFourierScore organizedFourierScore],1);


%%
close all;
plotBarGraph(disorganizedHaralickCorrelation,organizedHaralickCorrelation)
ylabel('Haralick Correlation Score');
%%
close all;
plotBarGraph(disorganizedHaralickContrast,organizedHaralickContrast)
ylabel('Haralick Contrast Score');
%%
close all;
plotBarGraph(disorganizedHaralickEnergy,organizedHaralickEnergy)
ylabel('Haralick Uniformity Score');
ax = gca;
ax.YAxis.Exponent = -3;
%%
close all;
plotBarGraph(disorganizedHaralickHomogeneity,organizedHaralickHomogeneity)
ylabel('Haralick Homogeneity Score');
ax = gca;
ax.YAxis.Exponent = -2;
%%
close all;
plotBarGraph(disorganizedHaralickVariance,organizedHaralickVariance)
ylabel('Variance Score');
%%
close all;
plotBarGraph(disorganizedGaborScore,organizedGaborScore)
ylabel('Gabor Filter Score');
%%
close all;
plotBarGraph(disorganizedFourierScore,organizedFourierScore)
ylabel('Fourier Transform Score');

%%
close all;
figure('units','pixels','position',[50 50 700 600]); hold on;
set(gca,'LineWidth',2,'Color','none','FontSize',16,'TickDir','out');
box off;
plot([0 1],[0 1],'Color',[0.5 0.5 0.5]);
p1 = plot(xCorrelation,yCorrelation,'LineWidth',2);
p2 = plot(xContrast,yContrast,'LineWidth',2);
p3 = plot(xEnergy,yEnergy,'LineWidth',2);
p4 = plot(xHomogeneity,yHomogeneity,'LineWidth',2);
p5 = plot(xVariance,yVariance,'LineWidth',2);
p6 = plot(xGabor,yGabor,'LineWidth',2);
p7 = plot(xFourier,yFourier,'LineWidth',2);
xlim([0 1]);
ylim([0 1]);
axis square;
xlabel('False Positive Rate');
ylabel('True Positive Rate');
legend([p1 p2 p3 p4 p5 p6 p7],{'Haralick Correlation, AUC = 0.97','Haralick Contrast, AUC = 0.97','Haralick Uniformity, AUC = 0.86','Haralick Homogeneity, AUC = 0.97','Haralick Variance, AUC = 0.72','Gabor Wavelet, AUC = 0.78','Fourier Transform, AUC = 0.90'},'Location','SouthEast','FontSize',12);
legend('boxoff');