function metrics = morph_texture_function_allmetrics(aact,mask,angleSet,spatialSet,d1_conversion,d2_conversion,rootname)
%% FUNCTION DESCRIPTION
% This function computes cell morphology features (cell area, elongation,
% circularity, eccentricity) and all potential sarcomere organization
% features (sarcomere organization score, sarcomere length, cell-sarcomere
% misalignment) for every cell in the image.
%
% INPUTS
% aact : alpha-actinin channel image
% mask : segmentation image with same dimensions as aact. Each pixel is labeled as belonging to a cell
%   (1 = cell #1, 2 = cell #2, etc.) or background = 0.
% angleSet : array of orientations at which to measure Haralick
%   correlation, typically 0:1:179
% spatialSet : array of pixel offsets at which to measure Haralick
%   correlation, typically 0:1:40
% d1_conversion : pixels to microns conversion for a single dimension
% d2_conversion : pixels to microns conversion for 2 dimensions (i.e. area)
% rootname : image file root name for aact and mask
%
% OUTPUTS
% metrics : structure with all measured features, each of which contains an
%   array, where each element is the measured feature for the indexed cell
%% Split image
metrics = regionprops(mask,aact,'all');

for iCell = 1:numel(metrics)
    metrics(iCell).RootName = rootname;
    
    metrics(iCell).Area_um2 = metrics(iCell).Area.*d2_conversion;
    metrics(iCell).MajorAxisLength_um = metrics(iCell).MajorAxisLength.*d1_conversion;
    metrics(iCell).MinorAxisLength_um = metrics(iCell).MinorAxisLength.*d1_conversion;
    
    iMask = (mask == iCell);
    
    % Crop image to only include entire cell
    iaact = aact(:,:,1);
    iaact(~any(iMask,2),:) = [];
    iaact(:,~any(iMask,1)) = [];
    
    iMask(~any(iMask,2),:) = [];
    iMask(:,~any(iMask,1)) = [];
    
    metrics(iCell).AACT_BoundingBox = iaact;
    
    % Calculate elongation and circularity
    metrics(iCell).Elongation = metrics(iCell).MajorAxisLength./metrics(iCell).MinorAxisLength;
    metrics(iCell).Circularity = (4.*pi.*metrics(iCell).Area)./(metrics(iCell).Perimeter)^2;
    
    % Set up for sub-pixel resolution interpolation
    spatialSetInterp = 0:0.1:max(spatialSet);
    
    % Create offset matrix for co-occurrence matrix calculation
    offsets = [zeros(numel(spatialSet),1) spatialSet'; -spatialSet' zeros(numel(spatialSet),1)];
    
    % Initialize m x n matrix, where m = number of angles and n =
    % number of interpolated spatial pairs
    correlationdata = zeros(numel(angleSet),numel(spatialSetInterp));
    contrastdata = zeros(numel(angleSet),numel(spatialSetInterp));
    energydata = zeros(numel(angleSet),numel(spatialSetInterp));
    homogeneitydata = zeros(numel(angleSet),numel(spatialSetInterp));
    % For efficiency, calculate together angles that are 90 degrees apart
    halfangleSet = angleSet(1:(numel(angleSet)/2));
    
    % Calculate correlation for each angle and spatial distance
    angleIndex = 1;
    for iAngle = halfangleSet
        rotateI = imrotate(iaact,-iAngle);
        rotateM = imrotate(iMask,-iAngle);
        
        rotateI(~rotateM) = NaN;
        
        glcm = graycomatrix(rotateI,'offset',offsets);
        stats = graycoprops(glcm);
        stats.Correlation(isnan(stats.Correlation)) = 0;
        stats.Contrast(isnan(stats.Contrast)) = 0;
        stats.Energy(isnan(stats.Energy)) = 0;
        stats.Homogeneity(isnan(stats.Homogeneity)) = 0;
        
        correlationdata(angleIndex,:) = interp1(spatialSet,stats.Correlation(1:numel(spatialSet)),spatialSetInterp,'spline');
        correlationdata(angleIndex+numel(angleSet)/2,:) = interp1(spatialSet,stats.Correlation(numel(spatialSet)+1:end),spatialSetInterp,'spline');
        
        contrastdata(angleIndex,:) = interp1(spatialSet,stats.Contrast(1:numel(spatialSet)),spatialSetInterp,'spline');
        contrastdata(angleIndex+numel(angleSet)/2,:) = interp1(spatialSet,stats.Contrast(numel(spatialSet)+1:end),spatialSetInterp,'spline');
        
        energydata(angleIndex,:) = interp1(spatialSet,stats.Energy(1:numel(spatialSet)),spatialSetInterp,'spline');
        energydata(angleIndex+numel(angleSet)/2,:) = interp1(spatialSet,stats.Energy(numel(spatialSet)+1:end),spatialSetInterp,'spline');
        
        homogeneitydata(angleIndex,:) = interp1(spatialSet,stats.Homogeneity(1:numel(spatialSet)),spatialSetInterp,'spline');
        homogeneitydata(angleIndex+numel(angleSet)/2,:) = interp1(spatialSet,stats.Homogeneity(numel(spatialSet)+1:end),spatialSetInterp,'spline');
        
        
        angleIndex = angleIndex + 1;
        
    end
    
    metrics(iCell).CorrelationData = correlationdata;
    metrics(iCell).ContrastData = contrastdata;
    metrics(iCell).EnergyData = energydata;
    metrics(iCell).HomogeneityData = homogeneitydata;
    
    amplitude = zeros(size(correlationdata,1),1);
    lengthIndex = zeros(size(correlationdata,1),1);
    
    for iAngle = 1:size(correlationdata,1)
        [pks,locs,~,p] = findpeaks(correlationdata(iAngle,:));
        
        if ~isempty(pks)
            [pMax,locMax] = max(p);
            amplitude(iAngle) = pMax;
            lengthIndex(iAngle) = spatialSetInterp(locs(locMax));
        end
    end
    
    [sarcOrgScore,sarcIndex] = max(amplitude);
    
    metrics(iCell).SarcomereOrganizationScore = sarcOrgScore;
    metrics(iCell).SarcomereOrientation = angleSet(sarcIndex);
    metrics(iCell).CellSarcomereAlignment = min(180-mod(metrics(iCell).Orientation-metrics(iCell).SarcomereOrientation,180),...
        mod(metrics(iCell).Orientation-metrics(iCell).SarcomereOrientation,180));
    
    if lengthIndex(sarcIndex) > 0
        metrics(iCell).SarcomereLength = lengthIndex(sarcIndex);
    else
        metrics(iCell).SarcomereLength = 0;
    end
    metrics(iCell).SarcomereLength_um = metrics(iCell).SarcomereLength.*d1_conversion;
    
    if sarcIndex+floor(numel(angleSet)./2) <= numel(angleSet)
        widthIndex = sarcIndex+floor(numel(angleSet)./2);
    else
        widthIndex = sarcIndex-floor(numel(angleSet)./2);
    end
    
    if lengthIndex(widthIndex) > 0
        metrics(iCell).SarcomereWidth = lengthIndex(widthIndex);
    else
        metrics(iCell).SarcomereWidth = 0;
    end
    metrics(iCell).SarcomereWidth_um = metrics(iCell).SarcomereWidth.*d1_conversion;
    
    % Contrast
    amplitude = zeros(size(contrastdata,1),1);
    lengthIndex = zeros(size(contrastdata,1),1);
    
    for iAngle = 1:size(contrastdata,1)
        [pks,locs,~,p] = findpeaks(contrastdata(iAngle,:));
        
        if ~isempty(pks)
            [pMax,locMax] = max(p);
            amplitude(iAngle) = pMax;
            lengthIndex(iAngle) = spatialSetInterp(locs(locMax));
        end
    end
    
    sarcOrgScore = max(amplitude);
    
    metrics(iCell).SarcomereOrganizationScoreContrast = sarcOrgScore;
    
    % Variance
    metrics(iCell).SarcomereOrganizationScoreVariance = mean(contrastdata(:,11));
    
    % Energy
    amplitude = zeros(size(energydata,1),1);
    lengthIndex = zeros(size(energydata,1),1);
    
    for iAngle = 1:size(energydata,1)
        [pks,locs,~,p] = findpeaks(energydata(iAngle,:));
        
        if ~isempty(pks)
            [pMax,locMax] = max(p);
            amplitude(iAngle) = pMax;
            lengthIndex(iAngle) = spatialSetInterp(locs(locMax));
        end
    end
    
    sarcOrgScore = max(amplitude);
    
    metrics(iCell).SarcomereOrganizationScoreEnergy = sarcOrgScore;
    
    % Homogeneity
    amplitude = zeros(size(homogeneitydata,1),1);
    lengthIndex = zeros(size(homogeneitydata,1),1);
    
    for iAngle = 1:size(homogeneitydata,1)
        [pks,locs,~,p] = findpeaks(homogeneitydata(iAngle,:));
        
        if ~isempty(pks)
            [pMax,locMax] = max(p);
            amplitude(iAngle) = pMax;
            lengthIndex(iAngle) = spatialSetInterp(locs(locMax));
        end
    end
    
    sarcOrgScore = max(amplitude);
    
    metrics(iCell).SarcomereOrganizationScoreHomogeneity = sarcOrgScore;
    
    
    % Gabor
    g = gabor(2:30,angleSet);
    xf = 2:30;
    mag = imgaborfilt(iaact,g);
    mags = squeeze(sum(sum(mag)))./sum(iMask(:));
    mags = reshape(mags,numel(xf),numel(angleSet));
    func = @(v) v(1).*xf.^2   +   v(2).*exp((-(xf-v(3)).^2)./v(4));
    mags_peaks = zeros(1,size(mags,2));
    for i = 1:size(mags,2)
        fitfunc = @(v) func(v) - mags(:,i)';
        x0 = [20 4000 10 20];
        lb = [0 0 0 0];
        ub = [inf inf 20 100];
        options = optimoptions('lsqnonlin','MaxFunctionEvaluations',1000);
        f = lsqnonlin(fitfunc,x0,lb,ub,options);
        mags_peaks(i) = f(2);
%         mags(:,i) =  mags(:,i) - flipud(undercut(flipud(mags(:,i))));
%         funcexp =  @(v) v(1).*xf.^2;
%         funcgauss = @(v) v(2).*exp((-(xf-v(3)).^2)./v(4));
%         figure; hold on;
%         scatter(xf,mags(:,i));
%         plot(xf,func(f));
%         plot(xf,funcexp(f));
%         plot(xf,funcgauss(f));
    end
    gaborScore = max(mags_peaks);
    
    metrics(iCell).GaborScore = gaborScore;
    
    % Fourier
    ft = fft2(iaact);
    ft = fftshift(ft);
    ft = abs(ft);
    ft = imgaussfilt(ft,2);
    
    pf = radialProfile(ft);
    smoothpf = smooth(pf);
    smoothpf = smoothpf./max(smoothpf(:));
%     aperiodic = undercut(smoothpf);
%     periodic = smoothpf - aperiodic;
%     periodic(periodic < 0) = 0;
%     f = fit([0:(numel(periodic)-1)]',periodic,'gauss1');
    xf = 0:(numel(smoothpf)-1);
%     yf = f(xf);
    
%     metrics(iCell).FourierScore = trapz(xf,yf);
    metrics(iCell).RadialProfile = smoothpf;
    metrics(iCell).FourierTransform = ft;
    
    
    func = @(v) v(1).*exp(-xf./v(2))   +   v(3).*exp(-xf./v(4))   +   v(5).*exp((-(xf-v(6)).^2)./v(7));
    fitfunc = @(v) func(v) - smoothpf';
    
    funcexp1 = @(v) v(1).*exp(-xf./v(2));
    funcexp2 = @(v) v(3).*exp(-xf./v(4));
    funcgauss = @(v) v(5).*exp((-(xf-v(6)).^2)./v(7));
    
%         Exp1    Exp2      Gauss
    x0 = [0.9 4   0.08 20   0.1 17 40];
    options = optimoptions('lsqnonlin','MaxFunctionEvaluations',1000);
    f = lsqnonlin(fitfunc,x0,zeros(size(x0)),[inf inf inf inf inf 30 50],options);
    
    metrics(iCell).FitParameters = f;
    metrics(iCell).FitExp1 = funcexp1(f);
    metrics(iCell).FitExp2 = funcexp2(f);
    metrics(iCell).FitGauss = funcgauss(f);
    
    metrics.FourierScore = trapz(xf,funcgauss(f));
    
end


metrics = orderfields(metrics);
end