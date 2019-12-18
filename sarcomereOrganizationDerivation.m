load data.mat;

%%
angleSet = 0:1:179;
spatialSet = 0:1:40;

% organized_aact = eCM.all(5).AACT_BoundingBox;
% organized_mask = eCM.all(5).Image;
% organized_rootname = eCM.all(5).RootName;
%%
organized_rootname = 'Organized 29';
organized_aact = imread('./Images/Metric Comparison/Organized 29.tif');
organized_aact = im2double(organized_aact(:,:,1));
organized_mask = imread('./Images/Metric Comparison/Organized 29 Mask.tif');
organized_mask = organized_mask(:,:,1);

organized_metrics = morph_texture_function(organized_aact,organized_mask,angleSet,spatialSet,d1_conversion,d2_conversion,organized_rootname);

% disorganized_aact = eCM.all(104).AACT_BoundingBox;
% disorganized_mask = eCM.all(104).Image;
% disorganized_rootname = eCM.all(104).RootName;
%%
disorganized_rootname = 'Disorganized 04';
disorganized_aact = imread('./Images/Metric Comparison/Disorganized 04.tif');
disorganized_aact = im2double(disorganized_aact(:,:,1));
disorganized_mask = imread('./Images/Metric Comparison/Disorganized 04 Mask.tif');
disorganized_mask = disorganized_mask(:,:,1);

disorganized_metrics = morph_texture_function(disorganized_aact,disorganized_mask,angleSet,spatialSet,d1_conversion,d2_conversion,disorganized_rootname);


%%
plotCorrelation(organized_metrics,d1_conversion,1)

%%
plotCorrelation(disorganized_metrics,d1_conversion,0);

%% Polar plot
organized_amplitude = zeros(size(organized_metrics.CorrelationData,1),1);
for iAngle = 1:size(organized_metrics.CorrelationData,1)
    [pks,locs,~,p] = findpeaks(organized_metrics.CorrelationData(iAngle,:));
    
    if ~isempty(pks)
        [pMax,locMax] = max(p);
        organized_amplitude(iAngle) = pMax;
    end
end
plotPolar(organized_amplitude);
%%
disorganized_amplitude = zeros(size(disorganized_metrics.CorrelationData,1),1);
for iAngle = 1:size(disorganized_metrics.CorrelationData,1)
    [pks,locs,~,p] = findpeaks(disorganized_metrics.CorrelationData(iAngle,:));
    
    if ~isempty(pks)
        [pMax,locMax] = max(p);
        disorganized_amplitude(iAngle) = pMax;
    end
end
plotPolar(disorganized_amplitude);