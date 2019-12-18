%% Master
% This script runs the segmentation, then computes cell morphology metrics
% and sarcomere organization metrics, and finally plots the data as as
% shown in the figures.

clear; close all; clc;

%% Texture analysis parameters
angleSet = 0:1:179;
spatialSet = 0:1:30;
d1_conversion = 160/1024;
d2_conversion = (160/1024).^2;

%% Endogenous CM
eCM_imageSourceFolder = '.\Images\Endogenous CM\';
eCM_imageWriteFolder = '.\Images\Endogenous CM Segmented\';
if ~exist(eCM_imageWriteFolder,'dir')
    mkdir(eCM_imageWriteFolder);
end
eCM_images = dir([eCM_imageSourceFolder '*C1.tif']);
eCM_images = {eCM_images.name};
eCM_images = cellfun(@(x) x(1:end-6),eCM_images,'UniformOutput',false);


[eCM_aactnppa_metrics,eCM_aacthcn4_metrics] = deal(cell(numel(eCM_images),1));
for iImage = 1:numel(eCM_images)
    % Read in four channels + merged image
    iDAPI = imread([eCM_imageSourceFolder eCM_images{iImage} 'C1.tif']);
    iNPPA = imread([eCM_imageSourceFolder eCM_images{iImage} 'C2.tif']);
    iHCN4 = imread([eCM_imageSourceFolder eCM_images{iImage} 'C3.tif']);
    iAACT = imread([eCM_imageSourceFolder eCM_images{iImage} 'C4.tif']);
    iMERGE = imread([eCM_imageSourceFolder eCM_images{iImage} 'Merge.tif']);
    
    % Convert to double format (0-1)
    iDAPI = im2double(iDAPI);
    iNPPA = im2double(iNPPA);
    iHCN4 = im2double(iHCN4);
    iAACT = im2double(iAACT);
    iMERGE = im2double(iMERGE);
    
    disp(['Analyzing endogenous image number: ' num2str(iImage) ' of ' num2str(numel(eCM_images))]);
    
    % Send images to segment function
    iSegmentedImage = eCM_segment_function(iDAPI,iNPPA,iHCN4,iAACT);

    % Write outlined image to file
    iR = iMERGE(:,:,1);
    iG = iMERGE(:,:,2);
    iB = iMERGE(:,:,3);
    
    iOutline_inner = imerode(iSegmentedImage.cell_label,strel('disk',2));
    iOutline_outer = imdilate(iSegmentedImage.cell_label,strel('disk',2));
    iOutline = (iOutline_outer > 0) & ~(iOutline_inner > 0);
    
    iR(iOutline) = 1;
    iG(iOutline) = 1;
    iB(iOutline) = 1;
    
    iMerge_outline = cat(3,iR,iG,iB);
    imwrite(iMerge_outline,[eCM_imageWriteFolder eCM_images{iImage} 'Outline.tif']);
    
    % Send images to morph & texture function
    eCM_aactnppa_metrics{iImage} = morph_texture_function(iAACT,iSegmentedImage.aactnppa_label,...
        angleSet,spatialSet,d1_conversion,d2_conversion,eCM_images{iImage});
    
    eCM_aacthcn4_metrics{iImage} = morph_texture_function(iAACT,iSegmentedImage.aacthcn4_label,...
        angleSet,spatialSet,d1_conversion,d2_conversion,eCM_images{iImage});
end

% Reorganize data into structure
eCM = struct();
eCM.aactnppa = vertcat(eCM_aactnppa_metrics{cellfun(@(x) ~isempty(x),eCM_aactnppa_metrics)});
eCM.aacthcn4 = vertcat(eCM_aacthcn4_metrics{cellfun(@(x) ~isempty(x),eCM_aacthcn4_metrics)});
eCM.all = [eCM.aactnppa; eCM.aacthcn4];
eCM = orderfields(eCM);

clear i* eCM_*_metrics

%% Reprogrammed CM
rCM_imageSourceFolder = '.\Images\Reprogrammed CM\';
rCM_imageWriteFolder = '.\Images\Reprogrammed CM Segmented\';
if ~exist(rCM_imageWriteFolder,'dir')
    mkdir(rCM_imageWriteFolder);
end
rCM_images = dir([rCM_imageSourceFolder '*C1.tif']);
rCM_images = {rCM_images.name};
rCM_images = cellfun(@(x) x(1:end-6),rCM_images,'UniformOutput',false);


[rCM_aact_metrics,rCM_aactnppa_metrics,rCM_hcn4_metrics,rCM_hcn4nppa_metrics,rCM_aacthcn4_metrics,rCM_aacthcn4nppa_metrics] = deal(cell(numel(rCM_images),1));
for iImage = 1:numel(rCM_images)
    % Read in four channels + merged image
    iDAPI = imread([rCM_imageSourceFolder rCM_images{iImage} 'C1.tif']);
    iNPPA = imread([rCM_imageSourceFolder rCM_images{iImage} 'C2.tif']);
    iHCN4 = imread([rCM_imageSourceFolder rCM_images{iImage} 'C3.tif']);
    iAACT = imread([rCM_imageSourceFolder rCM_images{iImage} 'C4.tif']);
    iMERGE = imread([rCM_imageSourceFolder rCM_images{iImage} 'Merge.tif']);
    
    % Convert to double format (0-1)
    iDAPI = im2double(iDAPI);
    iNPPA = im2double(iNPPA);
    iHCN4 = im2double(iHCN4);
    iAACT = im2double(iAACT);
    iMERGE = im2double(iMERGE);
    
    disp(['Analyzing reprogrammed image number: ' num2str(iImage) ' of ' num2str(numel(rCM_images))]);
    
    % Send images to segment function
    iSegmentedImage = rCM_segment_function(iDAPI,iNPPA,iHCN4,iAACT);
    
    iSegmentedImage.aact = iAACT;
    iSegmentedImage.rootname = rCM_images{iImage};
    
    % Write outlined image to file
    iR = iMERGE(:,:,1);
    iG = iMERGE(:,:,2);
    iB = iMERGE(:,:,3);
    
    iOutline_inner = imerode(iSegmentedImage.cell_label,strel('disk',2));
    iOutline_outer = imdilate(iSegmentedImage.cell_label,strel('disk',2));
    iOutline = (iOutline_outer > 0) & ~(iOutline_inner > 0);
    
    iR(iOutline) = 1;
    iG(iOutline) = 1;
    iB(iOutline) = 1;
    
    iMerge_outline = cat(3,iR,iG,iB);
    imwrite(iMerge_outline,[rCM_imageWriteFolder rCM_images{iImage} 'Outline.tif']);
    
    % Send images to morph & texture function
    rCM_aact_metrics{iImage} = morph_texture_function(iAACT,iSegmentedImage.aact_label,...
        angleSet,spatialSet,d1_conversion,d2_conversion,rCM_images{iImage});
    
    rCM_aactnppa_metrics{iImage} = morph_texture_function(iAACT,iSegmentedImage.aactnppa_label,...
        angleSet,spatialSet,d1_conversion,d2_conversion,rCM_images{iImage});
    
    rCM_hcn4_metrics{iImage} = morph_texture_function(iAACT,iSegmentedImage.hcn4_label,...
        angleSet,spatialSet,d1_conversion,d2_conversion,rCM_images{iImage});
    
    rCM_hcn4nppa_metrics{iImage} = morph_texture_function(iAACT,iSegmentedImage.hcn4nppa_label,...
        angleSet,spatialSet,d1_conversion,d2_conversion,rCM_images{iImage});
    
    rCM_aacthcn4_metrics{iImage} = morph_texture_function(iAACT,iSegmentedImage.aacthcn4_label,...
        angleSet,spatialSet,d1_conversion,d2_conversion,rCM_images{iImage});
    
    rCM_aacthcn4nppa_metrics{iImage} = morph_texture_function(iAACT,iSegmentedImage.aacthcn4nppa_label,...
        angleSet,spatialSet,d1_conversion,d2_conversion,rCM_images{iImage});
end

% Reorganize data into structure
rCM = struct();
rCM.aact = vertcat(rCM_aact_metrics{cellfun(@(x) ~isempty(x),rCM_aact_metrics)});
rCM.aactnppa = vertcat(rCM_aactnppa_metrics{cellfun(@(x) ~isempty(x),rCM_aactnppa_metrics)});
rCM.hcn4 = vertcat(rCM_hcn4_metrics{cellfun(@(x) ~isempty(x),rCM_hcn4_metrics)});
rCM.hcn4nppa = vertcat(rCM_hcn4nppa_metrics{cellfun(@(x) ~isempty(x),rCM_hcn4nppa_metrics)});
rCM.aacthcn4 = vertcat(rCM_aacthcn4_metrics{cellfun(@(x) ~isempty(x),rCM_aacthcn4_metrics)});
rCM.aacthcn4nppa = vertcat(rCM_aacthcn4nppa_metrics{cellfun(@(x) ~isempty(x),rCM_aacthcn4nppa_metrics)});
rCM.all = [rCM.aact; rCM.aactnppa; rCM.hcn4; rCM.hcn4nppa; rCM.aacthcn4; rCM.aacthcn4nppa];
rCM.all_aact = [rCM.aact; rCM.aactnppa; rCM.aacthcn4; rCM.aacthcn4nppa];
rCM = orderfields(rCM);

clear i* rCM_*_metrics

%% Generate Beeswarm Plots
close all;
generatePlots;