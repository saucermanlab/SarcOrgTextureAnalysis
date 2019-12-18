function segmentimages = eCM_segment_function(DAPI,NPPA,HCN4,AACT)
%% FUNCTION DESCRIPTION
% This function segments endogenous cells from the four channels.
%
% INPUTS
% DAPI : raw DAPI (Dapi, nuclear stain) channel image
% NPPA : raw NPPA (Nppa) channel image
% HCN4 : raw HCN4 (Hcn-4) channel image
% AACT : raw AACT (alpha-actinin) channel image
%
% OUTPUTS
% segmentimages : segmentation image with same dimensions as input
%   image. Each pixel is labeled as belonging to a cell
%   (1 = cell #1, 2 = cell #2, etc.) or background = 0.

% Convert images to double format [0 1]
dapi = im2double(DAPI);
nppa = im2double(NPPA);
hcn4 = im2double(HCN4);
aact = im2double(AACT);

% Correct Hcn4 and Nppa bleed-through into Dapi
dapi_correct = dapi;
dapi_correct((dapi_correct - hcn4*0.5) < 0) = 0;
dapi_correct((dapi_correct - nppa*0.5) < 0) = 0;

% Close and filter Dapi
dapi_close = imclose(dapi_correct,strel('disk',4));
dapi_gauss = imgaussfilt(dapi_close,4);

% Apply otsu threshold to filtered Dapi
otsu_dapi = graythresh(dapi_gauss);
dapi_bin = dapi_gauss > otsu_dapi/4;
dapi_fillholes = imfill(dapi_bin,'holes');

% Remove small objects from Dapi
dapi_removesmallobjects = bwareaopen(dapi_fillholes,1500);
dapi_distancetransform = bwdist(~dapi_removesmallobjects);
dapi_distancetransform = -dapi_distancetransform;
dapi_distancetransform(~dapi_removesmallobjects) = -Inf;
dapi_distancetransform = imhmin(dapi_distancetransform,2);
dapi_label = watershed(dapi_distancetransform);
dapi_label(~dapi_removesmallobjects) = 0;
dapi_label = dapi_label > 0;
dapi_label = bwareaopen(dapi_label,800);
dapi_label = bwlabel(dapi_label);
nuclei_label = dapi_label;

% Filter Hcn4
hcn4_gauss = imgaussfilt(hcn4);

% Apply otsu threshold to filtered Hcn4
otsu_hcn4 = graythresh(hcn4_gauss);
hcn4_bin = hcn4 > otsu_hcn4/4;
hcn4_bin_fillholes = imfill(hcn4_bin,'holes');
hcn4_holes = hcn4_bin_fillholes & ~hcn4_bin;
hcn4_bigholes = bwareaopen(hcn4_holes,800);
hcn4_smallholes = hcn4_holes & ~hcn4_bigholes;
hcn4_bin_fillsmallholes = hcn4_bin | hcn4_smallholes;
hcn4_bin_fillsmallholes_dilate = imdilate(hcn4_bin_fillsmallholes,strel('disk',8));
hcn4_bin_smooth = imerode(hcn4_bin_fillsmallholes_dilate,strel('disk',10));

% Classify Hcn4 nuclei. Nuclei fully within the Hcn4 region are considered positive
hcn4_dapi_label = dapi_label;
hcn4_dapi_label_dilate = imdilate(hcn4_dapi_label,strel('disk',6));
for iN = 1:max(hcn4_dapi_label_dilate(:))
    hcn4_nucleus = hcn4_dapi_label_dilate .* hcn4_bin_smooth;
    tCase = hcn4_nucleus(hcn4_dapi_label_dilate == iN);
    if sum(tCase(:)) ~= sum(hcn4_dapi_label_dilate(hcn4_dapi_label_dilate == iN))
        % Outside Hcn4 region, delete
        hcn4_dapi_label(hcn4_dapi_label == iN) = 0;
    end
    if sum(tCase(:)) ~= sum(hcn4_dapi_label_dilate(hcn4_dapi_label_dilate == iN)) && sum(tCase(:)) ~= 0
        dapi_label(dapi_label == iN) = 0;
    end
end

% Filter a-actinin
aact_gauss = imgaussfilt(aact,4);

% Subtract Hcn4 to get only a-actinin region
aact_minusHcn4 = aact_gauss - hcn4;

% Apply otsu threshold to filtered a-actinin
otsu_aact = graythresh(aact_minusHcn4);
aact_bin = aact_minusHcn4 > otsu_aact/2;
aact_bin = imdilate(aact_bin,strel('disk',4));
aact_bin_fillholes = imfill(aact_bin,'holes');
aact_holes = aact_bin_fillholes & ~aact_bin;
aact_bigholes = bwareaopen(aact_holes,4000);
aact_smallholes = aact_holes & ~aact_bigholes;
aact_bin_fillsmallholes = aact_bin | aact_smallholes;
aact_bin_fillsmallholes_dilate = imdilate(aact_bin_fillsmallholes,strel('disk',8));
aact_bin_smooth = imerode(aact_bin_fillsmallholes_dilate,strel('disk',10));

% Classify a-actinin nuclei. Nuclei fully within a-actinin region are considered positive
aact_dapi_label = dapi_label;
aact_dapi_label_dilate = imdilate(dapi_label,strel('disk',6));
for iN = 1:max(aact_dapi_label_dilate(:))
    aact_nucleus = aact_dapi_label_dilate .* aact_bin_smooth;
    tCase = aact_nucleus(aact_dapi_label_dilate == iN);
    if sum(tCase(:)) ~= sum(aact_dapi_label_dilate(aact_dapi_label_dilate == iN))
        % Outside a-actinin region, delete
        aact_dapi_label(aact_dapi_label == iN) = 0;
    end
end

% Create perinuclear ring for Nppa classification
dapi_label_dilate = imdilate(dapi_label,strel('disk',8));
dapi_label_erode = imerode(dapi_label,strel('disk',2));
dapi_label_ring = dapi_label_dilate - dapi_label_erode;

% Clasify Nppa nuclei. Only a-actinin nuclei are considered for Nppa classification.
aactnppa_dapi_label = aact_dapi_label;
for iN = 1:max(aactnppa_dapi_label(:))
    if sum(aactnppa_dapi_label(:) == iN) == 0
        continue;
    end
    nppaRing = nppa(dapi_label_ring == iN);
    if prctile(nppaRing,90) < 0.1
        aactnppa_dapi_label(aactnppa_dapi_label == iN) = 0;
    end
end

% Identify binucleates in the Hcn4 image
[hcn4_boundaries, hcn4_dapi_label_binucleates] = bwboundaries(hcn4_dapi_label);
hcn4_nucleiToJoin = [];
hcn4_coordinates = [];
if numel(hcn4_boundaries) > 1
    for b1 = 1:numel(hcn4_boundaries)
        for b2 = (b1+1):numel(hcn4_boundaries)
            b1pts = hcn4_boundaries{b1};
            b2pts = hcn4_boundaries{b2};
            
            rb1pts = repelem(b1pts,size(b2pts,1),1);
            rb2pts = repmat(b2pts,size(b1pts,1),1);
            dists = sqrt((rb1pts(:,1)-rb2pts(:,1)).^2+(rb1pts(:,2)-rb2pts(:,2)).^2);
            
            [minDist,minDistIndex] = min(dists);
            if minDist < 25
                hcn4_nucleiToJoin = [hcn4_nucleiToJoin; b1 b2];
                hcn4_coordinates = [hcn4_coordinates; [rb1pts(minDistIndex,:) rb2pts(minDistIndex,:)]];
            end
        end
    end
end

% Combine binucleates in the Hcn4 image
for iJoin = 1:size(hcn4_nucleiToJoin,1)
    iRow = hcn4_nucleiToJoin(end+1-iJoin,:);
    hcn4_dapi_label_binucleates(hcn4_dapi_label_binucleates == iRow(2)) = iRow(1);
    nPoints = ceil(sqrt((hcn4_coordinates(iJoin,1) - hcn4_coordinates(iJoin,3)).^2 + (hcn4_coordinates(iJoin,2) - hcn4_coordinates(iJoin,4)).^2)) + 1;
    xvalues = round(linspace(hcn4_coordinates(iJoin,1),hcn4_coordinates(iJoin,3),nPoints));
    yvalues = round(linspace(hcn4_coordinates(iJoin,2),hcn4_coordinates(iJoin,4),nPoints));
    hcn4_dapi_label_binucleates(sub2ind(size(hcn4_dapi_label_binucleates), xvalues, yvalues)) = iRow(1);
end

% Combine a-actinin and Hcn4 for cell segmentation
aact_gauss_gauss = imgaussfilt(aact_gauss,4);
hcn4_gauss = imgaussfilt(hcn4,4);
aacthcn4 = aact_gauss_gauss.*hcn4_gauss.*2;
aacthcn4_combinemask = aacthcn4*0.5 + 0.5*hcn4_bin_smooth;

% Calculate gradient and perform watershed transform to segment Hcn4 cells
hy = fspecial('sobel');
hx = hy';
Iy = imfilter(aacthcn4_combinemask, hy, 'replicate');
Ix = imfilter(aacthcn4_combinemask, hx, 'replicate');
hcn4_gradmag = sqrt(Ix.^2 + Iy.^2);

hcn4_distancetransform = bwdist(hcn4_bin_smooth);
hcn4_label = watershed(hcn4_distancetransform);
hcn4_backgroundmarker = hcn4_label == 0;

hcn4_gradmag2 = imimposemin(hcn4_gradmag,hcn4_backgroundmarker | hcn4_dapi_label_binucleates);
hcn4_label = watershed(hcn4_gradmag2);
hcn4_label(~hcn4_bin_smooth) = 0;
for iLabel = 1:max(hcn4_label(:))
    label_i = hcn4_label == iLabel;
    if sum(label_i(:).*(hcn4_dapi_label_binucleates(:) > 0)) == 0
        hcn4_label(hcn4_label == iLabel) = 0;
    end
end

% Fill holes on objects
for iLabel = 1:max(hcn4_label(:))
    label_i = hcn4_label == iLabel;
    label_i_fillholes = imfill(label_i,'holes');
    hcn4_label(label_i_fillholes) = iLabel;
end

% Remove objects touching edge of image
border = zeros(size(DAPI));
border(1,:) = 1; border(end,:) = 1;
border(:,1) = 1; border(:,end) = 1;
for iLabel = 1:max(hcn4_label(:))
    label_i = hcn4_label == iLabel;
    label_i_dilate = imdilate(label_i,strel('disk',1));
    label_i_ring = label_i_dilate - label_i;
    label_i_ring_border = label_i_ring & border;
    if sum(label_i_ring_border(:))/sum(label_i_ring(:)) > 0.1
        hcn4_label(label_i) = 0;
    end
end

% Create mask from Hcn4 cells
hcn4_cell_bin = hcn4_label > 0;
hcn4_cell_bin_dilate = imdilate(hcn4_cell_bin,strel('disk',4));
aact_bin_maskhcn4 = aact_bin_smooth;
aact_bin_maskhcn4(hcn4_cell_bin_dilate) = 0;
aact_combinemask = aact_gauss.*0.5 + 0.5.*aact_bin_maskhcn4;

% Identify binucleates in the a-actinin image
[aactnppa_boundaries, aactnppa_dapi_label_binucleates] = bwboundaries(aactnppa_dapi_label);
aactnppa_nucleiToJoin = [];
aactnppa_coordinates = [];
if numel(aactnppa_boundaries) > 1
    for b1 = 1:numel(aactnppa_boundaries)
        for b2 = (b1+1):numel(aactnppa_boundaries)
            b1pts = aactnppa_boundaries{b1};
            b2pts = aactnppa_boundaries{b2};
            
            rb1pts = repelem(b1pts,size(b2pts,1),1);
            rb2pts = repmat(b2pts,size(b1pts,1),1);
            dists = sqrt((rb1pts(:,1)-rb2pts(:,1)).^2+(rb1pts(:,2)-rb2pts(:,2)).^2);
            
            [minDist,minDistIndex] = min(dists);
            if minDist < 25
                aactnppa_nucleiToJoin = [aactnppa_nucleiToJoin; b1 b2];
                aactnppa_coordinates = [aactnppa_coordinates; [rb1pts(minDistIndex,:) rb2pts(minDistIndex,:)]];
            end
        end
    end
end

% Combine binucleates in the a-actinin image
for iJoin = 1:size(aactnppa_nucleiToJoin,1)
    iRow = aactnppa_nucleiToJoin(end+1-iJoin,:);
    aactnppa_dapi_label_binucleates(aactnppa_dapi_label_binucleates == iRow(2)) = iRow(1);
    nPoints = ceil(sqrt((aactnppa_coordinates(iJoin,1) - aactnppa_coordinates(iJoin,3)).^2 + (aactnppa_coordinates(iJoin,2) - aactnppa_coordinates(iJoin,4)).^2)) + 1;
    xvalues = round(linspace(aactnppa_coordinates(iJoin,1),aactnppa_coordinates(iJoin,3),nPoints));
    yvalues = round(linspace(aactnppa_coordinates(iJoin,2),aactnppa_coordinates(iJoin,4),nPoints));
    aactnppa_dapi_label_binucleates(sub2ind(size(aactnppa_dapi_label_binucleates), xvalues, yvalues)) = iRow(1);
end

% Calculate gradient and perform watershed transform to segment a-actinin cells
hy = fspecial('sobel');
hx = hy';
Iy = imfilter(aact_combinemask, hy, 'replicate');
Ix = imfilter(aact_combinemask, hx, 'replicate');
aact_gradmag = sqrt(Ix.^2 + Iy.^2);

aact_distancetransform = bwdist(aact_bin_smooth);
aact_label = watershed(aact_distancetransform);
aact_backgroundmarker = aact_label == 0;

aact_gradmag2 = imimposemin(aact_gradmag,aact_backgroundmarker | aactnppa_dapi_label_binucleates);
aactnppa_label = watershed(aact_gradmag2);
aactnppa_label(~aact_bin_smooth) = 0;
for iLabel = 1:max(aactnppa_label(:))
    label_i = aactnppa_label == iLabel;
    if sum(label_i(:).*(aactnppa_dapi_label_binucleates(:) > 0)) == 0
        aactnppa_label(aactnppa_label == iLabel) = 0;
    end
end

% Fill holes in objects
for iLabel = 1:max(aactnppa_label(:))
    label_i = aactnppa_label == iLabel;
    label_i_fillholes = imfill(label_i,'holes');
    aactnppa_label(label_i_fillholes) = iLabel;
end

% Remove a-actinin objects touching edge of image
border = zeros(size(DAPI));
border(1,:) = 1; border(end,:) = 1;
border(:,1) = 1; border(:,end) = 1;
for iLabel = 1:max(aactnppa_label(:))
    label_i = aactnppa_label == iLabel;
    label_i_erode = imerode(label_i,strel('disk',1));
    label_i_ring = label_i - label_i_erode;
    label_i_border = label_i & border;
    label_i_perim = label_i_border | label_i_ring;
    if sum(label_i_border(:))/sum(label_i_perim(:)) > 0.2
        aactnppa_label(label_i) = 0;
    end
end

% Remove Hcn4 objects touching edge of image
border = zeros(size(DAPI));
border(1,:) = 1; border(end,:) = 1;
border(:,1) = 1; border(:,end) = 1;
for iLabel = 1:max(hcn4_label(:))
    label_i = hcn4_label == iLabel;
    label_i_erode = imerode(label_i,strel('disk',1));
    label_i_ring = label_i - label_i_erode;
    label_i_border = label_i & border;
    label_i_perim = label_i_border | label_i_ring;
    if sum(label_i_border(:))/sum(label_i_perim(:)) > 0.2
        hcn4_label(label_i) = 0;
    end
end

% Reassign Hcn4 object numbers
hcn4_dapi_label_binucleates_bin = hcn4_dapi_label_binucleates > 0;
hcn4_nuclei_label = bwlabel(hcn4_dapi_label_binucleates_bin);

% Reassign a-actinin object numbers
aactnppa_dapi_label_binucleates_bin = aactnppa_dapi_label_binucleates > 0;
aactnppa_nuclei_label = bwlabel(aactnppa_dapi_label_binucleates_bin);

%% Reassign object numbers. Nuclei and cells will have same object number.
aact_label_bin = aactnppa_label > 0;
aact_label_bin = bwareaopen(aact_label_bin,10000);
aact_nuclei_label_bin = (aactnppa_nuclei_label > 0) & aact_label_bin;
aact_label = bwlabel(aact_label_bin);
aact_nuclei_label = aact_label;
aact_nuclei_label(~aact_nuclei_label_bin) = 0;

hcn4_label_bin = hcn4_label > 0;
hcn4_label_bin = bwareaopen(hcn4_label_bin,10000);
hcn4_nuclei_label_bin = (hcn4_nuclei_label > 0) & hcn4_label_bin;
hcn4_label = bwlabel(hcn4_label_bin);
hcn4_nuclei_label = hcn4_label;
hcn4_nuclei_label(~hcn4_nuclei_label_bin) = 0;

cell_label = aact_label;
cell_label = cell_label + max(cell_label(:)).*hcn4_label_bin + hcn4_label;

segmentimages = struct();

segmentimages.nuclei_label = nuclei_label;
segmentimages.cell_label = cell_label;

segmentimages.aactnppa_nuclei_label = aact_nuclei_label;
segmentimages.aactnppa_label = aact_label;

segmentimages.aacthcn4_nuclei_label = hcn4_nuclei_label;
segmentimages.aacthcn4_label = hcn4_label;
end