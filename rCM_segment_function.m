function segmentimages = rCM_segment_function(DAPI,NPPA,HCN4,AACT)
%% FUNCTION DESCRIPTION
% This function segments reprogrammed cells from the four channels.
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
%%
% Convert images to double format [0 1]
dapi = im2double(DAPI);
nppa = im2double(NPPA);
hcn4 = im2double(HCN4);
aact = im2double(AACT);

% Close and filter Dapi
dapi_close = imclose(dapi,strel('disk',4));
dapi_gauss = imgaussfilt(dapi_close,4);

% Apply otsu threshold to filtered Dapi
otsu_dapi = graythresh(dapi_gauss);
dapi_bin = dapi_gauss > otsu_dapi/2;
dapi_fillholes = imfill(dapi_bin,'holes');

% Remove small objects from Dapi
dapi_removesmallobjects = bwareaopen(dapi_fillholes,500);
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

% Filter a-actinin
aact_gauss = imgaussfilt(aact,4);

% Apply otsu threshold to filtered a-actinin
otsu_aact = graythresh(aact_gauss);
aact_bin = aact_gauss > otsu_aact/4;
aact_bin = imdilate(aact_bin,strel('disk',4));
aact_bin_fillholes = imfill(aact_bin,'holes');
aact_holes = aact_bin_fillholes & ~aact_bin;
aact_bigholes = bwareaopen(aact_holes,4000);
aact_smallholes = aact_holes & ~aact_bigholes;
aact_bin_fillsmallholes = aact_bin | aact_smallholes;
aact_bin_fillsmallholes_dilate = imdilate(aact_bin_fillsmallholes,strel('disk',8));
aact_bin_smooth = imerode(aact_bin_fillsmallholes_dilate,strel('disk',10));

% Create Hcn4 only mask
aact_bin_smooth_dilate = imdilate(aact_bin_smooth,strel('disk',8));
hcn4only_bin_smooth = hcn4_bin_smooth;
hcn4only_bin_smooth(aact_bin_smooth_dilate) = 0;
hcn4only_bin_smooth_fillholes = imfill(hcn4only_bin_smooth,'holes');
hcn4only_holes = hcn4only_bin_smooth_fillholes & ~hcn4_bin_smooth;
hcn4only_bigholes = bwareaopen(hcn4only_holes,10000);
hcn4only_smallholes = hcn4only_holes & ~hcn4only_bigholes;
hcn4only_bin_smooth_fillholes = hcn4only_bin_smooth | hcn4only_smallholes;
hcn4only_bin_smooth_fillholes = bwareaopen(hcn4only_bin_smooth_fillholes,5000);

% Classify Hcn4 nuclei. Nuclei fully within the Hcn4 region are considered positive
hcn4_dapi_label = dapi_label;
hcn4_dapi_label_dilate = imdilate(hcn4_dapi_label,strel('disk',1));
for iN = 1:max(hcn4_dapi_label_dilate(:))
    hcn4_nucleus = hcn4_dapi_label_dilate .* hcn4only_bin_smooth_fillholes;
    tCase = hcn4_nucleus(hcn4_dapi_label_dilate == iN);
    if sum(tCase(:)) ~= sum(hcn4_dapi_label_dilate(hcn4_dapi_label_dilate == iN))
        % Outside Hcn4 region, delete
        hcn4_dapi_label(hcn4_dapi_label == iN) = 0;
    end
    if sum(tCase(:)) ~= sum(hcn4_dapi_label_dilate(hcn4_dapi_label_dilate == iN)) && sum(tCase(:)) ~= 0
        dapi_label(dapi_label == iN) = 0;
    end
end

% Create a-actinin only mask
hcn4_bin_smooth_dilate = imdilate(hcn4_bin_smooth,strel('disk',8));
aactonly_bin_smooth = aact_bin_smooth;
aactonly_bin_smooth(hcn4_bin_smooth_dilate) = 0;
aactonly_bin_smooth_fillholes = imfill(aactonly_bin_smooth,'holes');
aactonly_holes = aactonly_bin_smooth_fillholes & ~aactonly_bin_smooth;
aactonly_bigholes = bwareaopen(aactonly_holes,10000);
aactonly_smallholes = aactonly_holes & ~aactonly_bigholes;
aactonly_bin_smooth_fillholes = aactonly_bin_smooth | aactonly_smallholes;
aactonly_bin_smooth_fillholes = bwareaopen(aactonly_bin_smooth_fillholes,5000);

% Classify a-actinin nuclei
aact_dapi_label = dapi_label;
aact_dapi_label_dilate = imdilate(aact_dapi_label,strel('disk',1));
for iN = 1:max(aact_dapi_label_dilate(:))
    aact_nucleus = aact_dapi_label_dilate .* aactonly_bin_smooth_fillholes;
    tCase = aact_nucleus(aact_dapi_label_dilate == iN);
    if sum(tCase(:)) ~= sum(aact_dapi_label_dilate(aact_dapi_label_dilate == iN))
        % Outside a-actinin region, delete
        aact_dapi_label(aact_dapi_label == iN) = 0;
    end
    if sum(tCase(:)) ~= sum(aact_dapi_label_dilate(aact_dapi_label_dilate == iN)) && sum(tCase(:)) ~= 0
        dapi_label(dapi_label == iN) = 0;
    end
end

% Create a-actinin and Hcn4 mask
aacthcn4_bin_smooth = aact_bin_smooth & hcn4_bin_smooth;
aacthcn4_bin_smooth_fillholes = imfill(aacthcn4_bin_smooth,'holes');
aacthcn4_holes = aacthcn4_bin_smooth_fillholes & ~aacthcn4_bin_smooth;
aacthcn4_bigholes = bwareaopen(aacthcn4_holes,10000);
aacthcn4_smallholes = aacthcn4_holes & ~aacthcn4_bigholes;
aacthcn4_bin_smooth_fillholes = aacthcn4_bin_smooth | aacthcn4_smallholes;
aacthcn4_bin_smooth_fillholes = bwareaopen(aacthcn4_bin_smooth_fillholes,5000);


% Classify a-actinin and Hcn4 nuclei
aacthcn4_dapi_label = dapi_label;
aacthcn4_dapi_label_dilate = imdilate(aacthcn4_dapi_label,strel('disk',1));
for iN = 1:max(aacthcn4_dapi_label_dilate(:))
    aacthcn4_nucleus = aacthcn4_dapi_label_dilate .* aacthcn4_bin_smooth_fillholes;
    tCase = aacthcn4_nucleus(aacthcn4_dapi_label_dilate == iN);
    if sum(tCase(:)) ~= sum(aacthcn4_dapi_label_dilate(aacthcn4_dapi_label_dilate == iN))
        % Outside a-actinin and Hcn4 region, delete
        aacthcn4_dapi_label(aacthcn4_dapi_label == iN) = 0;
    end
    if sum(tCase(:)) ~= sum(aacthcn4_dapi_label_dilate(aacthcn4_dapi_label_dilate == iN)) && sum(tCase(:)) ~= 0
        dapi_label(dapi_label == iN) = 0;
    end
end

%% Reclassify binucleates
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

% Reassign Hcn4 object numbers
hcn4_dapi_label_binucleates_bin = hcn4_dapi_label_binucleates > 0;
hcn4_nuclei_label = bwlabel(hcn4_dapi_label_binucleates_bin);

% Identify binucleates in the a-actinin image
[aact_boundaries, aact_dapi_label_binucleates] = bwboundaries(aact_dapi_label);
aact_nucleiToJoin = [];
aact_coordinates = [];
if numel(aact_boundaries) > 1
    for b1 = 1:numel(aact_boundaries)
        for b2 = (b1+1):numel(aact_boundaries)
            b1pts = aact_boundaries{b1};
            b2pts = aact_boundaries{b2};
            
            rb1pts = repelem(b1pts,size(b2pts,1),1);
            rb2pts = repmat(b2pts,size(b1pts,1),1);
            dists = sqrt((rb1pts(:,1)-rb2pts(:,1)).^2+(rb1pts(:,2)-rb2pts(:,2)).^2);
            
            [minDist,minDistIndex] = min(dists);
            if minDist < 25
                aact_nucleiToJoin = [aact_nucleiToJoin; b1 b2];
                aact_coordinates = [aact_coordinates; [rb1pts(minDistIndex,:) rb2pts(minDistIndex,:)]];
            end
        end
    end
end

% Combine binucleates in the a-actinin image
for iJoin = 1:size(aact_nucleiToJoin,1)
    iRow = aact_nucleiToJoin(end+1-iJoin,:);
    aact_dapi_label_binucleates(aact_dapi_label_binucleates == iRow(2)) = iRow(1);
    nPoints = ceil(sqrt((aact_coordinates(iJoin,1) - aact_coordinates(iJoin,3)).^2 + (aact_coordinates(iJoin,2) - aact_coordinates(iJoin,4)).^2)) + 1;
    xvalues = round(linspace(aact_coordinates(iJoin,1),aact_coordinates(iJoin,3),nPoints));
    yvalues = round(linspace(aact_coordinates(iJoin,2),aact_coordinates(iJoin,4),nPoints));
    aact_dapi_label_binucleates(sub2ind(size(aact_dapi_label_binucleates), xvalues, yvalues)) = iRow(1);
end

% Reassign a-actinin object numbers
aact_dapi_label_binucleates_bin = aact_dapi_label_binucleates > 0;
aact_nuclei_label = bwlabel(aact_dapi_label_binucleates_bin);

% Identify binucleates in the a-actinin and Hcn4 image
[aacthcn4_boundaries, aacthcn4_dapi_label_binucleates] = bwboundaries(aacthcn4_dapi_label);
aacthcn4_nucleiToJoin = [];
aacthcn4_coordinates = [];
if numel(aacthcn4_boundaries) > 1
    for b1 = 1:numel(aacthcn4_boundaries)
        for b2 = (b1+1):numel(aacthcn4_boundaries)
            b1pts = aacthcn4_boundaries{b1};
            b2pts = aacthcn4_boundaries{b2};
            
            rb1pts = repelem(b1pts,size(b2pts,1),1);
            rb2pts = repmat(b2pts,size(b1pts,1),1);
            dists = sqrt((rb1pts(:,1)-rb2pts(:,1)).^2+(rb1pts(:,2)-rb2pts(:,2)).^2);
            
            [minDist,minDistIndex] = min(dists);
            if minDist < 25
                aacthcn4_nucleiToJoin = [aacthcn4_nucleiToJoin; b1 b2];
                aacthcn4_coordinates = [aacthcn4_coordinates; [rb1pts(minDistIndex,:) rb2pts(minDistIndex,:)]];
            end
        end
    end
end

% Combine binucleates in the a-actinin and Hcn4 image
for iJoin = 1:size(aacthcn4_nucleiToJoin,1)
    iRow = aacthcn4_nucleiToJoin(end+1-iJoin,:);
    aacthcn4_dapi_label_binucleates(aacthcn4_dapi_label_binucleates == iRow(2)) = iRow(1);
    nPoints = ceil(sqrt((aacthcn4_coordinates(iJoin,1) - aacthcn4_coordinates(iJoin,3)).^2 + (aacthcn4_coordinates(iJoin,2) - aacthcn4_coordinates(iJoin,4)).^2)) + 1;
    xvalues = round(linspace(aacthcn4_coordinates(iJoin,1),aacthcn4_coordinates(iJoin,3),nPoints));
    yvalues = round(linspace(aacthcn4_coordinates(iJoin,2),aacthcn4_coordinates(iJoin,4),nPoints));
    aacthcn4_dapi_label_binucleates(sub2ind(size(aacthcn4_dapi_label_binucleates), xvalues, yvalues)) = iRow(1);
end

% Reassign a-actinin and Hcn4 object numbers
aacthcn4_dapi_label_binucleates_bin = aacthcn4_dapi_label_binucleates > 0;
aacthcn4_nuclei_label = bwlabel(aacthcn4_dapi_label_binucleates_bin);

%% Watershed segmentation
% Create mask
hcn4_combinemask = hcn4_gauss*0.5 + 0.5*hcn4only_bin_smooth_fillholes;
% Calculate gradient and perform watershed transform to segment Hcn4 cells
hy = fspecial('sobel');
hx = hy';
Iy = imfilter(hcn4_combinemask, hy, 'replicate');
Ix = imfilter(hcn4_combinemask, hx, 'replicate');
hcn4_gradmag = sqrt(Ix.^2 + Iy.^2);

hcn4_distancetransform = bwdist(hcn4_bin_smooth);
hcn4_label = watershed(hcn4_distancetransform);
hcn4_backgroundmarker = hcn4_label == 0;

hcn4_gradmag2 = imimposemin(hcn4_gradmag,hcn4_backgroundmarker | hcn4_nuclei_label);
hcn4_label = watershed(hcn4_gradmag2);
hcn4_label(~hcn4_bin_smooth) = 0;
for iLabel = 1:max(hcn4_label(:))
    label_i = hcn4_label == iLabel;
    if sum(label_i(:).*(hcn4_nuclei_label(:) > 0)) == 0
        hcn4_label(hcn4_label == iLabel) = 0;
    end
end

% Fill holes on objects
for iLabel = 1:max(hcn4_label(:))
    label_i = hcn4_label == iLabel;
    label_i_fillholes = imfill(label_i,'holes');
    label_i_holes = label_i_fillholes & ~label_i;
    label_i_bigholes = bwareaopen(label_i_holes,2000);
    label_i_smallholes = label_i_holes & ~ label_i_bigholes;
    label_i_fillholes = label_i | label_i_smallholes;
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
    label_i_border = label_i & border;
    label_i_perim = label_i_border | label_i_ring;
    if sum(label_i_border(:))/sum(label_i_perim(:)) > 0.2
        hcn4_label(label_i) = 0;
    end
end

% Create mask from Hcn4 cells
hcn4_cell_bin = hcn4_label > 0;
hcn4_cell_bin_dilate = imdilate(hcn4_cell_bin,strel('disk',4));
aactonly_bin_smooth_fillholes_maskhcn4 = aactonly_bin_smooth_fillholes;
aactonly_bin_smooth_fillholes_maskhcn4(hcn4_cell_bin_dilate) = 0;
aact_combinemask = aact_gauss.*0.5 + 0.5.*aactonly_bin_smooth_fillholes_maskhcn4;

% Calculate gradient and perform watershed transform to segment a-actinin cells
hy = fspecial('sobel');
hx = hy';
Iy = imfilter(aact_combinemask, hy, 'replicate');
Ix = imfilter(aact_combinemask, hx, 'replicate');
aact_gradmag = sqrt(Ix.^2 + Iy.^2);

aact_distancetransform = bwdist(aact_bin_smooth);
aact_label = watershed(aact_distancetransform);
aact_backgroundmarker = aact_label == 0;

aact_gradmag2 = imimposemin(aact_gradmag,aact_backgroundmarker | aact_nuclei_label);
aact_label = watershed(aact_gradmag2);
aact_label(~aact_bin_smooth) = 0;
for iLabel = 1:max(aact_label(:))
    label_i = aact_label == iLabel;
    if sum(label_i(:).*(aact_nuclei_label(:) > 0)) == 0
        aact_label(aact_label == iLabel) = 0;
    end
end

% Fill holes on objects
for iLabel = 1:max(aact_label(:))
    label_i = aact_label == iLabel;
    label_i_fillholes = imfill(label_i,'holes');
    label_i_holes = label_i_fillholes & ~label_i;
    label_i_bigholes = bwareaopen(label_i_holes,2000);
    label_i_smallholes = label_i_holes & ~ label_i_bigholes;
    label_i_fillholes = label_i | label_i_smallholes;
    aact_label(label_i_fillholes) = iLabel;
end

% Remove objects touching edge of image
border = zeros(size(DAPI));
border(1,:) = 1; border(end,:) = 1;
border(:,1) = 1; border(:,end) = 1;
for iLabel = 1:max(aact_label(:))
    label_i = aact_label == iLabel;
    label_i_dilate = imdilate(label_i,strel('disk',1));
    label_i_ring = label_i_dilate - label_i;
    label_i_border = label_i & border;
    label_i_perim = label_i_border | label_i_ring;
    if sum(label_i_border(:))/sum(label_i_perim(:)) > 0.2
        aact_label(label_i) = 0;
    end
end

% Create mask from Hcn4 or a-actinin cells
hcn4_cell_bin = hcn4_label > 0;
aact_cell_bin = aact_label > 0;
hcn4_cell_bin_dilate = imdilate(hcn4_cell_bin,strel('disk',4));
aact_cell_bin_dilate = imdilate(aact_cell_bin,strel('disk',4));
aacthcn4_bin_smooth_fillholes_maskcells = aacthcn4_bin_smooth_fillholes;
aacthcn4_bin_smooth_fillholes_maskcells(hcn4_cell_bin_dilate) = 0;
aacthcn4_bin_smooth_fillholes_maskcells(aact_cell_bin_dilate) = 0;
aacthcn4_combinemask = aact_gauss.*0.25 + hcn4_gauss.*0.25 + 0.5.*aacthcn4_bin_smooth_fillholes_maskcells;

% Calculate gradient and perform watershed transform to segment a-actinin cells
hy = fspecial('sobel');
hx = hy';
Iy = imfilter(aacthcn4_combinemask, hy, 'replicate');
Ix = imfilter(aacthcn4_combinemask, hx, 'replicate');
aacthcn4_gradmag = sqrt(Ix.^2 + Iy.^2);

aacthcn4_distancetransform = bwdist(aacthcn4_bin_smooth);
aacthcn4_label = watershed(aacthcn4_distancetransform);
aacthcn4_backgroundmarker = aacthcn4_label == 0;

aacthcn4_gradmag2 = imimposemin(aacthcn4_gradmag,aacthcn4_backgroundmarker | aacthcn4_nuclei_label);
aacthcn4_label = watershed(aacthcn4_gradmag2);
aacthcn4_label(~aacthcn4_bin_smooth) = 0;
for iLabel = 1:max(aacthcn4_label(:))
    label_i = aacthcn4_label == iLabel;
    if sum(label_i(:).*(aacthcn4_nuclei_label(:) > 0)) == 0
        aacthcn4_label(aacthcn4_label == iLabel) = 0;
    end
end

% Fill holes on objects
for iLabel = 1:max(aacthcn4_label(:))
    label_i = aacthcn4_label == iLabel;
    label_i_fillholes = imfill(label_i,'holes');
    label_i_holes = label_i_fillholes & ~label_i;
    label_i_bigholes = bwareaopen(label_i_holes,2000);
    label_i_smallholes = label_i_holes & ~ label_i_bigholes;
    label_i_fillholes = label_i | label_i_smallholes;
    aacthcn4_label(label_i_fillholes) = iLabel;
end

% Remove objects touching edge of image
border = zeros(size(DAPI));
border(1,:) = 1; border(end,:) = 1;
border(:,1) = 1; border(:,end) = 1;
for iLabel = 1:max(aacthcn4_label(:))
    label_i = aacthcn4_label == iLabel;
    label_i_dilate = imdilate(label_i,strel('disk',1));
    label_i_ring = label_i_dilate - label_i;
    label_i_border = label_i & border;
    label_i_perim = label_i_border | label_i_ring;
    if sum(label_i_border(:))/sum(label_i_perim(:)) > 0.2
        aacthcn4_label(label_i) = 0;
    end
end

%% Reassign object numbers. Nuclei and cells will have same object number.
aact_label_bin = aact_label > 0;
aact_nuclei_label_bin = (aact_nuclei_label > 0) & aact_label_bin;
aact_label = bwlabel(aact_label_bin);
aact_nuclei_label = aact_label;
aact_nuclei_label(~aact_nuclei_label_bin) = 0;

hcn4_label_bin = hcn4_label > 0;
hcn4_nuclei_label_bin = (hcn4_nuclei_label > 0) & hcn4_label_bin;
hcn4_label = bwlabel(hcn4_label_bin);
hcn4_nuclei_label = hcn4_label;
hcn4_nuclei_label(~hcn4_nuclei_label_bin) = 0;

aacthcn4_label_bin = aacthcn4_label > 0;
aacthcn4_nuclei_label_bin = (aacthcn4_nuclei_label > 0) & aacthcn4_label_bin;
aacthcn4_label = bwlabel(aacthcn4_label_bin);
aacthcn4_nuclei_label = aacthcn4_label;
aacthcn4_nuclei_label(~aacthcn4_nuclei_label_bin) = 0;

%% Nppa classification
aact_nuclei_label_dilate = imdilate(aact_nuclei_label,strel('disk',4));
aact_nuclei_label_erode = imerode(aact_nuclei_label,strel('disk',2));
aact_nuclei_label_ring = aact_nuclei_label_dilate;
aact_nuclei_label_ring(aact_nuclei_label_erode > 0) = 0;

aactnppa_nuclei_label = aact_nuclei_label;
aactnppa_label = aact_label;
for iN = 1:max(aact_nuclei_label(:))
    nppaRing = nppa(aact_nuclei_label_ring == iN);
    if prctile(nppaRing,90) > 0.1
        aact_nuclei_label(aact_nuclei_label == iN) = 0;
        aact_label(aact_label == iN) = 0;
    else
        aactnppa_nuclei_label(aactnppa_nuclei_label == iN) = 0;
        aactnppa_label(aactnppa_label == iN) = 0;
   end
end

hcn4_nuclei_label_dilate = imdilate(hcn4_nuclei_label,strel('disk',4));
hcn4_nuclei_label_erode = imerode(hcn4_nuclei_label,strel('disk',2));
hcn4_nuclei_label_ring = hcn4_nuclei_label_dilate;
hcn4_nuclei_label_ring(hcn4_nuclei_label_erode > 0) = 0;

hcn4nppa_nuclei_label = hcn4_nuclei_label;
hcn4nppa_label = hcn4_label;
for iN = 1:max(hcn4_nuclei_label(:))
    nppaRing = nppa(hcn4_nuclei_label_ring == iN);
    if prctile(nppaRing,90) > 0.1
        hcn4_nuclei_label(hcn4_nuclei_label == iN) = 0;
        hcn4_label(hcn4_label == iN) = 0;
    else
        hcn4nppa_nuclei_label(hcn4nppa_nuclei_label == iN) = 0;
        hcn4nppa_label(hcn4nppa_label == iN) = 0;
   end
end

aacthcn4_nuclei_label_dilate = imdilate(aacthcn4_nuclei_label,strel('disk',4));
aacthcn4_nuclei_label_erode = imerode(aacthcn4_nuclei_label,strel('disk',2));
aacthcn4_nuclei_label_ring = aacthcn4_nuclei_label_dilate;
aacthcn4_nuclei_label_ring(aacthcn4_nuclei_label_erode > 0) = 0;

aacthcn4nppa_nuclei_label = aacthcn4_nuclei_label;
aacthcn4nppa_label = aacthcn4_label;
for iN = 1:max(aacthcn4_nuclei_label(:))
    nppaRing = nppa(aacthcn4_nuclei_label_ring == iN);
    if prctile(nppaRing,90) > 0.1
        aacthcn4_nuclei_label(aacthcn4_nuclei_label == iN) = 0;
        aacthcn4_label(aacthcn4_label == iN) = 0;
    else
        aacthcn4nppa_nuclei_label(aacthcn4nppa_nuclei_label == iN) = 0;
        aacthcn4nppa_label(aacthcn4nppa_label == iN) = 0;
   end
end

%% Once again, reassign object numbers with Nppa classification
aact_label_bin = aact_label > 0;
aact_label_bin = bwareaopen(aact_label_bin,10000);
aact_nuclei_label_bin = (aact_nuclei_label > 0) & aact_label_bin;
aact_label = bwlabel(aact_label_bin);
aact_nuclei_label = aact_label;
aact_nuclei_label(~aact_nuclei_label_bin) = 0;

aactnppa_label_bin = aactnppa_label > 0;
aactnppa_label_bin = bwareaopen(aactnppa_label_bin,10000);
aactnppa_nuclei_label_bin = (aactnppa_nuclei_label > 0) & aactnppa_label_bin;
aactnppa_label = bwlabel(aactnppa_label_bin);
aactnppa_nuclei_label = aactnppa_label;
aactnppa_nuclei_label(~aactnppa_nuclei_label_bin) = 0;

hcn4_label_bin = hcn4_label > 0;
hcn4_label_bin = bwareaopen(hcn4_label_bin,10000);
hcn4_nuclei_label_bin = (hcn4_nuclei_label > 0) & hcn4_label_bin;
hcn4_label = bwlabel(hcn4_label_bin);
hcn4_nuclei_label = hcn4_label;
hcn4_nuclei_label(~hcn4_nuclei_label_bin) = 0;

hcn4nppa_label_bin = hcn4nppa_label > 0;
hcn4nppa_label_bin = bwareaopen(hcn4nppa_label_bin,10000);
hcn4nppa_nuclei_label_bin = (hcn4nppa_nuclei_label > 0) & hcn4nppa_label_bin;
hcn4nppa_label = bwlabel(hcn4nppa_label_bin);
hcn4nppa_nuclei_label = hcn4nppa_label;
hcn4nppa_nuclei_label(~hcn4nppa_nuclei_label_bin) = 0;

aacthcn4_label_bin = aacthcn4_label > 0;
aacthcn4_label_bin = bwareaopen(aacthcn4_label_bin,10000);
aacthcn4_nuclei_label_bin = (aacthcn4_nuclei_label > 0) & aacthcn4_label_bin;
aacthcn4_label = bwlabel(aacthcn4_label_bin);
aacthcn4_nuclei_label = aacthcn4_label;
aacthcn4_nuclei_label(~aacthcn4_nuclei_label_bin) = 0;

aacthcn4nppa_label_bin = aacthcn4nppa_label > 0;
aacthcn4nppa_label_bin = bwareaopen(aacthcn4nppa_label_bin,10000);
aacthcn4nppa_nuclei_label_bin = (aacthcn4nppa_nuclei_label > 0) & aacthcn4nppa_label_bin;
aacthcn4nppa_label = bwlabel(aacthcn4nppa_label_bin);
aacthcn4nppa_nuclei_label = aacthcn4nppa_label;
aacthcn4nppa_nuclei_label(~aacthcn4nppa_nuclei_label_bin) = 0;

cell_label = aact_label;
cell_label = cell_label + max(cell_label(:)).*aactnppa_label_bin + aactnppa_label;
cell_label = cell_label + max(cell_label(:)).*hcn4_label_bin + hcn4_label;
cell_label = cell_label + max(cell_label(:)).*hcn4nppa_label_bin + hcn4nppa_label;
cell_label = cell_label + max(cell_label(:)).*aacthcn4_label_bin + aacthcn4_label;
cell_label = cell_label + max(cell_label(:)).*aacthcn4nppa_label_bin + aacthcn4nppa_label;

segmentimages = struct();

segmentimages.nuclei_label = nuclei_label;
segmentimages.cell_label = cell_label;

segmentimages.aact_nuclei_label = aact_nuclei_label;
segmentimages.aact_label = aact_label;

segmentimages.aactnppa_nuclei_label = aactnppa_nuclei_label;
segmentimages.aactnppa_label = aactnppa_label;

segmentimages.hcn4_nuclei_label = hcn4_nuclei_label;
segmentimages.hcn4_label = hcn4_label;

segmentimages.hcn4nppa_nuclei_label = hcn4nppa_nuclei_label;
segmentimages.hcn4nppa_label = hcn4nppa_label;

segmentimages.aacthcn4_nuclei_label = aacthcn4_nuclei_label;
segmentimages.aacthcn4_label = aacthcn4_label;

segmentimages.aacthcn4nppa_nuclei_label = aacthcn4nppa_nuclei_label;
segmentimages.aacthcn4nppa_label = aacthcn4nppa_label;

end