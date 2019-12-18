%% Idealized cells
% This script generates the idealized images used in Figure 3.

clear; close all; clc;
rng(0);
%% Horizontal Stripes
n = 2^7;
I = zeros(n);
l = 4;
d = 4;

for iRow = 1:n
    for iCol = 1:n
        if (mod(iRow-1,l+d) >= d)
            I(iRow,iCol) = 1;
        else
            I(iRow,iCol) = 0;
        end
    end
end

% figure; imshow(I,'InitialMagnification',200);
% imwrite(I,'.\Images\Ideal\Ideal 01.tif');
metrics = morph_texture_function_nointerp(I,ones(size(I)),0:1:179,0:1:40,1,1,'IdealizedCell');
plotCorrelation_idealized(metrics,0,1);

%% Diagonal Stripes
n = 2^8;
I = zeros(n);
l = 4;
d = 4;

for iRow = 1:n
    for iCol = 1:n
        if (mod(iRow-1,l+d) >= d)
            I(iRow,iCol) = 1;
        else
            I(iRow,iCol) = 0;
        end
    end
end
I = imrotate(I,-45,'bilinear');
s = size(I);
I = I((s(1)/2-n/2/2):(s(1)/2+n/2/2-1),(s(2)/2-n/2/2):(s(2)/2+n/2/2-1));
% figure; imshow(I,'InitialMagnification',200);
% imwrite(I,'.\Images\Ideal\Ideal 02.tif');
metrics = morph_texture_function_nointerp(I,ones(size(I)),0:1:179,0:1:40,1,1,'IdealizedCell');
plotCorrelation_idealized(metrics,0,1);

%% Horizontal Stripes - Longer Wavelength
n = 2^7;
I = zeros(n);
l = 8;
d = 8;

for iRow = 1:n
    for iCol = 1:n
        if (mod(iRow-1,l+d) >= d)
            I(iRow,iCol) = 1;
        else
            I(iRow,iCol) = 0;
        end
    end
end

% figure; imshow(I,'InitialMagnification',200);
% imwrite(I,'.\Images\Ideal\Ideal 03.tif');
metrics = morph_texture_function_nointerp(I,ones(size(I)),0:1:179,0:1:40,1,1,'IdealizedCell');
plotCorrelation_idealized(metrics,0,1);

%% Horizontal Stripes - Noise Added
n = 2^7;
I = zeros(n);
l = 4;
d = 4;

for iRow = 1:n
    for iCol = 1:n
        if (mod(iRow-1,l+d) >= d)
            I(iRow,iCol) = 0.7+2*(rand-0.5)*0.3;
        else
            I(iRow,iCol) = 0.3+2*(rand-0.5)*0.3;
        end
    end
end

% figure; imshow(I,'InitialMagnification',200);
% imwrite(I,'.\Images\Ideal\Ideal 04.tif');
metrics = morph_texture_function_nointerp(I,ones(size(I)),0:1:179,0:1:40,1,1,'IdealizedCell');
plotCorrelation_idealized(metrics,0,1);

%% Horizontal Stripes - 1/3 Image
n = 2^7;
I = zeros(n);
l = 4;
d = 4;

for iRow = 1:n
    for iCol = 1:n
        if (mod(iRow-1,l+d) >= d) && (iRow <= n/3)
            I(iRow,iCol) = 1;
        else
            I(iRow,iCol) = 0;
        end
    end
end

% figure; imshow(I,'InitialMagnification',200);
% imwrite(I,'.\Images\Ideal\Ideal 05.tif');
metrics = morph_texture_function_nointerp(I,ones(size(I)),0:1:179,0:1:40,1,1,'IdealizedCell');
plotCorrelation_idealized(metrics,0,1);

%% Horizontal Stripes - Vertical Bands
n = 2^7;
I = zeros(n);
l = 4;
d = 4;
w = 8;
r = repelem(randi([0 l+d-1],1,n/w),w);
for iRow = 1:n
    for iCol = 1:n
        if (mod(iRow-1+r(iCol),l+d) >= d)
            I(iRow,iCol) = 1;
        else
            I(iRow,iCol) = 0;
        end
    end
end

% figure; imshow(I,'InitialMagnification',200);
% imwrite(I,'.\Images\Ideal\Ideal 06.tif');
metrics = morph_texture_function_nointerp(I,ones(size(I)),0:1:179,0:1:40,1,1,'IdealizedCell');
plotCorrelation_idealized(metrics,0,1);
