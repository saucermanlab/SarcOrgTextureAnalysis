function pf = radialProfile(I)
%% FUNCTION DESCRIPTION
% This function calculates the radial profile of the fourier transform.
%
% INPUTS
% I : Fourier transform image
%
% OUTPUTS
% pf : radial profile array

sz = size(I);
middleIndex = round(sz./2);
pf = zeros(1,floor(min(sz./2)));


for iAngle = 0:0.5:179.5
    II = imrotate(I,iAngle,'crop');
    pf = pf + II(middleIndex(1),middleIndex(2):(middleIndex(2)+numel(pf)-1));
end