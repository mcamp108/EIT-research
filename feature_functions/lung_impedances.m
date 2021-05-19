function lungZ = lung_impedances(imgs, triads)
% -------------------------------------------------------------------------
% DESCRIPTION:
%   lungZ = lung_impedances(imgs, triads)
% -------------------------------------------------------------------------
% PARAMETERS:
%   imgs (array):
%       images from calc_slices
%   triads (n x 3 array):
%       breath boundaries
% -------------------------------------------------------------------------   
% RETURNS:
%   lungZ (struct):
%       has fields EELI1, EILI, EELI2 which are the sum of lung impedances
%       at each point of the breath triad.
% -------------------------------------------------------------------------   
% AUTHOR:
%   Mark Campbell
%   Carleton University
%   markacampbell@sce.carleton.ca
% -------------------------------------------------------------------------
% VERSION:
%   1.0.0
% -------------------------------------------------------------------------
% (C) 2019-2021 Mark Campbell. License: GPL version 2 or version 3
% -------------------------------------------------------------------------
imgs    = -imgs; % impedance is high - lower values mean less air
roi     = horse_roi();
nBreaths= size(triads,1);
EELI1   = zeros(nBreaths,1);
EELI1_ltVent = zeros(nBreaths,1);
EELI1_ltDors = zeros(nBreaths,1);
EELI1_rtVent = zeros(nBreaths,1);
EELI1_rtDors = zeros(nBreaths,1);

EILI    = zeros(nBreaths,1);
EILI_ltVent = zeros(nBreaths,1);
EILI_ltDors = zeros(nBreaths,1);
EILI_rtVent = zeros(nBreaths,1);
EILI_rtDors = zeros(nBreaths,1);

EELI2   = zeros(nBreaths,1);
EELI2_ltVent = zeros(nBreaths,1);
EELI2_ltDors = zeros(nBreaths,1);
EELI2_rtVent = zeros(nBreaths,1);
EELI2_rtDors = zeros(nBreaths,1);

masks.ltVentMask  = find(roi.LeftVentral);
masks.ltDorsMask  = find(roi.LeftDorsal);
masks.rtVentMask  = find(roi.RightVentral);
masks.rtDorsMask  = find(roi.RightDorsal);
masks.bothLngMask = roi.BothLungs;

for i = 1:nBreaths
    endEx1  = imgs(:, :, triads(i, 1));
    endIn   = imgs(:, :, triads(i, 2));
    endEx2  = imgs(:, :, triads(i, 3));
    
    [EELI1(i), EELI1_ltVent(i), EELI1_ltDors(i), EELI1_rtVent(i), EELI1_rtDors(i)] = Lung_Z(endEx1, masks);
    
    [EILI(i), EILI_ltVent(i), EILI_ltDors(i), EILI_rtVent(i), EILI_rtDors(i)] = Lung_Z(endIn, masks);
    
    [EELI2(i), EELI2_ltVent(i), EELI2_ltDors(i), EELI2_rtVent(i), EELI2_rtDors(i)] = Lung_Z(endEx2, masks);

end

lungZ = struct;
lungZ.EELI1 = EELI1;
lungZ.EELI1_ltVent   = EELI1_ltVent;
lungZ.EELI1_ltDors   = EELI1_ltDors;
lungZ.EELI1_rtVent   = EELI1_rtVent;
lungZ.EELI1_rtDors   = EELI1_rtDors;

lungZ.EILI  = EILI;
lungZ.EILI_ltVent   = EILI_ltVent;
lungZ.EILI_ltDors   = EILI_ltDors;
lungZ.EILI_rtVent   = EILI_rtVent;
lungZ.EILI_rtDors   = EILI_rtDors;

lungZ.EELI2 = EELI2;
lungZ.EELI2_ltVent   = EELI2_ltVent;
lungZ.EELI2_ltDors   = EELI2_ltDors;
lungZ.EELI2_rtVent   = EELI2_rtVent;
lungZ.EELI2_rtDors   = EELI2_rtDors;

end % end function

function [both, ltVent, ltDors, rtVent, rtDors] = Lung_Z(frame, masks)
    both = sum(frame(masks.bothLngMask), 'all');
    ltVent = sum(frame(masks.ltVentMask), 'all');
    ltDors = sum(frame(masks.ltDorsMask), 'all');
    rtVent = sum(frame(masks.rtVentMask), 'all');
    rtDors = sum(frame(masks.rtDorsMask), 'all');
end