function endEx2Means = end_ex2_mean_Z(imgs, triads)
% -------------------------------------------------------------------------
% DESCRIPTION:
%   endEx2Means = end_ex2_mean_Z(imgs, triads)
% -------------------------------------------------------------------------
% PARAMETERS:
%   imgs (array):
%       images from calc_slices
%   triads (n x 3 array):
%       breath boundaries
% -------------------------------------------------------------------------   
% RETURNS:
%   lower values mean less air
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
imgs        = -imgs; % impedance is high - lower values mean less air
roi         = horse_roi();
nBreaths    = size(triads,1);
ltVent      = zeros(nBreaths,1);
ltDors      = zeros(nBreaths,1);
rtVent      = zeros(nBreaths,1);
rtDors      = zeros(nBreaths,1);
bothLungs   = zeros(nBreaths,1);

ltVentMask  = find(roi.LeftVentral);
ltDorsMask  = find(roi.LeftDorsal);
rtVentMask  = find(roi.RightVentral);
rtDorsMask  = find(roi.RightDorsal);
bothLngMask = find(roi.BothLungs);

for i = 1:nBreaths
    triad       = triads(i, :);
    endEx2      = imgs(:,:,triad(3));
    ltVent(i)   = sum(endEx2(ltVentMask), 'all') / length(ltVentMask);
    ltDors(i)   = sum(endEx2(ltDorsMask), 'all') / length(ltDorsMask);
    rtVent(i)   = sum(endEx2(rtVentMask), 'all') / length(rtVentMask);
    rtDors(i)   = sum(endEx2(rtDorsMask), 'all') / length(rtDorsMask);
    bothLungs(i)= sum(endEx2(bothLngMask), 'all') / length(bothLngMask);
end

endEx2Means          = struct;
endEx2Means.ltVent   = ltVent;
endEx2Means.ltDors   = ltDors;
endEx2Means.rtVent   = rtVent;
endEx2Means.rtDors   = rtDors;
endEx2Means.bothLungs= bothLungs;

end % end function