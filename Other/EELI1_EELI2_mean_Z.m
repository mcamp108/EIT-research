function meanEELI = EELI1_EELI2_mean_Z(imgs, triads)
% -------------------------------------------------------------------------
% DESCRIPTION:
%
% -------------------------------------------------------------------------
% PARAMETERS:
% 
% -------------------------------------------------------------------------   
% RETURNS:
% 
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
nBreaths    = size(triads, 1);
ltVent      = zeros(nBreaths, 1);
ltDors      = zeros(nBreaths, 1);
rtVent      = zeros(nBreaths, 1);
rtDors      = zeros(nBreaths, 1);
bothLungs   = zeros(nBreaths, 1);

ltVentMask  = find(roi.LeftVentral);
ltDorsMask  = find(roi.LeftDorsal);
rtVentMask  = find(roi.RightVentral);
rtDorsMask  = find(roi.RightDorsal);
bothLngMask = find(roi.BothLungs);

for i = 1: nBreaths
    triad = triads(i, :);
    eeli_1 = imgs(:, :, triad(1));
    eeli_2 = imgs(:, :, triad(3));
    mean_eeli = (eeli_1 + eeli_2) ./ 2;
    ltVent(i) = sum(mean_eeli(ltVentMask), 'all') / length(ltVentMask);
    ltDors(i) = sum(mean_eeli(ltDorsMask), 'all') / length(ltDorsMask);
    rtVent(i) = sum(mean_eeli(rtVentMask), 'all') / length(rtVentMask);
    rtDors(i) = sum(mean_eeli(rtDorsMask), 'all') / length(rtDorsMask);
    bothLungs(i) = sum(mean_eeli(bothLngMask), 'all') / length(bothLngMask);
end

meanEELI = struct;
meanEELI.ltVent = ltVent;
meanEELI.ltDors = ltDors;
meanEELI.rtVent = rtVent;
meanEELI.rtDors = rtDors;
meanEELI.bothLungs = bothLungs;

end % end function