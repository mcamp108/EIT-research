function deltas = mean_end_ex_deltas(imgs, triads)
% -------------------------------------------------------------------------
% DESCRIPTION:
%   deltas = mean_end_ex_deltas(imgs, triads)
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

for i = 1: nBreaths
    triad = triads(i, :);
    endEx1 = imgs(:, :, triad(1));
    endEx2 = imgs(:, :, triad(3));
    mean_end_ex = (endEx1 + endEx2) ./ 2;
    
    ltVent(i) = sum(mean_end_ex(ltVentMask), 'all');
    ltDors(i) = sum(mean_end_ex(ltDorsMask), 'all');
    rtVent(i) = sum(mean_end_ex(rtVentMask), 'all');
    rtDors(i) = sum(mean_end_ex(rtDorsMask), 'all');
    bothLungs(i) = sum(mean_end_ex(bothLngMask), 'all');
end

% calc deltas
for i = 1: nBreaths
    ltVent(i) = ltVent(i) - ltVent(1);
    ltDors(i) = ltDors(i) - ltDors(1);
    rtVent(i) = rtVent(i) - rtVent(1);
    rtDors(i) = rtDors(i) - rtDors(1);
    bothLungs(i) = bothLungs(i) - bothLungs(1);
end

deltas = struct;
deltas.bothLungs = bothLungs;
deltas.ltVent = ltVent;
deltas.ltDors = ltDors;
deltas.rtVent = rtVent;
deltas.rtDors = rtDors;

end % end function