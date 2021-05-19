function [onOpening, onClosing] = pixels_not_exceeding_pct_lung_mean(imgs, triads, pct)
% -------------------------------------------------------------------------
% DESCRIPTION:
%   [onOpening, onClosing] = pixels_not_exceeding_pct_lung_mean(imgs, triads, pct)
% -------------------------------------------------------------------------
% PARAMETERS:
% 
% -------------------------------------------------------------------------   
% RETURNS:
%   Percent of pixels in roi that did not exceed pct (%) of the max mean
%   pixel change from end ex1 to end in (onOpening) or from end in to end
%   ex2 (onClosing)
% -------------------------------------------------------------------------   
% AUTHOR:
%   Mark Campbell
%   Carleton University
%   mark@analyti.ca
% -------------------------------------------------------------------------
% VERSION:
%   1.0.0
% -------------------------------------------------------------------------
% (C) 2019-2021 Mark Campbell. License: GPL version 2 or version 3
% -------------------------------------------------------------------------
roi             = horse_roi();
openImgs        = -imgs;    % impedance is high
closeImgs       = imgs;     % impedance is low
closeTriads     = triads;
closeTriads(:,1)= triads(:,2);
closeTriads(:,2)= triads(:,3);

onOpening = do_pixels_not_exceeding_pct_lung_mean(openImgs, triads, roi, pct);
onClosing = do_pixels_not_exceeding_pct_lung_mean(closeImgs, closeTriads, roi, pct);

end % end function


function pctPix = do_pixels_not_exceeding_pct_lung_mean(imgs, triads, roi, pct)
    nBreaths    = size(triads,1);
    bothLungs   = (roi.RightLung + roi.LeftLung);
    mask        = find(bothLungs);
    nPixPerDim  = size(imgs,1);
    
    rsImgs      = reshape(imgs, [nPixPerDim^2, size(imgs,3)]);
    ltVent      = zeros(nBreaths,1);
    ltDors      = zeros(nBreaths,1);
    rtVent      = zeros(nBreaths,1);
    rtDors      = zeros(nBreaths,1);

    ltVentMask  = find(roi.LeftVentral);
    ltDorsMask  = find(roi.LeftDorsal);
    rtVentMask  = find(roi.RightVentral);
    rtDorsMask  = find(roi.RightDorsal);

    for i = 1:nBreaths
        triad       = triads(i, :);
        ltVent(i)   = find_valid_lung_pixels_in_image(imgs, triad, ltVentMask, pct);
        ltDors(i)   = find_valid_lung_pixels_in_image(imgs, triad, ltDorsMask, pct);
        rtVent(i)   = find_valid_lung_pixels_in_image(imgs, triad, rtVentMask, pct);
        rtDors(i)   = find_valid_lung_pixels_in_image(imgs, triad, rtDorsMask, pct);
    end
    
    pctPix                = struct;
    pctPix.LeftVentral    = ltVent;
    pctPix.LeftDorsal     = ltDors;
    pctPix.RightVentral   = rtVent;
    pctPix.RightDorsal    = rtDors;
end

function [included] = find_valid_lung_pixels_in_image(imgs, triad, mask, pct)
    breathDeltaZ        = calc_breath_delta_z(imgs, triad); % impedance change from end ex to end in
    LungDeltaZ          = breathDeltaZ(mask); % isolate lung delta Z
    pctMeanLungDeltaZ   = mean(LungDeltaZ, 'all') * (pct/100);
    included            = sum(LungDeltaZ >= pctMeanLungDeltaZ) / length(mask);
end