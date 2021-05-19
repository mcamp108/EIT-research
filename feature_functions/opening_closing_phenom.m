function [openingPhenom, closingPhenom] = opening_closing_phenom(imgs, triads, fs, pct)
% -------------------------------------------------------------------------
% DESCRIPTION:
%   [openingPhenom, closingPhenom] = opening_closing_phenom(imgs, triads) 
% 
%   Assumes animal in imgs is in dorsal recumbancy.
% 
%   For each pixel, the relative time taken to cross the threshold
%   impedance value (thresh) starting from end exhilation is calculated
%   (estimator of compliance). The mean time of pixels is then calculated
%   per region of interest.
% 
%   After end expiration, the moment at which an individual waveform
%   reaches 10% of its maximum tidal change (or overall change determined
%   through a low-flow inflation maneuver) is identified and the
%   corresponding value of instantaneous airway pressure is used to compose
%   a functional image. Pixels exhibiting tidal changes less than half the
%   average were excluded.
% -------------------------------------------------------------------------
% PARAMETERS:
%   imgs:
%       EIT images for full time series
%   triads:
%       n x 3 array of indices in images for triad (end_ex1, end_in,
%       end_ex2)
% -------------------------------------------------------------------------   
% RETURNS:
%   out (struct):
%       average time per region to reach threshold impedance from average
%       end_ex to end_in for left / right, dorsal / ventral
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
roi             = horse_roi();
thresh          = 0.1;
openImgs        = -imgs;    % impedance is high
closeImgs       = imgs;     % impedance is low
closeTriads     = triads;
closeTriads(:,1)= triads(:,2);
closeTriads(:,2)= triads(:,3);

openingPhenom   = do_phenom_calc(openImgs, triads, roi, thresh, fs, pct);
closingPhenom   = do_phenom_calc(closeImgs, closeTriads, roi, 1-thresh, fs, pct);

end % end function

function ocPhenom = do_phenom_calc(imgs, triads, roi, thresh, fs, pct)
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
        triad               = triads(i, :);
        pctMeanLungDeltaZ   = calc_Z_thresh(imgs, triad, mask, pct);
        threshZ     = calc_Z_thresh(imgs, triad, mask, thresh);
        pixWaveform = rsImgs( mask, triads(i,1):triads(i,2) );
        pixWaveform = pixWaveform - pixWaveform(:,1); % subtract baseline
        resultsImg  = nan(nPixPerDim,nPixPerDim);
        
        for j = 1:length(mask)
            pixel       = pixWaveform(j, :);
            isValid     = find(pixel >= pctMeanLungDeltaZ, 1, 'first');
            if isempty(isValid)
                continue
            end
            crossThresh = find(pixel >= threshZ, 1, 'first');
            if isempty(crossThresh)
                continue
            end
            resultsImg(mask(j)) = crossThresh;
        end    

        % index results img for each roi and calc mean
        ltVent(i) = mean(resultsImg(ltVentMask), 'omitnan') / fs;
        ltDors(i) = mean(resultsImg(ltDorsMask), 'omitnan') / fs;
        rtVent(i) = mean(resultsImg(rtVentMask), 'omitnan') / fs;
        rtDors(i) = mean(resultsImg(rtDorsMask), 'omitnan') / fs;
    end % end for
    
    ocPhenom                = struct;
    ocPhenom.LeftVentral    = ltVent;
    ocPhenom.LeftDorsal     = ltDors;
    ocPhenom.RightVentral   = rtVent;
    ocPhenom.RightDorsal    = rtDors;
end

function pctMeanLungDeltaZ = calc_Z_thresh(imgs, triad, mask, pct)
    breathDeltaZ        = calc_breath_delta_z(imgs, triad); % impedance change from end ex to end in
    LungDeltaZ          = breathDeltaZ(mask); % isolate lung delta Z
    pctMeanLungDeltaZ   = mean(LungDeltaZ, 'all')  * (pct/100);
end

