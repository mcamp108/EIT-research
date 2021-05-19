function ETC = expiratory_time_constant(imgs, triads, fs, mode)
% -------------------------------------------------------------------------
% DESCRIPTION:
%   Reference: Regional expiratory time constants in severe respiratory
%   failure estimated by electrical impedance tomography: a feasibility
%   study
% 
%   Time constants were calculated by exponential fitting for every single
%   pixel. If the single pixel had no exponential fit, e.g. pixel from the
%   heart or atelectasis without ventilation, this pixel was excluded from
%   further analysis. Furthermore, the global impedance (?Z(t)global)
%   signal was calculated as the sum of all signals at any point of time.
%   The start and end of expiration were determined, respectively from the
%   global EIT signal. Then within this time sequence, for each single
%   pixel the 75% amplitude was calculated as the start of the curve
%   fitting (Fig. 2b) and a local minimum was used to define the end of
%   expiration. However, the global signal is only used to define the time
%   sequence of the expiration but start and end was defined pixelwise.
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
if nargin < 4
    mode = 'regional';
end    
imgs = -imgs;
nBreaths = size(triads, 1);
roi = horse_roi();
g = fittype('a*exp(-x/tau)+c');
nPixPerDim = size(imgs, 1);
regional = false;
pixelwise = false;
if strcmp(mode, 'regional')
    regional = true;
    ltVentMask = find(roi.LeftVentral);
    ltDorsMask = find(roi.LeftDorsal);
    rtVentMask = find(roi.RightVentral);
    rtDorsMask = find(roi.RightDorsal);
    glblMask = find(roi.BothLungs);
    
    ltVent = zeros(nBreaths, 1);
    ltDors = zeros(nBreaths, 1);
    rtVent = zeros(nBreaths, 1);
    rtDors = zeros(nBreaths, 1);
    glbl = zeros(nBreaths, 1);
elseif strcmp(mode, 'pixelwise')
    pixelwise = true;
    ETC = zeros(nPixPerDim, nPixPerDim, nBreaths);
else
    keyboard;
end

etcImg = zeros(nPixPerDim, nPixPerDim);
etcImgCounts = zeros(nPixPerDim, nPixPerDim);

for i = 1:nBreaths
    thisTriad = triads(i, :);
    breathZ = calc_breath_delta_z(imgs, thisTriad);
    breathZ = breathZ .* roi.BothLungs;
    includedPixelIdx = find(breathZ > 0);
    expirationImgs = imgs(:, :, thisTriad(2) : thisTriad(3));
    rsExpirationImgs = reshape(expirationImgs, [nPixPerDim^2, size(expirationImgs,3)]);
    tauImg = nan(size(breathZ));
    for j = 1:length(includedPixelIdx)
        pixelIdx = includedPixelIdx(j);
        try
            tauImg(pixelIdx) = calc_pixel_tau(g, rsExpirationImgs(pixelIdx, :), fs);
        catch
            continue
        end
    end
    if pixelwise
        ETC(:,:,i) = tauImg;
    elseif regional
        tempMask = ltVentMask(ismember(ltVentMask, includedPixelIdx));
        ltVent(i) = mean(tauImg(tempMask), 'omitnan');
        
        tempMask = ltDorsMask(ismember(ltDorsMask, includedPixelIdx));
        ltDors(i) = mean(tauImg(tempMask), 'omitnan');
        
        tempMask = rtVentMask(ismember(rtVentMask, includedPixelIdx));
        rtVent(i) = mean(tauImg(tempMask), 'omitnan');
        
        tempMask = rtDorsMask(ismember(rtDorsMask, includedPixelIdx));
        rtDors(i) = mean(tauImg(tempMask), 'omitnan');
        
        tempMask = glblMask(ismember(glblMask, includedPixelIdx));
        glbl(i) = mean(tauImg(tempMask), 'omitnan');
        
        etcImg(tempMask) = etcImg(tempMask) + tauImg(tempMask);
        etcImgCounts(tempMask) = etcImgCounts(tempMask) + 1;
%         try
%             glbl(i) = calc_pixel_tau(g, sum(rsExpirationImgs(includedPixelIdx, :), 1), fs);
%         catch
%             glbl(i) = nan;
%         end
    end
end

if regional
    ETC.glbl = glbl;
    ETC.ltVent = ltVent;
    ETC.ltDors = ltDors;
    ETC.rtVent = rtVent;
    ETC.rtDors = rtDors;
end

etc_image(etcImg ./ etcImgCounts);

end % end function


function tau = calc_pixel_tau(g, expiration, fs)
    expiration_ = (expiration - min(expiration)) ./ (max(expiration) - min(expiration));
    expiration = expiration_;
    start = find(expiration >= 0.75);
    if length(start) > 1
        start = start(end);
    end
    expiration = expiration(start : end);

    endExp = find(expiration == min(expiration));
    
    expiration = expiration(1 : endExp);
    x_ = (1 : endExp);
    % fit exponential curve
    f0 = fit(x_', expiration', g, 'StartPoint', [.75, 0, fs * 1.25]);
    tau = f0.tau / fs;
    assert(tau <= 5);
%     if tau > (endExp / fs)
%         clf();
%         subplot(2,1,1)
%         plot(expiration_);
%         hold on;
%         xline(start);
%         xline(start + endExp - 1);
%         hold off;
%         
%         subplot(2,1,2)
%         xx = linspace(1, length(expiration), 500);
%         plot(expiration); hold on;
%         plot(xx, f0(xx), 'r-');
%         hold off;
%         keyboard
%     end
end


function etc_image(etcImg)
    fig = figure();
    fig.Units='normalized';
    fig.OuterPosition = [0.5 0.25 .5 .75];
    clf(); 
    fig.PaperOrientation = 'landscape';
    etcImg(isnan(etcImg)) = 0;
    imagesc(etcImg); axis image; axis off; axis equal; axis tight;    
    colormap('jet');
    cmap = colormap; %get current colormap
    cmap(1, :) = 0;
    colormap(cmap); % apply new colormap
    cb = colorbar();
    caxis([0 3]);
    cb.Label.String = 'Expiratory Time Constant (s)';
    cb.FontSize = 12;
    cb.Label.FontSize = 18;
end