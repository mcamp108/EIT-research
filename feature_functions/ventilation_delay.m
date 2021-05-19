function RVD = ventilation_delay(imgs, triads)
% -------------------------------------------------------------------------
% DESCRIPTION:
%   RVD = ventilation_delay(imgs, triads)
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
roi = horse_roi();
imgs = -imgs;    % impedance is high
isValidThresh = 0.15;
VDthresh = 0.4;
RVD = do_rvd_calc(imgs, triads, roi, isValidThresh, VDthresh);

end % end function

function RVD = do_rvd_calc(imgs, triads, roi, isValidThresh, VDthresh)
    nBreaths    = size(triads,1);
    nPixPerDim  = size(imgs,1);
    
    rsImgs      = reshape(imgs, [nPixPerDim^2, size(imgs,3)]);
    ltVent      = zeros(nBreaths,1);
    ltDors      = zeros(nBreaths,1);
    rtVent      = zeros(nBreaths,1);
    rtDors      = zeros(nBreaths,1);
    
    ventDelayImg = zeros(nPixPerDim, nPixPerDim);
    ventDelayImgCounts = zeros(nPixPerDim, nPixPerDim);
    
    for i = 1:nBreaths
        triad = triads(i, :);
        breathDeltaZ = calc_breath_delta_z(imgs, triad); % impedance change from end ex to end in
        breathDeltaZ(roi.BothLungs == 0) = 0;
        maxThisBreath = max(breathDeltaZ, [], 'all');
        pixelThresh = maxThisBreath * isValidThresh;
        breathDeltaZ(breathDeltaZ < pixelThresh) = 0;
        
        thisROI = breathDeltaZ > 0;
        thisROIIdx = find(thisROI);
        ltVentMask = find(thisROI + roi.LeftVentral == 2);
        ltDorsMask = find(thisROI + roi.LeftDorsal == 2);
        rtVentMask = find(thisROI + roi.RightVentral == 2);
        rtDorsMask = find(thisROI + roi.RightDorsal == 2);
        
%         thisBreathImgs = imgs(:, :, triad(1) : triad(2));
        thisRsImgs = rsImgs(thisROIIdx, triad(1) : triad(2));
        
        globalZ = nan(size(thisRsImgs, 2), 1);
        for j = 1 : size(thisRsImgs, 2)
            globalZ(j) = sum(thisRsImgs(:, j), 'all');
        end
        t1 = find(globalZ == min(globalZ));
        tN = find(globalZ == max(globalZ));
        thisBreathLen = tN - t1 + 1;
%         thisBreathImgs = thisBreathImgs(:, :, t1 : tN);
        thisRsImgs = thisRsImgs(:, t1 : tN);
        
        ventDelayVals = nan(length(thisROIIdx), 1);
        for j = 1 : size(thisRsImgs, 1)
           thisRsImgs(j, :) = normalize(thisRsImgs(j, :), 'range');
           ventDelayVals(j) = find(thisRsImgs(j, :) >= VDthresh, 1) / thisBreathLen;
        end
        
        thisVentDelayImg = zeros(nPixPerDim, nPixPerDim);
        thisVentDelayImg(thisROIIdx) = ventDelayVals;
        
        ventDelayImg(thisROIIdx) = ventDelayImg(thisROIIdx) + ventDelayVals;
        ventDelayImgCounts(thisROIIdx) = ventDelayImgCounts(thisROIIdx) + 1;
        
        % index results img for each roi and calc mean
        ltVent(i) = mean(thisVentDelayImg(ltVentMask), 'omitnan');
        ltDors(i) = mean(thisVentDelayImg(ltDorsMask), 'omitnan');
        rtVent(i) = mean(thisVentDelayImg(rtVentMask), 'omitnan');
        rtDors(i) = mean(thisVentDelayImg(rtDorsMask), 'omitnan');
        
        % RVD methods figure
%         rvd_methods_fig(thisVentDelayImg, thisRsImgs, thisROIIdx, ltVentMask, ltDorsMask, rtVentMask, rtDorsMask, thisBreathLen, roi);
        
    end % end for
    
    RVD = struct;
    RVD.LeftVentral = ltVent;
    RVD.LeftDorsal = ltDors;
    RVD.RightVentral = rtVent;
    RVD.RightDorsal = rtDors;
    
    bkg = isnan(imgs(:, :, 1));
    vdImg = ventDelayImg ./ ventDelayImgCounts;
    vdImg(ventDelayImgCounts == 0) = 0;
    vdImg(bkg) = nan;
    rvd_image(vdImg);
    
end % end function


function hrImg = make_high_res(img)
    inside = img == 0;
    bkg = isnan(img);
    lungs = img > 0;

    nCols = size(img, 1);
    nRows = size(img, 2);
    [X, Y] = meshgrid(1 : nRows, 1 : nCols);
    [X2, Y2] = meshgrid(1 : 0.25 : nRows, 1 : 0.25 : nCols);
    
    hrInside = interp2(X, Y, inside.*1, X2, Y2, 'linear');
    hrInside(hrInside < 0.5) = 0;
    hrInside(hrInside >= 0.5) = 1;
    
    hrBkg = interp2(X, Y, bkg.*1, X2, Y2, 'linear');
    hrBkg(hrBkg < 0.5) = 0;
    hrBkg(hrBkg >= 0.5) = 1;
    
    hrLung = interp2(X, Y, lungs.*1, X2, Y2, 'linear');
    hrLung(hrLung < 0.5) = 0;
    hrLung(hrLung >= 0.5) = 1;
    
    insideLungOvlp = hrLung + hrInside == 2;
    hrLung(insideLungOvlp) = 0;
    
    bkgInsideOvlp = hrInside + hrBkg == 2;
    hrBkg(bkgInsideOvlp) = 0;
    
    hrImg = interp2(X, Y, img, X2, Y2, 'linear');
    hrImg(hrBkg == 1) = 0;
    hrImg(hrInside == 1) = 1;
end

function rvd_image(ventDelayImg)
    fig=figure();
    fig.Units='normalized';
    fig.OuterPosition=[0.5 0.25 .5 .75];
    clf(); 
    fig.PaperOrientation='landscape';
    
    inside = ventDelayImg == 0;
    bkg = isnan(ventDelayImg);
    lungs = ventDelayImg > 0;
    vdImg = ventDelayImg;
    vdImg(isnan(vdImg)) = 0;
    vdImg(inside) = 1;
    
    betterImg = make_high_res(ventDelayImg);
    

    
    % show
%     imagesc(vdImg);
    imagesc(betterImg);
    axis equal; axis off; axis image;axis tight;
    cb = colorbar();
%     caxis([0 1]);
    
    % adjust labels
    nLbls = length(cb.TickLabels);
    lbls = linspace(0, 100, nLbls);
%     num2str(str2double(cb.TickLabels) * 100);
    cb.TickLabels = lbls(1 : end);
    cb.Label.String = 'Regional ventilation delay (%)';
    cb.FontSize = 12;
    cb.Label.FontSize = 18;
end


function rvd_methods_fig(thisVentDelayImg, thisRsImgs, thisROIIdx, ltVentMask, ltDorsMask, rtVentMask, rtDorsMask, thisBreathLen, roi)
    lv = mean(thisRsImgs(ismember(thisROIIdx, ltVentMask), :), 1);
    ld = mean(thisRsImgs(ismember(thisROIIdx, ltDorsMask), :), 1);
    rv = mean(thisRsImgs(ismember(thisROIIdx, rtVentMask), :), 1);
    rd = mean(thisRsImgs(ismember(thisROIIdx, rtDorsMask), :), 1);
    tmpImg = thisVentDelayImg;
    xax = linspace(0, 100, thisBreathLen);
    clf();
    subplot(1,2,1);
    tmpImg(roi.Inside == 0) = nan;
    tmpImg(tmpImg==0) = 0.02;
    tmpImg(roi.Inside == 0) = 0;
    imagesc(tmpImg);
    axis equal; axis image; axis tight; axis off;
    cmap = colormap('jet');
    cmap(1,:) = 1;
    cmap(2,:) = .9;
    colormap(cmap);
    cb = colorbar;
    caxis([0 1]);
    lbls = num2str(str2double(cb.TickLabels) * 100);
    cb.TickLabels = lbls;
    cb.Label.String = 'Regional ventilation delay (%)';
    cb.FontSize = 11;
    hold on;
    xline(17.5, 'linewidth', 4);
    yline(13.5, 'linewidth', 4);
    text(2,2,'Right ventral', 'fontsize', 12);
    text(31,2,'Left ventral', 'horizontalalignment','right', 'fontsize', 12);
    text(2,31,'Right dorsal', 'fontsize', 12);
    text(31,31,'Left dorsal', 'horizontalalignment','right', 'fontsize', 12);
    hold off;

    ax2 = subplot(1,2,2);
    title('Ventilation delay image');
    LW = 2;
    plot(xax, lv, 'linewidth', LW); hold on;
    plot(xax, ld, 'linewidth', LW);
    plot(xax, rv, 'linewidth', LW);
    plot(xax, rd, 'linewidth', LW);
    xlv = find(abs(lv-0.4) == min(abs(lv-0.4)));
    xld = find(abs(ld-0.4) == min(abs(ld-0.4)));
    xrv = find(abs(rv-0.4) == min(abs(rv-0.4)));
    xrd = find(abs(rd-0.4) == min(abs(rd-0.4)));

    plot([xax(xlv), xax(xlv)], [0, lv(xlv)], '-.', 'Color','#0072BD', 'linewidth', LW);
    plot([xax(xld), xax(xld)], [0, ld(xld)], '-.', 'Color','#D95319', 'linewidth', LW);
    plot([xax(xrv), xax(xrv)], [0, rv(xrv)], '-.', 'Color','#EDB120', 'linewidth', LW);
    plot([xax(xrd), xax(xrd)], [0, rd(xrd)], '-.', 'Color','#7E2F8E', 'linewidth', LW);         
    yline(0.4, '--k', 'linewidth', LW); hold off
    ax2.Position(4) = .615;
    ax2.Position(2) = .21;
    ax2.OuterPosition = [0.5343, 0.1270, 0.4097, 0.7546];
    ax2.XLabel.String = 'Regional ventilation delay (%)';
    ax2.YLabel.String = 'Normalized impedance';
    lg = legend();
    lg.String = {'Left ventral','Left dorsal','Right ventral','Right dorsal'};
    lg.Location = 'northwest';
end

