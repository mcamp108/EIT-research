function ventDist = ventilation_distribution(imgs, triads)
% -------------------------------------------------------------------------
% DESCRIPTION:
%   ventDist = ventilation_distribution(imgs, triads)
%   Assumes animal in imgs is in dorsal recumbancy
% -------------------------------------------------------------------------
% PARAMETERS:
%   imgs:
%       EIT images for full time series
%   triads:
%       n x 3 array of indices in images for triad (end_ex1, end_in,
%       end_ex2)
% -------------------------------------------------------------------------   
% RETURNS:
%   ventDist (struct):
%       relative ventilation distribution of region determined by relative
%       distrubution of impedance change from average end_ex to end_in left
%       / right, dorsal / ventral
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

imgs        = -imgs;
roi         = horse_roi();
ventDist    = struct;
nBreaths    = size(triads,1);
ltVent      = zeros(nBreaths,1);
ltDors      = zeros(nBreaths,1);
rtVent      = zeros(nBreaths,1);
rtDors      = zeros(nBreaths,1);
bothLungs   = roi.BothLungs;

ventDistImg = zeros(32, 32);
ventDistImgCounts = zeros(32, 32);

for i = 1 : nBreaths
    triad = triads(i, :);
    breathZ = calc_breath_delta_z(imgs, triad);
    lungZ = breathZ(bothLungs==1);
    lungZ(lungZ < 0) = 0;

    totalZ = sum(lungZ,'all');
    img = zeros(32,32);
    img(bothLungs == 1) = lungZ / totalZ; % percentage of total lung impedance represented by each pixel
    mask = find(img > 0);
    ventDistImg(mask) = ventDistImg(mask) + img(mask);
    ventDistImgCounts(mask) = ventDistImgCounts(mask) + 1;
    
    [LV, LD, RV, RD] = calc_vent_distribution(img, roi);
    ltVent(i) = LV;
    ltDors(i) = LD;
    rtVent(i) = RV;
    rtDors(i) = RD;
    rowMeans = calc_norm_lung_row_means(img, roi);
    if ~exist('rowMeanRes', 'var')
        rowMeanRes = zeros(nBreaths, length(rowMeans));
    end
    rowMeanRes(i,:) = rowMeans;
end % end for

ventDist.ltVent = ltVent;
ventDist.ltDors = ltDors;
ventDist.rtVent = rtVent;
ventDist.rtDors = rtDors;

for i = 1:size(rowMeanRes, 2)
    fn = sprintf('rowMeanRes%.0f', i);
    ventDist.(fn) = rowMeanRes(:,i);
end

% make figure
bkg = isnan(imgs(:, :, 1));
vdImg = ventDistImg ./ ventDistImgCounts;
vdImg(ventDistImgCounts == 0) = 0;
vdImg(bkg) = nan;
vent_dist_image(vdImg);

end % end function

function [LV,LD,RV,RD] = calc_vent_distribution(img, roi)
    ltVentMask = find(roi.LeftVentral);
    ltDorsMask = find(roi.LeftDorsal);
    rtVentMask = find(roi.RightVentral);
    rtDorsMask = find(roi.RightDorsal);
    % TODO: these should be the sums, rather than the division and the
    % summing then more division!
    zLV     = sum(img(ltVentMask), 'all') ./ sum(img(ltVentMask) > 0, 'all');
    zLD     = sum(img(ltDorsMask), 'all') ./ sum(img(ltDorsMask) > 0, 'all');
    zRV     = sum(img(rtVentMask), 'all') ./ sum(img(rtVentMask) > 0, 'all');
    zRD     = sum(img(rtDorsMask), 'all') ./ sum(img(rtDorsMask) > 0, 'all');

    zTotal  = zLV + zLD + zRV + zRD;
    LV      = (zLV / zTotal);
    LD      = (zLD / zTotal);
    RV      = (zRV / zTotal);
    RD      = (zRD / zTotal);
    
    assert( round(LV + LD + RV + RD, 10) == 1 );
    assert( sum([LV, LD, RV, RD] > 0) == 4 );
end % end function


function hrImg = make_high_res(img, cbot, clim)
    if nargin == 1
        cbot = 0;
        clim = 1;
    end
    inside = img == 0;
    bkg = isnan(img);
    lungs = img > 0;
    
    nCols = size(img, 1);
    nRows = size(img, 2);
    [X, Y] = meshgrid(1 : nRows, 1 : nCols);
    [X2, Y2] = meshgrid(1 : 0.25 : nRows, 1 : 0.25 : nCols);
    
    % find lungs
    hrImg = interp2(X, Y, img, X2, Y2, 'linear');
    hrLung = hrImg >= min(img(lungs), [], 'omitnan');
    
    % find inside, not lung
    hrInside = interp2(X, Y, inside.*1, X2, Y2, 'linear');
    hrInside(hrInside < 0.5) = 0;
    hrInside(hrInside >= 0.5) = 1;
    
    % find background
    hrBkg = interp2(X, Y, bkg.*1, X2, Y2, 'linear');
    hrBkg(hrBkg < 0.5) = 0;
    hrBkg(hrBkg >= 0.5) = 1;
    
    % overlap between inside and background becomes inside
    bkgInsideOvlp = hrInside + hrBkg == 2;
    hrBkg(bkgInsideOvlp) = 0;
    
    % overlap between inside and lung becomes lung
    insideLungOvlp = hrLung + hrInside == 2;
    hrInside(insideLungOvlp) = 0;
    
    % all unassigned pixels are inside
    all = hrLung + hrInside + hrBkg;
    hrInside(all == 0) = 1;
    
    hrImg(hrBkg == 1) = clim;
    hrImg(hrInside == 1) = cbot;
end


function vent_dist_image(ventDistImg)
    fig = figure();
    fig.Units = 'normalized';
    fig.OuterPosition = [0.5 0.25 .5 .75];
    clf(); 
    fig.PaperOrientation = 'landscape';
    
    clim = 0.0175;
    lungs = ventDistImg > 0;
    ctop = max(ventDistImg(lungs));
    if ctop >= clim
       error();
    end    
    betterImg = make_high_res(ventDistImg, 0, clim);  

    % adjust cmap
    colormap('jet');
    cmap = colormap; %get current colormap
    cmap(1, :) = 0;
    cmap(end, :) = 1;
    colormap(cmap); % apply new colormap
    
    % show
    imagesc(betterImg);
    axis equal; axis off; axis image;axis tight;
    cb = colorbar();
    lbls = num2str(str2double(cb.TickLabels) * 100);
    cb.TickLabels = lbls;
    cb.Label.String = 'Regional ventilation distribution (%)';
    cb.FontSize = 12;
    cb.Label.FontSize = 18;
end