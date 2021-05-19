function perfDist = perfusion_distribution(eitFile)
% -------------------------------------------------------------------------
% DESCRIPTION:
%   perfDist = perfusion_distribution(imgr, eitFile)
% -------------------------------------------------------------------------
% PARAMETERS:
%   imgr (EIDORS reconstructed image object)
%   eitFile (EIT file object with perfStart, perfEnd, and imgr fields)
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

NPOINTS     = 5;
roi         = horse_roi();
ltVentMask  = find(roi.LeftVentral);
ltDorsMask  = find(roi.LeftDorsal);
rtVentMask  = find(roi.RightVentral);
rtDorsMask  = find(roi.RightDorsal);
bothLungs   = roi.BothLungs;
lungMask    = find(bothLungs);

refXVal     = 47;
start       = eitFile.perfStart - refXVal;
stop        = eitFile.perfEnd;
nFrames     = stop-start+1;

perfDist    = struct;
ltVent      = zeros(NPOINTS,1);
ltDors      = zeros(NPOINTS,1);
rtVent      = zeros(NPOINTS,1);
rtDors      = zeros(NPOINTS,1);

imgr = eitFile.imgr;
perfImgr = imgr;
perfImgr.elem_data = imgr.elem_data(:,start:stop);
clim = max(perfImgr.elem_data,[],'all');
perfImgr.calc_colours.ref_level  = 0;
perfImgr.calc_colours.clim       = clim;
perfImgs = calc_slices(perfImgr);

timePoints = linspace(refXVal, nFrames, NPOINTS+2);
timePoints = round(timePoints(2:NPOINTS+1));

if strcmp(eitFile.name, '296439 IE11 Perfusion breathing')
    timePoints(3) = 1380 - eitFile.perfStart + 1;
elseif strcmp(eitFile.name, '296439 IE13 Perfusion breathing')
    timePoints(3) = 2888 - eitFile.perfStart + 1;
end

for i=1:NPOINTS
    perfZ = perfImgs(:,:,timePoints(i));
    perfZ = perfZ - perfImgs(:, :, 1);
    lungZ = perfZ(bothLungs == 1);
    lungZ(lungZ < 0) = 0;
    
    totalZ = sum(lungZ, 'all');
    slc = zeros(32,32);
    slc(bothLungs==1) = lungZ / totalZ; % proportion of total lung impedance represented by each pixel;
    
    zLV = sum(slc(ltVentMask), 'all') ./ sum(slc(ltVentMask) > 0, 'all');
    zLD = sum(slc(ltDorsMask), 'all') ./ sum(slc(ltDorsMask) > 0, 'all');
    zRV = sum(slc(rtVentMask), 'all') ./ sum(slc(rtVentMask) > 0, 'all');
    zRD = sum(slc(rtDorsMask), 'all') ./ sum(slc(rtDorsMask) > 0, 'all');
    
%     zLV         = sum(slc(ltVentMask), 'all') ./ length(ltVentMask);
%     zLD         = sum(slc(ltDorsMask), 'all') ./ length(ltDorsMask);
%     zRV         = sum(slc(rtVentMask), 'all') ./ length(rtVentMask);
%     zRD         = sum(slc(rtDorsMask), 'all') ./ length(rtDorsMask);
    
    zTotal      = zLV + zLD + zRV + zRD;
    ltVent(i)   = (zLV / zTotal);
    ltDors(i)   = (zLD / zTotal);
    rtVent(i)   = (zRV / zTotal);
    rtDors(i)   = (zRD / zTotal);
    
    rowMeans    = calc_norm_lung_row_means(slc, roi);
    if ~exist('rowMeanRes', 'var')
        rowMeanRes = zeros(NPOINTS, length(rowMeans));
    end % end if
    rowMeanRes(i,:) = rowMeans;
end % end for

perfDist.ltVent = ltVent;
perfDist.ltDors = ltDors;
perfDist.rtVent = rtVent;
perfDist.rtDors = rtDors;

for i=1:size(rowMeanRes, 2)
    fn = sprintf('rowMeanRes%.0f', i);
    perfDist.(fn) = rowMeanRes(:,i);
end
    
end % end function