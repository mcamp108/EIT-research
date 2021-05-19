function normMeanRowSum = calc_norm_lung_row_means(img, roi)
% -------------------------------------------------------------------------
% DESCRIPTION:
%   rowMeans = calc_norm_lung_row_means(img, roi)
% -------------------------------------------------------------------------
% PARAMETERS:
%   img (array): 
%       reconstructed image
%   roi (struct): 
%       horse ROI with lung segmentations
% -------------------------------------------------------------------------   
% RETURNS:
%   rowSums (array): 
%       row sums, length equal to number of rows in the lung segmentation.
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
bothLungs = roi.BothLungs;
img(bothLungs==0) = 0; % zero everything outside of lungs
img(img < 0) = 0;
rowSums = sum(img, 2); % total percentage in each row
lungRows = sum(bothLungs,2) ~= 0;

% rowDenom = sum(bothLungs(lungRows, :), 2);
rowDenom = sum(img > 0, 2);

meanRowSums = rowSums(lungRows) ./ rowDenom(lungRows);
meanRowSums(isnan(meanRowSums)) = 0;
normMeanRowSum = meanRowSums ./ sum(meanRowSums);
assert( round(sum(normMeanRowSum), 10) == 1);

end % end function