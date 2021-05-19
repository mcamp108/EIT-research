function cov = centre_of_ventilation(imgs, triads)
% -------------------------------------------------------------------------
% DESCRIPTION:
%   cov = centre_of_ventilation(imgs, triads)
% -------------------------------------------------------------------------
% PARAMETERS:
%   img:
%       Breath difference image, where impedance is high, conductivity is
%       low. 2D array of image pixels.
%   lung_roi:
%       Logical matrix showing lung segmentation.
% -------------------------------------------------------------------------   
% RETURNS:
%   cov:
%       struct with fields:
%       x:
%           x position (image column number) of centre of ventilation
%       y:
%           y position (image row number) of centre of ventilation
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
n_breaths = size(triads, 1);
cov_x = nan(n_breaths, 1);
cov_y = nan(n_breaths, 1);
imgs = -imgs;

for i = 1 : n_breaths
    this_breath = calc_breath_delta_z(imgs, triads(i, :));
    [cov_x(i), cov_y(i)] = do_cov_calc(this_breath, roi.BothLungs);
end

cov = struct;
cov.x = cov_x;
cov.y = cov_y;

end

function [cov_x, cov_y] = do_cov_calc(img, lung_roi)
    % remove image components not related to respiration
    img(isnan(img)) = 0;        % zero NaN elements
    img_ = (img .* lung_roi);   % zero out pixels outside of lungs
    img_(img_ < 0) = 0;         % zero out pixels with no ventilation
    img_ = 100 * img_ ./ sum(img_, 'all'); % convert pixel values to % of total

    [n_rows, n_cols] = size(img_);

    % find x position
    col_sums = sum(img_, 1);
    c = cumsum(col_sums);
    cn = find(c <=50, 1, 'last');

    if isempty(cn)
        cov_x = nan;
        cov_y = nan;
        return
    end % end if

    k = (50 - sum(col_sums(1 : cn))) / cn;
    cov_x = round( ((cn + k + 0.5) / (n_cols + 1)) * n_cols);

    % find y position
    row_sums = sum(img_, 2);
    r = cumsum(row_sums);
    rn = find(r <= 50, 1, 'last');
    k = (50 - sum(row_sums(1 : rn))) / rn;
    cov_y = round( ((rn + k + 0.5) / (n_rows + 1)) * n_rows); % Y coor
end % end function