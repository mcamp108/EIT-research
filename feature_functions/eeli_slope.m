function slope = eeli_slope(eitFile)
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
imgs = calc_slices(eitFile.imgr);
triads = eitFile.triads;
fs = eitFile.fs;
Zimgs = -imgs;
roi = horse_roi();
eeliLocs = [triads(:, 1); triads(end, 3)];
eeliSums = nan( length(eeliLocs), 1);
for i = 1 : length(eeliLocs)
    img = squeeze(Zimgs(:, :, eeliLocs(i)));
    eeliSums(i) = sum(img(roi.BothLungs), 'all');
end

L1 = polyfit(eeliLocs, eeliSums, 1);
slope = L1(1) * fs * 60;
        
end % end function