function Dxy = localization_err(img1, img2)
% -------------------------------------------------------------------------
% DESCRIPTION:
%   [x,y] = localization_acc(img1, img2)
%   Ayati 2015
%
% -------------------------------------------------------------------------
% PARAMETERS:
%   img1:
%
%   img2:
% -------------------------------------------------------------------------   
% RETURNS:
%   x:
%
%   y:
% -------------------------------------------------------------------------   
% AUTHOR:
%   Mark Campbell
%   Carleton University
%   markacampbell@cmail.carleton.ca
% -------------------------------------------------------------------------
% VERSION:
%   1.0.0
% -------------------------------------------------------------------------

[x1,y1] = weighted_sum(img1);
[x2,y2] = weighted_sum(img2);

Dxy = sqrt((x2 - x1)^2 + (y2 - y1)^2);

% divide Dxy by image radius / dimensions (and mismatch between image
% sizes?) Need this because dimensions of each image are not the same

end % end function

function [x,y] = weighted_sum(img)
    szX     = size(img, 2);
    szY     = size(img, 1);
    sumX    = sum(img .* (1:szX), 'all');
    sumY    = sum(img .* (1:szY)', 'all');
    x       = (sumX / sum(img(:))) / szX;
    y       = (sumY / sum(img(:))) / szY;
end % end function
