function breathZ = calc_breath_delta_z(imgs, triad)
% -------------------------------------------------------------------------
% DESCRIPTION:
%
% -------------------------------------------------------------------------
% PARAMETERS:
%   imgs (n x n x m array)
%       full array of imgs from EIT time series
%   triad (1 x 3 array)
%       a single breath triad (end_ex, end_in, end_ex)
% -------------------------------------------------------------------------   
% RETURNS:
%   breathZ (array):
%       impedance change of breath from EELI to EILI.
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

end_ex  = imgs(:,:,triad(1));
end_in  = imgs(:,:,triad(2));
breathZ = end_in - end_ex;

% end_ex2 = imgs(:,:,triad(3));
% end_ex  = mean(end_ex1, end_ex2, 3);

end % end function