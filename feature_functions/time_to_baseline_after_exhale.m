function ttb = time_to_baseline_after_exhale(triads, fs)
% -------------------------------------------------------------------------
% DESCRIPTION:
%   tau = time_to_baseline_after_exhale(imgs, triads, fs)
% 
%   Calculate the time for impedance signal to reach breath baseline from
%   end_inhalation landmark
% -------------------------------------------------------------------------
% PARAMETERS:
%   triad:
%       end_ex end_in end_ex breath landmarks for data
%   fs:
%       sampling frequency
% -------------------------------------------------------------------------   
% RETURNS:
%   tau:
%       array of time constants for each breath in data specified in triad
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
if size(triads, 2) == 3
    end_in = triads(:,2);
    end_ex = triads(:,3);
else
    error('triad must contain 3 landmarks for each breath')
end % end if

ttb = (end_ex - end_in) ./ fs;

end % end function