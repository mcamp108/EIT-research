function NV_startup()

% -------------------------------------------------------------------------
% DESCRIPTION:
%
%   add paths for NV work
%
% -------------------------------------------------------------------------
% PARAMETERS:
% 
%
%
% -------------------------------------------------------------------------   
% RETURNS:
% 
%
%
% -------------------------------------------------------------------------   
% AUTHOR:
%   Mark Campbell
%   Carleton University
%   markacampbell@cmail.carleton.ca
% -------------------------------------------------------------------------
%   01.Oct.2019

if ~exist(did_NV_startup, 'var')
    addpath('C:\Users\Mark\Documents\Neurovine\Scripts\data-analysis');
    addpath('C:\Users\Mark\Documents\Neurovine\data\Mark\2019-09-26');
    addpath('C:\Users\Mark\Documents\Neurovine\data\Marwan\2019-09-26');
    addpath('C:\Users\Mark\Documents\Neurovine\Scripts\data-analysis');
    addpath('C:\Users\Mark\Documents\Neurovine\data\Mark\2019-09-18');
    addpath('C:\Users\Mark\Documents\Neurovine\data\Mark\2019-09-26');
    addpath('C:\Users\Mark\Documents\GraduateStudies\LAB\Scripts');
end % end if

end % end function