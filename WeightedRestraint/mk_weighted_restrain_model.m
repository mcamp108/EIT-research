function mdl= mk_weighted_restrain_model(volfile)

% -------------------------------------------------------------------------
% DESCRIPTION:
%
%
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
% VERSION:
%   1.0.0
% -------------------------------------------------------------------------

cd 'C:\Users\Mark\Documents\GraduateStudies\LAB\WeightedRestraint\Models';
numElec= 32;
volFilename1= 'wholeSkullAllSegMsh';
name1= "HumanHead";

stim_pattern= [];
z_contact= 0.01;
e_rad= 1; 
e_height= 0;
radius= 0.25; % - requested weighting matrix  (recommend 0.2 for 16 electrodes)
weight= 1; % - weighting matrix (weighting of noise vs signal). Can be empty options.noise_figure is specified
opt.imgsz= imgsize;
volfilename= 'total_w_electrode';


end % end function