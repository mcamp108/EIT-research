function [fmdl, imdl]= mk_weighted_restraint_model()

% -------------------------------------------------------------------------
% DESCRIPTION:
%   [fmdl, imdl]= mk_weighted_restraint_model()
% -------------------------------------------------------------------------
% PARAMETERS:
% -------------------------------------------------------------------------   
% RETURNS:
%   fmdl:
%       EIDORS forward model of human thorax with 32 electrodes arranged in
%       square pattern
%   imdl:
%       EIDORS inverse model of fmdl with keep_intermediate_results= true.
% -------------------------------------------------------------------------   
% AUTHOR:
%   Mark Campbell
%   Carleton University
%   markacampbell@cmail.carleton.ca
% -------------------------------------------------------------------------
% VERSION:
%   1.0.0
% -------------------------------------------------------------------------

starting_dir= cd;
cd 'C:\Users\Mark\Documents\GraduateStudies\LAB\WeightedRestraint\Models';

if exist('weighted_restraint_fmdl.mat', 'file') == 2
    mdl= load('weighted_restraint_fmdl.mat');
    fmdl= mdl.fmdl;
    mdl= load('weighted_restraint_imdl.mat');
    imdl= mdl.imdl;
else
    thorax = shape_library('get','adult_male','boundary');
    shape = { 2, {thorax}, [4,40], 0.5};

    % Set up the electrodes with 1:16 on bottom and 17:32 on top
    elec_pos = [16,0,0.75,1.25];
    elec_spec = [0.05];
    % elec_spec = [0.5, 0, 0.1]; % Radius of circular electrodes
    fmdl = ng_mk_extruded_model(shape, elec_pos, elec_spec);
    row1= [13:16, 1:12]; 
    row2= [29:32, 17:28]; 
    idx= [row1;row2]';
    fmdl.electrode(idx)= fmdl.electrode(:);

    % Set up square pattern
    idx = reshape([32, 1:31],2,[])'; %matrix of 32 elem numbered from 1:32, shaped with row length= L then transposed so numbering goes across rows
    idx(2:2:end,:) = fliplr(idx(2:2:end,:)); % flips every second row to get proper square electrode arrangement
    fmdl.electrode(idx) = fmdl.electrode(:);

    % Stim pattern
    [fmdl.stimulation, fmdl.meas_select]=mk_stim_patterns(32,1,[0,5],[0,5],{'no_meas_current_next2'},1); % Skip 4 
    fmdl.normalize = 0;

    % imdl
    vopt.imgsz = [32 32];
    vopt.square_pixels = true;
    vopt.zvec = linspace(0.5,1.5,4);
    vopt.save_memory = 1;
    opt.noise_figure = 1.0;
    opt.keep_intermediate_results= true;
    [imdl_t, opt.distr] = GREIT3D_distribution(fmdl, vopt);
    imdl= mk_GREIT_model(imdl_t, 0.2, [], opt);
    
    save('weighted_restraint_fmdl.mat','fmdl');
    save('weighted_restraint_imdl.mat','imdl');
end % end if

cd(starting_dir);

end % end function