function [fmdl, imdl] = mk_horse_model(shiftElecs)
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

mdl = load('C:\Users\Mark\Documents\EIT\uMontreal\recon\horse.mat');
fmdl = mdl.mdl;
if nargin == 1
    fmdl = shift_electrodes(fmdl, shiftElecs);
end % end if

% stim pattern
[fmdl.stimulation, fmdl.meas_select] = mk_stim_patterns(32,1,[0,5],[0,5],{'no_meas_current_next2'}, 0.3);

% adjust element conductivity
img = mk_image(fmdl, 1);

% img.elem_data([fmdl.mat_idx{2}; fmdl.mat_idx{3}]) = 0.3;

img.fwd_model.normalize_measurements = 0;

radius  = 0.2;
opt.noise_figure = 0.5;
opt.imgsz = [32 32];
opt.keep_intermediate_results = true;
opt.square_pixels = true;
opt.save_memory = 1;
imdl = mk_GREIT_model(img, radius, [], opt);

end % end function