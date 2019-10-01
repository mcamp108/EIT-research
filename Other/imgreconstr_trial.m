
nelecs= 32;
nrings= 4;
elec_per_ring= 8;
forearm_length= 25;
av_ring_rad= 4.02;
num_electrode_rings= 4;
elec_elec_dist= 5;
% load EIT data.
[vv, auxdata, stim]= eidors_readdata('StandingArmLoweredEx-30PU.eit', 'LQ4');
% first column of vv is 32x32= 1024 rows. second column is related to
% number of frames.
% forward model using measured forearm length and electrode positions
fmdl= ng_mk_cyl_models([forearm_length, av_ring_rad], ... 
    [elec_per_ring, linspace(5, 20, num_electrode_rings)], 0.05);
% fmdl with rectangular electrode arrangement
fmdl = jan12_stims_rearrange(fmdl, nrings);
% create inverse model
imdl= mk_common_model('b3cr', [32,4]);
imdl.fwd_model= fmdl;
img = mk_image(imdl); vh = fwd_solve(img);
show_3d_slices(img);

%EVERYTHING ABOVE THIS LINE IS WORKING.

% create imgsz= 32x32 matrix with square pixels= true. 
opt = struct('imgsz', [32 32],'square_pixels', true);
%vopt= volume optimization
vopt= opt; vopt.save_memory = 1;
%  NOTE: mk_voxel_volume is slow and uses lots of memory To reduce the 
%memory footprint, try vopt.save_memory = 1; %(or try 10)
%zvec will be 7 evenly spaced planes (one at each electrode row and one in
%between).
vopt.zvec = linspace(5,20,7);

%Output: 
% imdl   - GREIT inverse model
% weight - value of the weight paramater chosed to satisfy the prescribed
% noise figure (NF). See options.noise_figure below. 
[imdl,opt.distr] = GREIT3D_distribution(fmdl, vopt);
opt.noise_figure = 1.0;
imdl= mk_GREIT_model(imdl, av_ring_rad, [], opt);

% solve
% img= inv_solve(imdl, vv, vh);

