function imdl=jan12_GREIT_cyl(L); % L is layers

    fmdl= ng_mk_cyl_models([4,1.0,.5],[32/L,linspace(1.6,2.4,L)],[0.05]);
    % forward model = netgen make cylindrical model. height of 4, radius of
    % 1, 0.5 = max size of mesh elements. Second argument is electrode
    % position. 32 electrodes in 4 L layers (in my case 4), evenly spaced 
    %in z plane between positions z= 1.6 -> z= 2.4. electrode radius= 0.05
    %m
    fmdl = jan12_stims_rearrange(fmdl,L);
    %
    opt = struct('imgsz', [32 32],'square_pixels', true);
    % create 32x32 matrix. vopt= volume optimization
    vopt= opt; vopt.save_memory = 1;
    %  NOTE: mk_voxel_volume is slow and uses lots of memory To reduce the 
    %memory footprint, try vopt.save_memory = 1; %(or try 10)
    vopt.zvec = linspace(-1,1,10)*1.125+2;
    %add z vector, weird values though...
    [imdl,opt.distr] = GREIT3D_distribution(fmdl, vopt);
    opt.noise_figure = 1.0;
    imdl= mk_GREIT_model(imdl, 0.20, [], opt);
