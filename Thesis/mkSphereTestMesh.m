function [fmdl, img]= mkSphereTestMesh(maxsz, maxh)
% This function was written by :
%                               Mark Campbell
%                               Carleton University    
% Make spherical model
cd 'C:\Users\Mark\Documents\GraduateStudies\LAB\Thesis\Models\SphereTest';
if exist('sphereTestImg.mat', 'file') == 2
    mdl= load('sphereTestImg.mat');
    img= mdl.img;
    mdl= load('sphereTestFmdl.mat');
    fmdl= mdl.scalpSkullBrain;
else
    numElec= 32;
    volFilename1= '80msh';
    volFilename2= '75msh';
    volFilename3= '70msh';
    centres= [];
    stim_pattern= [];
    z_contact= 0.01;
    e_rad= 1; 
    e_height= 0;
    radius= 0.25; % - requested weighting matrix  (recommend 0.2 for 16 electrodes)
    weight= 1; % - weighting matrix (weighting of noise vs signal). Can be empty options.noise_figure is specified
%     opt.imgsz= imgsize;
    elTheta= [ -92 -74 -92 -60 -32 -72 -92 -46 -32 -72 -92 -60 46  -74 -92 92 ...
                92  74  60  92  72  32  46  92  72  32  60  92  74  92  46  0 ];        
    elPhi=   [ -72 -65 -36 -51 -45 -21 0   0   45  21  36  51  -90 65  72  -90 ...
               -72 -65 -51 -36 -21 -45 0   0   21  45  51  36  65  72  90  0 ];    
    elPhi= elPhi+ 180;
    
    [scalp, mat_indices1]= gmsh_mk_fwd_model(volFilename1, centres, stim_pattern, z_contact);
    [skull, mat_indices2]= gmsh_mk_fwd_model(volFilename2, centres, stim_pattern, z_contact);
    [brain, mat_indices3]= gmsh_mk_fwd_model(volFilename3, centres, stim_pattern, z_contact);
    
    ctr = mean(scalp.nodes(scalp.boundary,:));
    scalp.nodes= scalp.nodes- ctr; % center the model
    skull.nodes= skull.nodes- ctr; % center the model
    brain.nodes= brain.nodes- ctr; % center the model
    T= [1 1 -1];
    scalp.nodes= scalp.nodes .* T; %flip z
    skull.nodes= skull.nodes .* T; %flip z
    brain.nodes= brain.nodes .* T; %flip z
    scalp.Nz=  [0   77  -6];
    scalp.Iz=  [0   -77 -30 ];
    scalp.PAR= [77  0   -21 ];
    scalp.PAL= [-77 0   -21 ];
    elecPos= getElecPos(scalp, elTheta, elPhi); % find electrode positions
    mdl1WithElec= place_elec_on_surf(scalp, elecPos, [e_rad, e_height, maxsz], [], maxh);
    skullBrain= merge_meshes(skull, brain);
    fmdl= merge_meshes(mdl1WithElec, skullBrain);
    img = mk_image(fmdl, 1); % Background conductivity is scalp
    img.elem_data([fmdl.mat_idx{1}]) = 0.41;    %   1: scalp       
    img.elem_data([fmdl.mat_idx{2}]) = 0.016;   %   2: skull      
    img.elem_data([fmdl.mat_idx{3}]) = 0.47;    %   6: grey matter
    save('sphereTestImg.mat','img');
    save('sphereTestFmdl.mat','scalpSkullBrain');
end % end if
end % end function