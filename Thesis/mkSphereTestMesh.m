function [fmdl, imdl] = mkSphereTestMesh()
% This function was written by :
%                               Mark Campbell
%                               Carleton University    
fmdl_file= 'sphereFmdl.mat';
imdl_file= 'sphereImdl.mat';
dataDir = 'C:\Users\Mark\Documents\GraduateStudies\LAB\EIT-neuroimaging\Models\sphere';
cd(dataDir);

if exist(horzcat(fmdl_file), 'file') == 2
    fmdl= load(fmdl_file);
    fmdl= fmdl.fmdl;
    imdl= load(imdl_file);
    imdl= imdl.imdl;
else
    % electrode angles from origin in 3D space
    elTheta = [ -92 -74 -92 -60 -32 -72 -92 -46 -32 -72 -92 -60 46  -74 -92 92 ...
                92  74  60  92  72  32  46  92  72  32  60  92  74  92  46  0 ];
    elPhi   = [ -72 -65 -36 -51 -45 -21 0   0   45  21  36  51  -90 65  72  -90 ...
               -72 -65 -51 -36 -21 -45 0   0   21  45  51  36  65  72  90  0 ];    
    elPhi= elPhi + 180;
    inrFileName     = 'sphereSegImg.inr';
    elecFileName    = 'sphereSeg_elec.txt';
    meshFileName    = 'sphereMesh';
    fmdlFileName    = 'sphereFmdl';
    imdlFileName    = 'sphereImdl';

    % Make spherical model
    DIAM    = 100;
    RAD     = DIAM/2;
    npers   = 1;
    spc     = 1/npers; 
    rlim    = RAD; % specify spacing and radius limit
    skLim   = 0.9*RAD;
    brLim   = 0.8*RAD;
    [x,y,z] = meshgrid(-rlim:spc:rlim,-rlim:spc:rlim, -rlim:spc:rlim); % from: by: to.
    img = sqrt(x.^2 + y.^2 + z.^2) < rlim;
    sk = sqrt(x.^2 + y.^2 + z.^2) < skLim;
    br = sqrt(x.^2 + y.^2 + z.^2) < brLim;
    V.data = img*1 + sk*1 + br*1;

    % find electrode locations
    f     = size(img,1);    % full model dimension
    h     = f/2;            % half model dimension
    V.Nz  = [h 1 h];
    V.Iz  = [h f h];
    V.PAR = [1 h h];
    V.PAL = [f h h];
    elecPos = getElecPos(V, elTheta, elPhi); % find electrode positions
    writetable( array2table(elecPos), elecFileName, 'WriteVariableNames', 0 ); % save electrode locations to file
    saveinr_EIT(uint8(V.data), inrFileName, [1,1,1]); % save the segmentation to inr for the mesher

    % Run the mesher
    factor = 1;
    P = getmesherparam;
    % facet sizes
    P.cell_coarse_size_mm   = 0.2*DIAM;
    P.facet_distance_mm     = P.cell_coarse_size_mm;
    P.cell_fine_size_mm     = P.cell_coarse_size_mm/8;
    
    % electrode size
    P.refine.electrodes         = 1;
    P.electrode_radius_mm       = P.cell_fine_size_mm;
    P.cell_size_electrodes_mm   = P.cell_fine_size_mm/4;
    % Turn off all optimisations
    P.opt.exude_opt     = 0;
    P.opt.lloyd_opt     = 0;
    P.opt.odt_opt       = 0;
    P.opt.perturb_opt   = 0;
    % save the output to csv to load into matlab
    P.save.save_nodes_tetra = 1;
    P.save.save_vtk         = 0; % save VTK to view it paraview   
    % Move the electrode positions to ensure they are placed on the surface,
    % the positions are not exact in the first place
    P.move.move_electrodes  = 1;
    P.move.outermost_tissue = 1;
    % write parameter file
    writemesherparam('mesh_params.txt',P);
    runmesher(inrFileName, elecFileName, 'mesh_params.txt','output/',meshFileName);
    % Load into EIDORS and check electrodes
    Mesh         = loadmesh(sprintf('output/%s',meshFileName));
    %remember to run eidors startup
    %based on mk_fmdl_from_nodes
    MDL          = eidors_obj('fwd_model','test');
    MDL.nodes    = Mesh.Nodes;
    MDL.elems    = Mesh.Tetra;
    [srf, idx]   = find_boundary(Mesh.Tetra);
    MDL.boundary = srf;
    %find ground node
    gnd_dist        = sqrt(sum((Mesh.Nodes - Mesh.gnd_pos).^2,2));
    [~, gnd_node]   = min(gnd_dist);
    MDL.gnd_node    = gnd_node;
    MDL.mat_idx     = cell( 1, length( unique(Mesh.mat_ref) ) ); % format mat_idx
    for i=1:size(MDL.mat_idx, 2)
       MDL.mat_idx{i} = find(Mesh.mat_ref == i);
    end % end for

    % Electrodes
    % find nodes on the suface within electrode radius
    z_contact   = 0.01;
    elec_radius = P.electrode_radius_mm / 1000; %in meters
    srf_idx     = unique(srf(:));
    srf_nodes   = Mesh.Nodes(srf_idx,:);

    for iElec = 1:size(Mesh.elec_pos,1)
        electrodes(iElec).z_contact = z_contact;
        edist                       = sqrt(sum((srf_nodes - Mesh.elec_pos(iElec,:)).^2,2));
        enodes                      = srf_idx(edist <= elec_radius);
        curNodes                    = enodes';
        curNodes                    = unique(curNodes);

        % only use nodes which can be used in triangulation as this confuses
        % find_electrode_bdy.m, as it thinks any unused node is a point electrode.
        % This is adapted from find_electrode_bdy

        % count how many times each node is used in the surface
        bdy_els = zeros(size(srf,1),1);
        for nd= curNodes(:)'
            bdy_els = bdy_els + any(srf==nd,2);
        end
        % Nodes used three times are part of a triangulation
        ffb = find(bdy_els == size(srf,2));
        % only use nodes which appear 3 times 
        curNodes=intersect(curNodes,unique(srf(ffb,:)));
        electrodes(iElec).nodes = curNodes;
    end

    MDL.electrode =     electrodes;
    MDL.solve =         @fwd_solve_1st_order;
    MDL.jacobian =      @jacobian_adjoint;
    MDL.system_mat =    @system_mat_1st_order;
    MDL.normalize_measurements = 0;
    [MDL.stimulation, MDL.meas_select] = mk_stim_patterns(32,1,[0,5],[0,5],{'no_meas_current_next2'},1); % try next3
    fmdl = MDL;

    % center model
    ctr = mean(fmdl.nodes);
    fmdl.nodes= fmdl.nodes - ctr; % center the model

%     figure(1);
%     sgtitle('Tissue Segmentation Meshes');
%     subplot(3,4,[1:3,5:7,9:11]);    show_fem(fmdl, [0 1 0]); axis equal; title('FEM: Front'); view(0,45);
%     subplot(3,4,4);  wireframe(fmdl,1); axis equal; view(-30,30); title('Scalp');
%     subplot(3,4,8);  wireframe(fmdl,2); axis equal; view(-30,30); title('Skull');
%     subplot(3,4,12);  wireframe(fmdl,3); axis equal; view(-30,30); title('Brain');

    % Check FWD is ok
    if valid_fwd_model(fmdl)
        disp('Forward model is ok!');
    end
%     keyboard;
%     saveas(gcf,'C:\Users\Mark\Documents\GraduateStudies\LAB\EIT-neuroimaging\Figures\Sphere model.svg');

    % Add conductivity and solve


    % Make 3D distribution
    nPix = 32;
    vopt.imgsz          = [nPix nPix nPix];
    vopt.square_pixels  = true;
%     vopt.zvec           = linspace(0, max(fmdl.nodes(:,3)), nPix + 1); % try 6 next time. targetsize= 0.17 for 10 0.20 for 6
    vopt.save_memory    = 1;
    [imdl_t, opt.distr] = GREIT3D_distribution(fmdl, vopt);

    radius = 0.2; % - requested weighting matrix  (recommend 0.2 for 16 electrodes)
    weight = 0.5; % - weighting matrix (weighting of noise vs signal). Can be empty options.noise_figure is specified
    opt.keep_intermediate_results = true;
    imdl = mk_GREIT_model(imdl_t, radius, weight, opt);
    save(fmdlFileName,'fmdl');
    save(imdlFileName,'imdl');
end % end if
    
end % end function

% ======================================================================= %
%                              SUBFUNCTIONS
% ======================================================================= %

function elecPos = getElecPos(mdl1, elTheta, elPhi)
% Find electrode xyz coordinates from anatomical landmarks (must be fields
% in the model)

% find the polar coordinates of origin and radius of imaginary sphere from 
% nasion and preauriculars landmarks
Nz  = mdl1.Nz;
Iz  = mdl1.Iz;
PAL = mdl1.PAL;
PAR = mdl1.PAR;
ap  = (Iz-Nz)/norm(Iz-Nz); % anterior/posterior unit vector
lr  = (PAL-PAR)/norm(PAL-PAR); % left/right unit vector
si  = cross(lr, ap); % superior/ inferior unit vector
Or  = size(mdl1.data)/2; % location of origin of imaginary sphere
rad = norm(Or- Nz);
T   = [lr; ap; si; Or]; % compose transformation matrix (X Y Z origin)

% find Cartesian coordinates and do transformation to RAS coordinates
elX= rad* sind(elTheta).* cosd(elPhi);
elY= rad* sind(elTheta).* sind(elPhi);
elZ= rad* cosd(elTheta);
elX1= elX.* T(1,1)+ elY* T(1,2)+ elZ* T(1,3)+ T(4,1);
elY1= elX.* T(2,1)+ elY* T(2,2)+ elZ* T(2,3)+ T(4,2);
elZ1= elX.* T(3,1)+ elY* T(3,2)+ elZ* T(3,3)+ T(4,3);
elecPos = [elX1; elY1; elZ1]';
end % end function

% ======================================================================= %

function wireframe(fmdl, material)
nodeIdx = fmdl.elems( fmdl.mat_idx{material},: );
hold on;
h = trimesh(nodeIdx, fmdl.nodes(:,1),fmdl.nodes(:,2),fmdl.nodes(:,3));
h.FaceAlpha = 0;
hold off;
end % end function
