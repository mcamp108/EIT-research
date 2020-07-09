function [fmdl, imdl] = get_pig_mdl(pig)
% pig = '12-2';
switch pig
    case '8-2'
        fmdl_file= '8-2_fmdl.mat';
        imdl_file= '8-2_imdl.mat';
        elec_z_plane= 106; % average of reference dot z planes from elec_loc_refs
    case '9-2'
        fmdl_file= '9-2_fmdl.mat';
        imdl_file= '9-2_imdl.mat';
        elec_z_plane= 106; % average of reference dot z planes from elec_loc_refs
    case '10-2'
        fmdl_file= '10-2_fmdl.mat';
        imdl_file= '10-2_imdl.mat';
        elec_z_plane= 97; % average of reference dot z planes from elec_loc_refs
    case '11-2'
        pig = '11-2';
        fmdl_file= '11-2_fmdl.mat';
        imdl_file= '11-2_imdl.mat';
        elec_z_plane= 96; % average of reference dot z planes from elec_loc_refs
    case '12-2'
        pig = '12-2';
        fmdl_file= '12-2_fmdl.mat';
        imdl_file= '12-2_imdl.mat';
        elec_z_plane= 116; % average of reference dot z planes from elec_loc_refs
end % end switch

cd(horzcat('C:\Users\Mark\Documents\GraduateStudies\LAB\HamburgBrain\Models\', pig,'\mesh'));

if exist(fmdl_file, 'file') == 2
    fmdl= load(fmdl_file);
    fmdl= fmdl.fmdl;
    imdl= load(imdl_file);
    imdl= imdl.imdl;
else
    seg3dOutMatFile = horzcat(pig, '_seg3D_out');
    inrFileName     = horzcat(pig,'_SegImg.inr');
    elecFileName    = horzcat(pig,'_Seg_elec.txt');
    meshFileName    = horzcat(pig,'_Mesh');
    fmdlFileName    = horzcat(pig,'_fmdl');
    imdlFileName    = horzcat(pig,'_imdl');
    nElec = 32;

    load(seg3dOutMatFile);
    V = scirunnrrd.data;
    V = rot90(V, 1);    % looking through front of pig (anatomical position)      

    pixelScale = 1;
    volRes = 1 / pixelScale;
    xyzRes = [1 1 1] * volRes;

    elecPlaneImg = squeeze( V(:, :, elec_z_plane) );
    clf(); imagesc(elecPlaneImg); axis equal
    keyboard;
    rad = -linspace(-2*pi, 0, nElec+1); % first elec is bottom of pig, then proceed clockwise
    elec_pos = rad(1:32)';
    elecCoors = find_elec_locations(elecPlaneImg, nElec, elec_pos, elec_z_plane);

    % save electrode locations to file
    writetable( array2table(elecCoors), elecFileName, 'WriteVariableNames', 0 );
    % save the segmentation to inr for the mesher
    saveinr_EIT(uint8(V), inrFileName, xyzRes);

    % Run the mesher
    P = getmesherparam;

    P.facet_distance_mm = 0.4;
    P.cell_fine_size_mm = 0.5;
    P.cell_coarse_size_mm = 4;

    P.refine.electrodes = 1;
    P.electrode_radius_mm = 3;
    P.cell_size_electrodes_mm = .5;

    % Turn off all optimisations
    P.opt.exude_opt = 0;
    P.opt.lloyd_opt = 0;
    P.opt.odt_opt = 0;
    P.opt.perturb_opt = 0;

    % save the output to csv to load into matlab
    P.save.save_nodes_tetra = 1;
    % save VTK to view it paraview
    P.save.save_vtk = 0;

    % Move the electrode positions to ensure they are placed on the surface,
    % the positions are not exact in the first place
    P.move.move_electrodes = 1;
    P.move.outermost_tissue = 1;

    % write parameter file
    writemesherparam('mesh_params.txt',P);

    runmesher(inrFileName, elecFileName, 'mesh_params.txt','output/',meshFileName);

    % Load into EIDORS and check electrodes
    Mesh = loadmesh(sprintf('output/%s',meshFileName));
    %remember to run eidors startup
    %based on mk_fmdl_from_nodes
    MDL=eidors_obj('fwd_model','test');
    MDL.nodes = Mesh.Nodes;
    MDL.elems = Mesh.Tetra;
    [srf, idx] = find_boundary(Mesh.Tetra);
    MDL.boundary = srf;

    %find ground node
    gnd_dist=sqrt(sum((Mesh.Nodes - Mesh.gnd_pos).^2,2));
    [~, gnd_node] = min(gnd_dist);
    MDL.gnd_node = gnd_node;

    % format mat_idx
    MDL.mat_idx = cell( 1, length( unique(Mesh.mat_ref) ) );
    for i=1:size(MDL.mat_idx, 2)
       MDL.mat_idx{i} = find(Mesh.mat_ref == i);
    end % end for

    % Electrodes
    % find nodes on the suface within electrode radius
    z_contact = 0.01;
    elec_radius = .002; %in meters

    srf_idx = unique(srf(:));
    srf_nodes = Mesh.Nodes(srf_idx,:);

    for iElec = 1:size(Mesh.elec_pos,1)
        electrodes(iElec).z_contact= z_contact;
        edist = sqrt(sum((srf_nodes - Mesh.elec_pos(iElec,:)).^2,2));
        enodes = srf_idx(edist <= elec_radius);

        curNodes = enodes';
        curNodes = unique(curNodes);

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
    
    % TRANSFORM
    x = fmdl.nodes(:,1);
    y = fmdl.nodes(:,2);
    fmdl.nodes(:,1) = y;
    fmdl.nodes(:,2) = -x; % rotate x/y for EIDORS
    
    % center model
    ctr = mean(fmdl.nodes);
    fmdl.nodes= fmdl.nodes - ctr; % center the model
    
    figure(1);
    sgtitle(sprintf('Tissue Segmentation Meshes for %s',pig));
    subplot(4,5,1:16);
    show_fem(fmdl, [0 1 0]); view(0,45);
    subplot(4,4,17); wireframe(fmdl,1); axis equal;
    subplot(4,4,18); wireframe(fmdl,2); axis equal;
    subplot(4,4,19); wireframe(fmdl,3); axis equal;
    subplot(4,4,20); wireframe(fmdl,4); axis equal;
    % Check FWD is ok
    if valid_fwd_model(fmdl)
        disp('Forward model is ok!');
    end
    
    keyboard;
    saveas(gcf,sprintf('%s model.svg',pig));

    % Add conductivity and solve

    img = mk_image(fmdl, 0.41); % Background conductivity is scalp
    img.elem_data([fmdl.mat_idx{1}]) = 0.41;    %   1: scalp        0.41
    img.elem_data([fmdl.mat_idx{2}]) = 0.016;   %   2: skull        0.016
    img.elem_data([fmdl.mat_idx{3}]) = 0.47;    %   3: grey matter  0.47
    img.elem_data([fmdl.mat_idx{4}]) = 0.0001;  %   4: air          0.0001

    % Set stim patterns
    imgsize = [64 64];
    radius= 0.2; % - requested weighting matrix  (recommend 0.2 for 16 electrodes)
    weight= 1; % - weighting matrix (weighting of noise vs signal). Can be empty options.noise_figure is specified
    opt.noise_figure = [];
    opt.imgsz= imgsize;
    opt.keep_intermediate_results= true;
    opt.save_memory = 1;
    img.fwd_model.normalize_measurements = 0;

    imdl = mk_GREIT_model(img, radius, weight, opt);

    save(fmdlFileName,'fmdl');
    save(imdlFileName,'imdl');
end % end if
end % end function
%%
function elecCoors = find_elec_locations(elecPlaneImg, numElec, elec_pos, zPlane)

elecPlaneImg(elecPlaneImg~=0) = 1;
cY = mean(find(sum(elecPlaneImg,2) ~= 0));
cX = mean(find(sum(elecPlaneImg,1) ~= 0));
mdlCenter = [cX, cY];
elecCoors = zeros(numElec, 3);
elecCoors(:,3) = zPlane;

for i =1:length(elec_pos)
    iter = 0;
    theta = elec_pos(i);
    vx = cos(theta);
    vy = sin(theta);
    newloc = round(mdlCenter);
    while elecPlaneImg( newloc(1),newloc(2) ) ~=0
        iter = iter + 1;
        trueloc = newloc;
        newloc = round(mdlCenter + [vx vy] * iter );
    end % end while
    elecCoors(i,[1,2]) = trueloc;
end % end for

end % end function

function wireframe(mesh, material)
nodeIdx = mesh.elems( mesh.mat_idx{material}, : );
hold on;
h = trimesh(nodeIdx, mesh.nodes(:,1), mesh.nodes(:,2), mesh.nodes(:,3));
h.FaceAlpha = 0;
hold off;
end % end function

function scat3(mesh, material)
nodes = mesh.nodes( mesh.elems( mesh.mat_idx{material},: ), : );
hold on;
scatter3(nodes(:,1), nodes(:,2), nodes(:,3));
hold off;
end % end function


% rotate 90 about z axis
% x = fmdl.nodes(:,1);
% y = fmdl.nodes(:,2);
% fmdl.nodes(:,1) = y;
% fmdl.nodes(:,2) = -x;
% ctr = mean(fmdl.nodes(fmdl.boundary,:));
% fmdl.nodes= fmdl.nodes - ctr; % center the model
