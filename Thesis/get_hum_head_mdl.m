function [fmdl, imdl] = get_hum_head_mdl()
% -------------------------------------------------------------------------
% DESCRIPTION:
% [fmdl, imdl] = get_hum_head_mdl()
% -------------------------------------------------------------------------
% PARAMETERS:
% -------------------------------------------------------------------------   
% RETURNS:
% fmdl:
% imdl:
% -------------------------------------------------------------------------   
% AUTHOR:
%   Mark Campbell
%   Carleton University
%   markacampbell@cmail.carleton.ca
% -------------------------------------------------------------------------
% VERSION:
%   1.0.0
% -------------------------------------------------------------------------
fmdl_file= 'humHeadFmdl.mat';
imdl_file= 'humHeadImdl.mat';
dataDir = 'C:\Users\Mark\Documents\GraduateStudies\LAB\EIT-neuroimaging\Models\HumanHead\mesh\';
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

    elPhi =   [ -72 -65 -36 -51 -45 -21 0   0   45  21  36  51  -90 65  72  -90 ...
               -72 -65 -51 -36 -21 -45 0   0   21  45  51  36  65  72  90  0 ];    
    elPhi= elPhi + 180;
        
    seg3dOutMatFile = 'humHeadSeg';
    inrFileName     = 'humHeadSegImg.inr';
    elecFileName    = 'humHeadSeg_elec.txt';
    meshFileName    = 'humHeadMesh';
    fmdlFileName    = 'humHeadFmdl';
    imdlFileName    = 'humHeadImdl';
%     nElec           = 32;

    % load segmentation
    I           = load(horzcat(dataDir,seg3dOutMatFile));
    I           = I.scirunnrrd;
    V.data      = I.data;
    pixelScale  = [I.axis(1).spacing, I.axis(2).spacing, I.axis(3).spacing];
    xyzRes      = 1 ./ pixelScale;
    
    % find electrode locations
    V.Nz  = [80  3   24];
    V.Iz  = [80  195 11 ];
    V.PAR = [5   117 11 ];
    V.PAL = [156 117 11 ];
    elecPos = getElecPos(V, elTheta, elPhi); % find electrode positions
    writetable( array2table(elecPos), elecFileName, 'WriteVariableNames', 0 ); % save electrode locations to file
    saveinr_EIT(uint8(V.data), inrFileName, xyzRes); % save the segmentation to inr for the mesher

    % Run the mesher
    factor = 1;
    P = getmesherparam;
    
    % facet sizes
    P.facet_distance_mm     = 0.2 * factor;
    P.cell_fine_size_mm     = 2 * factor;
    P.cell_coarse_size_mm   = 4 * factor;
    
    % electrode size
    P.refine.electrodes         = 1;
    P.electrode_radius_mm       = 10 * factor;
    P.cell_size_electrodes_mm   = .5 * factor;
    
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
    
    % TRANSFORM
    x = fmdl.nodes(:,1);
    y = fmdl.nodes(:,2);
    fmdl.nodes(:,1) = y;
    fmdl.nodes(:,2) = -x; % rotate x/y for EIDORS
    
    % center model
    ctr = mean(fmdl.nodes);
    fmdl.nodes= fmdl.nodes - ctr; % center the model

    figure(1);
    sgtitle('Tissue Segmentation Meshes');
    subplot(4,4,[1,2,5,6,]);    show_fem(fmdl, [0 1 0]); view(-90,45); axis equal; title('FEM: Front');
    subplot(4,4,[9,10,13,14,]); show_fem(fmdl, [0 1 0]); view(90,45);  axis equal; title('FEM: Back');
    subplot(4,4,3);  wireframe(fmdl,1); axis equal; view(-30,30); title('Scalp');
    subplot(4,4,4);  wireframe(fmdl,2); axis equal; view(-30,30); title('Skull');
    subplot(4,4,7);  wireframe(fmdl,3); axis equal; view(-30,30); title('CSF');
    subplot(4,4,8);  wireframe(fmdl,4); axis equal; view(-30,30); title('Grey Matter');
    subplot(4,4,11); wireframe(fmdl,5); axis equal; view(-30,30); title('White Matter');
    subplot(4,4,12); wireframe(fmdl,6); axis equal; view(-30,30); title('Diploe');
    subplot(4,4,15); wireframe(fmdl,7); axis equal; view(-30,30); title('Air');

    % Check FWD is ok
    if valid_fwd_model(fmdl)
        disp('Forward model is ok!');
    end
    keyboard;
    saveas(gcf,'C:\Users\Mark\Documents\GraduateStudies\LAB\EIT-neuroimaging\Figures\Human model.svg');

    % Add conductivity and solve

    img = mk_image(fmdl, 0.41); % Background conductivity is scalp
    img.elem_data([fmdl.mat_idx{1}]) = 0.41;    %   1: scalp
    img.elem_data([fmdl.mat_idx{2}]) = 0.016;   %   2: skull
    img.elem_data([fmdl.mat_idx{3}]) = 1.71;    %   3: CSF
    img.elem_data([fmdl.mat_idx{4}]) = 0.47;    %   4: grey matter
    img.elem_data([fmdl.mat_idx{5}]) = 0.22;    %   5: white matter
    img.elem_data([fmdl.mat_idx{6}]) = 0.7;     %   6: diploe
    img.elem_data([fmdl.mat_idx{7}]) = 0.0001;  %   7: air

    % Set stim patterns
    imgsize = [64 64];
    radius = 0.2; % - requested weighting matrix  (recommend 0.2 for 16 electrodes)
    weight = 1; % - weighting matrix (weighting of noise vs signal). Can be empty options.noise_figure is specified
    opt.noise_figure = [];
    opt.imgsz = imgsize;
    opt.keep_intermediate_results = true;
    img.fwd_model.normalize_measurements = 0;
    imdl = mk_GREIT_model(img, radius, weight, opt);
    save(fmdlFileName,'fmdl');
    save(imdlFileName,'imdl');
end % end if

%%
% clf();
% wireframe(fmdl,4); axis equal; view(-30,30); title('Scalp');
%%

end % end function



% ======================================================================= %
%                              SUBFUNCTIONS
% ======================================================================= %

function wireframe(mesh, material)
nodeIdx = mesh.elems( mesh.mat_idx{material}, : );
hold on;
h = trimesh(nodeIdx, mesh.nodes(:,1), mesh.nodes(:,2), mesh.nodes(:,3));
h.FaceAlpha = 0;
hold off;
end % end function

% ======================================================================= %

function scat3(mesh, material)
nodes = mesh.nodes( mesh.elems( mesh.mat_idx{material},: ), : );
hold on;
scatter3(nodes(:,1), nodes(:,2), nodes(:,3));
hold off;
end % end function

% ======================================================================= %

function elecPos = getElecPos(mdl1, elTheta, elPhi)
% Find electrode xyz coordinates from anatomical landmarks (must be fields
% in the model)

% find the polar coordinates of origin and radius of imaginary sphere from 
% nasion and preauriculars landmarks
alph= 65; % 25 degree angle between nasion-inion line and level plane
Nz  = mdl1.Nz;
Iz  = mdl1.Iz;
PAL = mdl1.PAL;
PAR = mdl1.PAR;
ap  = (Iz-Nz)/norm(Iz-Nz); % anterior/posterior unit vector
apM = (Nz+ Iz)/ 2; % midpoint between nasion and inion
apL = norm(Nz- Iz); % distance between nasion and inion
lr  = (PAL-PAR)/norm(PAL-PAR); % left/right unit vector
si  = cross(lr, ap); % superior/ inferior unit vector
Or  = apM + (apL* si)/ (2*tand(alph)); % location of origin of imaginary sphere
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
% mdl = prepare_surf_model(mdl1);

% % Find electrode xyz coordinates
% V1= mdl.nodes(mdl.boundary(:, 1), :);
% V2= mdl.nodes(mdl.boundary(:, 2), :);
% V3= mdl.nodes(mdl.boundary(:, 3), :);
% elecPos= zeros(32, 3);
% for e= 1:32
%     elV= elecPosABC(e,:);
%     elecPos(e,:)= getEPos(elV, Or, mdl, V1, V2, V3);
% end % end for
% figure; show_fem(mdl); hold on
% % plot3([Or(1) elecPosABC(1,1)], [Or(2) elecPosABC(1,2)], [Or(3) elecPosABC(1,3)]);
% % plot3(Vs(:,1), Vs(:,2), Vs(:,3));
% plot3(elecPos(:,1), elecPos(:,2), elecPos(:,3), 'o'); hold off
% keyboard;
end % end function

% ======================================================================= %

function [I,check]=plane_line_intersect(n,V0,P0,P1)
%plane_line_intersect computes the intersection of a plane and a segment(or
%a straight line)
% Inputs: 
%       n: normal vector of the Plane 
%       V0: any point that belongs to the Plane 
%       P0: end point 1 of the segment P0P1
%       P1:  end point 2 of the segment P0P1
%
%Outputs:
%      I    is the point of interection 
%     Check is an indicator:
%      0 => disjoint (no intersection)
%      1 => the plane intersects P0P1 in the unique point I
%      2 => the segment lies in the plane
%      3=>the intersection lies outside the segment P0P1
%
% Example:
% Determine the intersection of following the plane x+y+z+3=0 with the segment P0P1:
% The plane is represented by the normal vector n=[1 1 1]
% and an arbitrary point that lies on the plane, ex: V0=[1 1 -5]
% The segment is represented by the following two points
% P0=[-5 1 -1]
%P1=[1 2 3]   
% [I,check]=plane_line_intersect([1 1 1],[1 1 -5],[-5 1 -1],[1 2 3]);
%This function is written by :
%                             Nassim Khaled
%                             Wayne State University
%                             Research Assistant and Phd candidate
%If you have any comments or face any problems, please feel free to leave
%your comments and i will try to reply to you as fast as possible.

I=[0 0 0];
u = P1-P0;
w = P0 - V0;
D = dot(n,u);
N = -dot(n,w);
check=0;
if abs(D) < 10^-7        % The segment is parallel to plane
        if N == 0           % The segment lies in plane
            check=2;
            return
        else
            check=0;       %no intersection
            return
        end
end
%compute the intersection parameter
sI = N / D;
I = P0+ sI.*u;
if (sI < 0 || sI > 1)
    check= 3;          %The intersection point  lies outside the segment, so there is no intersection
else
    check=1;
end % end if

end % end function

% ======================================================================= %

function ePos= getEPos(p, Or, mdl, V1, V2, V3)

fc= findClosestFace(p, mdl, V1, V2, V3); 
if length(fc)> 1
    keyboard;
end % end if
fcVert= mdl.boundary(fc,:);
v1= mdl.nodes(fcVert(1),:);
v2= mdl.nodes(fcVert(2),:);
v3= mdl.nodes(fcVert(3),:);
n= cross(v1- v2, v1-v3);
ePos= plane_line_intersect(n, v1, p, Or);

% double check we had the correct face
fc2= findClosestFace(ePos, mdl, V1, V2, V3);
if fc~= fc2
    ePos= getEPos(ePos, Or, mdl, V1, V2, V3);
end % end if

end % end function

% ======================================================================= %

function fc= findClosestFace(p, mdl, V1, V2, V3)

p0= repmat(p, size(mdl.boundary, 1), 1);
dif= sum([(V1- p0).^2, (V2- p0).^2, (V3- p0).^2], 2); % find the closest face
fc= find(dif==min(dif));

end % end function

% ======================================================================= %