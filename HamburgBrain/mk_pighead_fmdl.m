function [fmdl, imdl]= mk_pighead_fmdl(maxsz, maxh, imgsize)
% This function was written by :
%                               Mark Campbell
%                               Carleton University

% cd 'C:\Users\Mark\Documents\GraduateStudies\LAB\HamburgBrain\Models\HamburgPig';
cd 'C:\Users\Mark\Documents\GraduateStudies\LAB\HamburgBrain\Models\10-2';

if exist('fmdl_10_2.mat', 'file') == 2
    fmdl= load('fmdl_10_2.mat');
    fmdl= fmdl.fmdl;
    imdl= load('imdl_10_2.mat');
    imdl= imdl.imdl;
else    
    numElec= 32;
    elec_z_plane= 116; % average of reference dot z planes from elec_loc_refs
    degrees= linspace(0, 360, numElec+1);
    elec_pos= [degrees(1:32)', repmat(elec_z_plane, numElec, 1)];
    e_rad= 1;
    e_height= 0;
%         volFilename1= 'Airchill_3D_sbde6_msh_ascii4';
%         name= 'airchill3D';
%         volFilename2= 'Airchill_3D_sbde6_noScalp_ascii2';
%         name2= 'airchillNoScalp';
    volFilename1= '10-2 simplified segmentations msh';
    name= '10-2';
    volFilename2= '10-2simplifiednoScalpmsh';
    name2= '10-2NoScalp';
    stim_pattern= [];
    z_contact= 0.01;
    radius= 0.2; % - requested weighting matrix  (recommend 0.2 for 16 electrodes)
    weight= []; % - weighting matrix (weighting of noise vs signal). Can be empty options.noise_figure is specified
    opt.noise_figure = 1.0;
    opt.imgsz= imgsize;
    opt.keep_intermediate_results= true;
    
    % Create model from all segmentations
    [mdl1, mat_indices]= gmsh_mk_fwd_model(volFilename1, name, [], stim_pattern, z_contact);
    
    % Add electrodes to boundary (removes internal segmentations)
    mdl1WithElec= place_elec_on_surf(mdl1, elec_pos, [e_rad, e_height, maxsz], [], maxh);
    
    % re-arrange so electrode 1 is at bottom of head
    temp= mdl1WithElec; 
    for i= 1: numElec
        j= i- 8;
        if j< 1
            j= j+ 32;
        end % end if
        mdl1WithElec.electrode(i).nodes= temp.electrode(j).nodes;
    end % end for
    
    % Create model from all segmentations except scalp
    [mdl2, mat_indices2]= gmsh_mk_fwd_model(volFilename2, name2, [], stim_pattern, z_contact);
    
    % Merge meshes into one model
    fmdl= merge_meshes(mdl1WithElec, mdl2);
    
    % Set conductivity parameters
    img = mk_image(fmdl, 0.41); % Background conductivity is scalp
%         img.elem_data([fmdl.mat_idx{1}]) = 0.44;    %   1: scalp        0.44 (soft tissue)
%         img.elem_data([fmdl.mat_idx{2}]) = 0.018;   %   2: skull        0.018
%         img.elem_data([fmdl.mat_idx{3}]) = 0.3;     %   3: white matter 0.3
%         img.elem_data([fmdl.mat_idx{4}]) = 0.15;    %   4: grey matter  0.15
%         img.elem_data([fmdl.mat_idx{5}]) = 0.0001;  %   5: air          0.0001
%         img.elem_data([fmdl.mat_idx{6}]) = 1.79;    %   6: CSF          1.79
%         img.elem_data([fmdl.mat_idx{7}]) = 0.7;     %   7: diploe       0.7
    img.elem_data([fmdl.mat_idx{1}]) = 0.41;    %   1: scalp        0.41
    img.elem_data([fmdl.mat_idx{2}]) = 0.0001;  %   2: air          0.0001
    img.elem_data([fmdl.mat_idx{3}]) = 0.47;    %   3: grey matter  0.47
    img.elem_data([fmdl.mat_idx{4}]) = 0.016;   %   4: skull        0.016

    % Set stim patterns
    [img.fwd_model.stimulation, img.fwd_model.meas_select] = mk_stim_patterns(32,1,[0,5],[0,5],{'no_meas_current_next2'},1);
    img.fwd_model.normalize_measurements = 0;
    imdl = mk_GREIT_model(img, radius, weight, opt);
    save('fmdl_10_2.mat','fmdl');
    save('imdl_10_2.mat','imdl');
    
end % end if
    
end % end function

function [fwd_mdl, new_electrode_centers]= map_electrodes_to_surface(fwd_mdl, electrode_centers)
% Using an array of Nx3 electrode center coordinates located outside of the
% model, find where the center of the electrodes should be on the surface
% of the model
% 
% fwd_model: and EIDORS forward model object. 
% electrode centers: An Nx3 array specifying the coordinates of electrodes
% in space outside of the model. Electrode intersections will be calculated
% assuming a planar ring.
%
% Returns an array of electrode centers located on the surface of the mesh.
% This function is written by :
%                               Mark Campbell
%                               Carleton University
faces= fwd_mdl.boundary;
nodes= fwd_mdl.nodes;
num_elec= size(electrode_centers, 1);
elec_z= electrode_centers(1, 3);
new_electrode_centers= electrode_centers;
vert1= nodes(faces(:, 1), :);
vert2= nodes(faces(:, 2), :);
vert3= nodes(faces(:, 3), :);
elec_node_idx_in_fmdl= zeros(32, 1);
for e= 1:num_elec
    p0= repmat(electrode_centers(e, :), size(vert1, 1), 1);
    p1= [0 0 elec_z];
    a= sum([(vert1- p0).^2, (vert2- p0).^2, (vert3- p0).^2], 2); % find the closest face
    closest_face= find(a==min(a));
    v1= vert1(closest_face, :);
    v2= vert2(closest_face, :);
    v3= vert3(closest_face, :);
    n= cross(v1- v2, v1-v3);
    new_node= plane_line_intersect(n, v1, p0(1,:), p1);
    new_electrode_centers(e, :)= new_node;
    [fwd_mdl, node_idx]= move_closest_to_e_pos(fwd_mdl, v1, v2, v3, faces(closest_face, :), new_node);
    elec_node_idx_in_fmdl(e)= node_idx;
end % end for

num_nodes= length(fwd_mdl.nodes);
for i= 1:num_elec
    idx= elec_node_idx_in_fmdl(i);
    fwd_mdl.electrode(i).nodes= idx;
    fwd_mdl.electrode(i).z_contact= 0.01;
end % end for

end % end function

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

function [fmdl, node_idx]= move_closest_to_e_pos(fwd_model, v1, v2, v3, v_idx, new_node)
    fmdl= fwd_model;
    a= sum([(v1- new_node).^2; (v2- new_node).^2; (v3- new_node).^2], 2); % find the closest face
    node_idx= v_idx(find(a==min(a)));
    fmdl.nodes(node_idx, :)= new_node;
end % end function