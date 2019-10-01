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
%This function is written by :
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
    fwd_mdl= add_node_to_mdl(fwd_mdl, faces(closest_face, :), new_node);
end % end for

fwd_mdl= fix_model(fwd_mdl);
num_nodes= length(fwd_mdl.nodes);
for i= 1:num_elec
    idx= num_nodes- (num_elec- i);
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

function fmdl= add_node_to_mdl(fwd_model, closest_face, new_node)
    fmdl= fwd_model;
    fmdl.nodes= [fmdl.nodes; new_node];
    new_node_num= length(fmdl.nodes);
    fmdl.elems= [fmdl.elems; [closest_face, new_node_num]];
end % end function