function elecPos= getElecPos(mdl1, elTheta, elPhi)
% Find electrode xyz coordinates from anatomical landmarks (must be fields
% in the model)

% find the polar coordinates of origin and radius of imaginary sphere from 
% nasion and preauriculars landmarks
alph= 65; % 25 degree angle between nasion-inion line and level plane
Nz= mdl1.Nz;
Iz= mdl1.Iz;
PAL= mdl1.PAL;
PAR= mdl1.PAR;
ap= (Iz-Nz)/norm(Iz-Nz); % anterior/posterior unit vector
apM= (Nz+ Iz)/ 2; % midpoint between nasion and inion
apL= norm(Nz- Iz); % distance between nasion and inion
lr= (PAL-PAR)/norm(PAL-PAR); % left/right unit vector
si= cross(lr, ap); % superior/ inferior unit vector
Or= apM + (apL* si)/ (2*tand(alph)); % location of origin of imaginary sphere
rad= norm(Or- Nz);
T= [lr; ap; si; Or]; % compose transformation matrix

% find Cartesian coordinates and do transformation to RAS coordinates
elX= rad* sind(elTheta).* cosd(elPhi);
elY= rad* sind(elTheta).* sind(elPhi);
elZ= rad* cosd(elTheta);
elX1= elX.* T(1,1)+ elY* T(1,2)+ elZ* T(1,3)+ T(4,1);
elY1= elX.* T(2,1)+ elY* T(2,2)+ elZ* T(2,3)+ T(4,2);
elZ1= elX.* T(3,1)+ elY* T(3,2)+ elZ* T(3,3)+ T(4,3);
elecPosABC= [elX1; elY1; elZ1]';
mdl = prepare_surf_model(mdl1);

% Find electrode xyz coordinates
V1= mdl.nodes(mdl.boundary(:, 1), :);
V2= mdl.nodes(mdl.boundary(:, 2), :);
V3= mdl.nodes(mdl.boundary(:, 3), :);
elecPos= zeros(32, 3);
for e= 1:32
    elV= elecPosABC(e,:);
    elecPos(e,:)= getEPos(elV, Or, mdl, V1, V2, V3);
end % end for
% figure; show_fem(mdl); hold on
% % plot3([Or(1) elecPosABC(1,1)], [Or(2) elecPosABC(1,2)], [Or(3) elecPosABC(1,3)]);
% % plot3(Vs(:,1), Vs(:,2), Vs(:,3));
% plot3(elecPos(:,1), elecPos(:,2), elecPos(:,3), 'o'); hold off
% keyboard;
end % end function

% Extract a nice surface model from the one given
function mdl = prepare_surf_model(mdl)

try mdl = rmfield(mdl,'boundary');  end
mdl = fix_model(mdl);
mdl = orient_boundary(mdl);
mdl.elems = mdl.boundary;
mdl.faces = mdl.boundary;
mdl.face_centre = mdl.face_centre(mdl.boundary_face,:);
mdl.normals = mdl.normals(mdl.boundary_face,:);
mdl = rmfield(mdl, 'inner_normal');
mdl = rmfield(mdl, 'boundary_face');
idx = nchoosek(1:3, 2);
elem_sorted = sort(mdl.elems,2);
[mdl.edges ib ia] = unique(reshape(elem_sorted(:,idx),[],2),'rows');
D = mdl.nodes(mdl.edges(:,1),:) - mdl.nodes(mdl.edges(:,2),:);
mdl.edge_len = sqrt(sum(D.^2,2));

end % end function

function mdl = orient_boundary(mdl)

% consistently orient boundary elements
flip = mdl.elem2face(logical(mdl.boundary_face(mdl.elem2face).*mdl.inner_normal));
mdl.faces(flip,:) = mdl.faces(flip,[1 3 2]);
mdl.normals(flip,:) = -mdl.normals(flip,:);
mdl.boundary = mdl.faces(mdl.boundary_face,:);

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

function ePos= getEPos(p, Or, mdl, V1, V2, V3, depth)
if nargin< 7
    depth= 0;
end % end if
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
    depth= depth+ 1;
    if depth< 10
        ePos= getEPos(ePos, Or, mdl, V1, V2, V3, depth);
    else
        disp("WARNING: Recursion depth limit reached!");
        keyboard;
    end % end if
end % end if

end % end function

function fc= findClosestFace(p, mdl, V1, V2, V3)

p0= repmat(p, size(mdl.boundary, 1), 1);
dif= sum([(V1- p0).^2, (V2- p0).^2, (V3- p0).^2], 2); % find the closest face
fc= find(dif==min(dif));

end % end function