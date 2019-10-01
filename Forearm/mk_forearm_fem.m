function [fmdl, mat_idx]= mk_forearm_fem(elec_config, forearm_geo)

% elec_config     - A cell aray containing [number of electrodes, number of
%                   rings]
% forearm_geo     - A cell array containing: [forearm length, distance
%                   between rings, ringRad1, ringRad2, rinRad3, ringRad4]
% Will add 0.5 of an electrode spacing to the top and bottom of the mesh to
% allow reconstruction above and below planes of the electrodes, also to
% allow positioning of electrodes on the mesh.

numElecs= elec_config(1); numRings= elec_config(2);
elecPerRing= numElecs/numRings;
forearmLength= forearm_geo(1); ringDist= forearm_geo(2);
ringRads= linspace(forearm_geo(3), forearm_geo(end), numRings);
ring_zPlanes= linspace(0, forearmLength, numRings)+ringDist/2; 

slope= forearmLength/(ringRads(1)- ringRads(4));
top_extenstion_rad= ringRads(4)- (ringDist/2)/slope;
bot_extenstion_rad= ringRads(1)+ (ringDist/2)/slope;

body_geometry.cone= struct;
body_geometry.cone.bottom_center= [0; 0; 0];
body_geometry.cone.bottom_radius= bot_extenstion_rad;
body_geometry.cone.top_center= [0; 0; forearmLength+ringDist];
body_geometry.cone.top_radius= top_extenstion_rad;
body_geometry.cone.complement_flag= false;

% ** TODO ** Construct forearm by stacking multiple cylinders?

theta = linspace(0, 2*pi, elecPerRing+1); theta(end) = [];
elec= 1;
for r= 1:numRings
    for i = 1:elecPerRing     
%         electrode_geometry{elec}.point = [cos(theta(i))*ringRads(r) sin(theta(i))*ringRads(r) ring_zPlanes(r)];
        electrode_geometry{elec}.sphere.center = [cos(theta(i))*ringRads(r) sin(theta(i))*ringRads(r) ring_zPlanes(r)];
        electrode_geometry{elec}.sphere.radius = 0.1;
        elec= elec+ 1;
    end
end
[fmdl, mat_idx] = ng_mk_geometric_models(body_geometry, electrode_geometry);
end