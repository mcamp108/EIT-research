function seg = get_brain_segmentation_for_pig(pig)

% -------------------------------------------------------------------------
% DESCRIPTION:
%   seg = get_brain_segmentation_for_pig(seq)
% -------------------------------------------------------------------------
% PARAMETERS:
%   pig:        
%       Name of a pig.
% -------------------------------------------------------------------------   
% RETURNS:
%   seg: 
%       A binary matrix where 0 is not brain, and what is brain is numbered
%       1-6, where each number is a 60 degree zone starting from bottom
%       going counterclockwise from the viewer's perspective. Area 1 is
%       pig's bottom left, area 2 is middle left, and so on.
% -------------------------------------------------------------------------   
% AUTHOR:
%   Mark Campbell
%   Carleton University
%   markacampbell@cmail.carleton.ca
% -------------------------------------------------------------------------

switch pig
    case '8-2'
        elec_z_plane= 106; % average of reference dot z planes from elec_loc_refs
    case '9-2'
        elec_z_plane= 106; % average of reference dot z planes from elec_loc_refs
    case '10-2'
        elec_z_plane= 97; % average of reference dot z planes from elec_loc_refs
    case '11-2'
        elec_z_plane= 96; % average of reference dot z planes from elec_loc_refs
    case '12-2'
        elec_z_plane= 116; % average of reference dot z planes from elec_loc_refs
end % end switch

cd(horzcat('E:\University\Masters\HamburgBrain\Models\', pig,'\mesh'));

seg3dOutMatFile = horzcat(pig, '_seg3D_out');
load(seg3dOutMatFile);
V = scirunnrrd.data;
V = rot90(V, 1); % looking through front of pig

elecPlaneImg = squeeze( V(:, :, elec_z_plane) );
rowsum = sum(elecPlaneImg, 2);
colsum = sum(elecPlaneImg, 1);
elecPlaneImg = elecPlaneImg(rowsum>0, colsum>0);
elecPlaneImg(elecPlaneImg~=3) = 0;
seg = imresize(elecPlaneImg, [64 64]);
seg(seg>0) = 1;

% find brain centroid
rows = find(sum(seg, 2)>0);
cols = find(sum(seg, 1)>0);
cX = mean(cols);
cY = mean(rows);
r=40;
cnt = 0;

for i = -90:60:210
    cnt = cnt + 1;
    t = linspace(-i,-i-60,128);
    x = [cX, cX+r*cosd(t), cX];
    y = [cY, cY+r*sind(t), cY];
    bw = poly2mask( x, y, size(seg,1),size(seg,2));
    seg(bw) = seg(bw)*cnt;
end

end % end function
