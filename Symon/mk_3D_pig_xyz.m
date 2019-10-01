function [] = mk_3D_pig_xyz()

file_name= 'UCL06_pighead.mat';
external_name= 'UCL06_pighead_external.XYZ';

% Using many equations for lung, esophagus and heart boundaries 
% Generated from CT segmentation using a watershed algorithm 
% Create a 3D EIT inverse and forward model.
addpath('C:\Users\Mark\Documents\GraduateStudies\LAB\HamburgBrain\Data');

% STEP 1 - get the segmented boundaries for each row
% adapted from image processing course:
% This could take a while - save these as a separate file once done
if exist(file_name, 'file') == 2
    bounds = load(file_name);
    bounds = bounds.bounds;
    bounds.length= length(bounds.exterior);
else
    bounds = ct_segmentation_pig_head(); % Only run if we have to... Takes ages!!
    cd 'C:\Users\Mark\Documents\GraduateStudies\LAB\HamburgBrain\Data';
    save(file_name,'bounds'); % Save it to save time
end % end if

xyres= 0.1895;
zres= 0.45;
xyscaling_factor= zres/ xyres;
% For each layer downsample the plots and shift to be centered on 0
ext = [];
n = 1; % Keep every nth element only
% Scale and offset - the largest dimension should be between 0 and 1 so it fits within a 1 unit cube
offset = 0; % Centre on 0
bounds.exterior= center_model(bounds.exterior);

for i=1:bounds.length % For each layer
    temp = bounds.exterior{i};
    if ~isempty(temp)
        temp = temp(1:n:end,:);
        temp = (temp-offset)./ (xyscaling_factor* bounds.length);
        temp(:,3) = i/ bounds.length;
        ext = [ext; temp];
    end % end if
end % end for
% Do a nice nice plot
figure(1);
clf;
scatter3(ext(:,1),ext(:,2),ext(:,3),'.');
% STEP 2 - Save pointcloud to .XYZ files
cd 'C:\Users\Mark Campbell\Documents\GraduateStudies\LAB\HamburgBrain\Data';
dlmwrite(external_name, ext, 'delimiter', ' ', 'precision', 5);
end % end function

function boundsexterior= center_model(boundsexterior)

xyz_minmax= find_xyz_minmax(boundsexterior);
x_center= round( (xyz_minmax.x.max- xyz_minmax.x.min)/ 2);
y_center= round( (xyz_minmax.y.max- xyz_minmax.y.min)/ 2);

for layer= 1:length(boundsexterior)
    if isempty(boundsexterior{layer})
        continue
    else
        boundsexterior{layer}(:, 1)= (boundsexterior{layer}(:, 1)- x_center);
        boundsexterior{layer}(:, 2)= (boundsexterior{layer}(:, 2)- y_center);
    end % end if
end % end for

end % end function


function xyz_minmax= find_xyz_minmax(boundsexterior)

xyz_minmax= struct;
got_ref= 0;
for layer= 1:length(boundsexterior)
    if isempty(boundsexterior{layer})
        continue
    end % end if
    if got_ref== 0
        got_ref= 1;
        xyz_minmax.x.min= min(boundsexterior{layer}(:, 1));
        xyz_minmax.x.max= max(boundsexterior{layer}(:, 1));
        xyz_minmax.y.min= min(boundsexterior{layer}(:, 2));
        xyz_minmax.y.max= max(boundsexterior{layer}(:, 2));
    end % end if
    if min(boundsexterior{layer}(:, 1)) < xyz_minmax.x.min
        xyz_minmax.x.min= min(boundsexterior{layer}(:, 1));
    end % end if
    if max(boundsexterior{layer}(:, 1)) > xyz_minmax.x.max
        xyz_minmax.x.max= max(boundsexterior{layer}(:, 1));
    end % end if
    if min(boundsexterior{layer}(:, 2)) < xyz_minmax.y.min
        xyz_minmax.y.min= min(boundsexterior{layer}(:, 2));
    end % end if
    if max(boundsexterior{layer}(:, 2)) > xyz_minmax.y.max
        xyz_minmax.y.max= max(boundsexterior{layer}(:, 2));
    end % end if
end % end for

end % end function