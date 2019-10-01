cd 'C:\EIDORS\eidors';
run 'startup.m';
%%
% body_geometry.cone= struct;
% body_geometry.cone(1).top_center = [0 0 2];
% body_geometry.cone(1).bottom_center = [0 0 0];
% body_geometry.cone.bottom_radius= 2;
% body_geometry.cone.top_radius= 1;
% electrode_geometry= {};
% fmdl = ng_mk_geometric_models(body_geometry, electrode_geometry);
% rmdl = mk_grid_model([],-2:2,-2:2,0:2);
% c2f  = mk_grid_c2f(fmdl,rmdl);
% h = show_fem(fmdl); set(h,'EdgeColor','r','LineWidth',0.1)
% hold on
% h = show_fem(rmdl); set(h,'EdgeColor','b','LineWidth',2);
% hold off
%%
% % 3D cylinder with radius 1. One plane of 16 electrodes with radius 0.1
% body_geometry.cylinder = struct;
% n_elect = 16;
% theta = linspace(0, 2*pi, n_elect+1); theta(end) = [];
% for i = 1:n_elect
%     electrode_geometry{i}.sphere.center = [cos(theta(i)) sin(theta(i)) 0.5];
%     electrode_geometry{i}.sphere.radius = 0.1;
%  end
% fmdl = ng_mk_geometric_models(body_geometry, electrode_geometry);
%%
% 3D cone with radius 1. One plane of 16 electrodes with radius 0.1
body_geometry.cone= struct;
body_geometry.cone(1).top_center = [0 0 1];
body_geometry.cone(1).bottom_center = [0 0 0];
n_elect = 16;
theta = linspace(0, 2*pi, n_elect+1); theta(end) = [];
for i = 1:n_elect
    electrode_geometry{i}.sphere.center= [cos(theta(i))*0.75 sin(theta(i))*0.75 0.5];
%     elec_pos= [elec_pos; electrode_geometry{i}.sphere.center];
    electrode_geometry{i}.sphere.radius = 0.01;
end
fmdl = ng_mk_geometric_models(body_geometry, electrode_geometry);
clf;
show_fem(fmdl);
img= mk_image(fmdl, 1);

show_fem(fmdl);  
[fmdl.stimulation,fmdl.meas_select] = mk_stim_patterns(16,1,[0,5],[0,5],{},1);

vopt.imgsz = [32 32];
vopt.zvec = linspace( 0,1,8);
vopt.save_memory = 1;
opt.imgsz = [32 32];
opt.square_pixels = true;
opt.noise_figure = 0.5;
opt.target_size= 0.2;

[imdl_t,opt.distr] = GREIT3D_distribution(fmdl, vopt); % Makes a tet model
imdl_t.interp_mesh.n_points = 1;
% inv_mdl= select_imdl(imdl, {'TV solve dif'});
imdl= mk_GREIT_model(imdl_t, 0.20, [], opt);


% % [cmdl,coarse2fine]= mk_grid_model(fmdl, linspace(-1,1,32), linspace(-1,1,32), linspace(0,1,16));
% 
% [OUT FMDL] = mk_voxel_volume(fmdl, vopt);
% OUT.electrode= fmdl.electrode;
% OUT.stimulation= fmdl.stimulation;
% OUT.meas_select= fmdl.meas_select;
% % OUT.electrode= fmdl.electrode;


