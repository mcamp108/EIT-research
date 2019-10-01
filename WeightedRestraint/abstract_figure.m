run('myStartup.m');
% eidors_startup();
% Test the complex 3D mesh to add electrodes

% Load the volume file
ng_vol_filename = 'C:\Users\MarkCampbell\Documents\GraduateStudies\LAB\RestraintBreathing\Models\external2.vol';
% set parameters required to run forward model maker - will be overwritten later
cent = [1,1,1];
stim = [0];

% Generate the fmdl
fmdl = ng_mk_fwd_model(ng_vol_filename, cent, 'ng',stim,50);
sep = 360/16;
row_1 = zeros(32,2);
row_2 = row_1;

% BLAIR & MARK ELECTRODE CONFIG
%row_1([1:16]*2-1,1) = angles'; 
%row_1([1:16]*2-1,2) = 9;
%row_2([1:16]*2,1)   = angles'; 
%row_2([1:16]*2,2) = 6;
%elec_pos = row_1+row_2; 
% Square electrode configuration 
%elec_pos = reshape(1:32,2,[])';
%idx(2:2:end,:) = fliplr(idx(2:2:end,:)); % TODO check this fits with the setup used

% SYMON ELECTRODE CONFIG
angles = [1:16]*sep-190-22.5; %-sep/2;
% angles = [1:16]*sep-190-22.5; %-sep/2;
row_1([3,4,7,8,11,12,15,16,19,20,23,24,27,28,31,32],1) = fliplr(angles)'; 
row_1([3,4,7,8,11,12,15,16,19,20,23,24,27,28,31,32],2) = 9;
next_angles = angles;
next_angles(1) = [];
next_angles(16) = angles(1);
row_2([1,2,5,6,9,10,13,14,17,18,21,22,25,26,29,30],1) = fliplr(next_angles)'; 
row_2([1,2,5,6,9,10,13,14,17,18,21,22,25,26,29,30],2) = 6;
elec_pos = row_1+row_2; 
elec_spec = [0.5, 0, 0.1]; % Radius of circular electrodes
ng_opt_file = []; % use the default
maxh = 0.5; % max spacing on the remesh maybe? 0.5 and larger does not work - nohelpful error message
fmdl= fix_model(fmdl);
fmdl2 = place_elec_on_surf(fmdl, elec_pos, elec_spec, ng_opt_file, maxh);
[stim,msel]=mk_stim_patterns(32,1,[0,5],[0,5],{'no_meas_current_next2'},1); % Skip 4 
fmdl2.stimulation = stim;
fmdl2.meas_select = msel;
fmdl2.normalize = 0;

% Show the forward model to confirm
%figure(1); clf;
%show_fem(fmdl2,2)
%axis off

% Load the data
cd 'C:\Users\Mark\Documents\GraduateStudies\LAB\RestraintBreathing\data\subj1';
files = filenames();
file = files(1,:); % Manually select a file - Prone 50lbs no exercise
[vv, fs, auxdata, stim_x] = data_loader(file);
time = linspace(0,length(vv)/fs,length(vv));

figure('position',[0,0,1100,350])
axes('Position', [.07,0.6,0.9,0.3])
% Get the points
%plot(sum(vv))
%ginput(10)
end_ex = [169,1705,2953,4177,5455];
end_in = [267,1798,3015,4257,5535];
plot(time,sum(vv),'Color',[0.63,0.63,0.63])
hold on;
for i=1:5
    plot([end_ex(i)/fs, end_ex(i)/fs],[0.116,0.128],'Color',[0 0.4470 0.7410],'LineWidth',2);
    plot([end_in(i)/fs, end_in(i)/fs],[0.116,0.128],'Color',[0 0.4470 0.7410],'LineWidth',2);
end
h_xlabel = xlabel('time (seconds)');
h_ylabel = ylabel('\Delta z');

set(gca,'FontName','Helvetica');
set(gca, ...
  'Box'       ,'off'    , ...
  'TickDir'   ,'out'    , ...
  'TickLength',[.007 .007], ...
  'LineWidth' , 1 );
set([h_xlabel, h_ylabel], 'FontName' , 'AvantGarde', 'FontSize', 12);
hold off

% Make the inverse model 
vopt.imgsz = [32 32]; % Image size?
vopt.square_pixels = true;
vopt.zvec = linspace(6,8,6); % Target reconstruction planes
vopt.save_memory = 1;
opt.noise_figure = 2.0; % The higher this is the more likely it converges to the correct noise figure
[imdl,opt.distr] = GREIT3D_distribution(fmdl2, vopt);
imdl3= mk_GREIT_model(imdl, 0.20, [], opt);

% Now the images...
labels = ['A','B','C','D','E'];
for i=1:5
    % Select a small segment
    vd = vv(:,end_in(i));
    vh = vv(:,end_ex(i));
    % Reconstruct some data
    img = inv_solve(imdl3, vh, vd);
    ishift = 0.18*(i-1);
    axes('Position', [ishift-0.07,0.05,0.42,0.42])
    img.calc_slices.clim = 5;
    show_3d_slices(img)
    axis off
    %annotation('textbox',[ishift-0.07,0.05,0.35,0.35],'String',labels(i),'EdgeColor','none')
end
print('ventilation_images.png')



