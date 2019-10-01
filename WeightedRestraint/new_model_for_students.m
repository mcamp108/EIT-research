run('myStartup.m');
% Test the complex 3D mesh to add electrodes

% Generate the fmdl
thorax = shape_library('get','adult_male','boundary');
shape = { 2,
	  {thorax},
	  [4,40],
	  1};

% Set up the electrodes
elec_pos = [16,0,0.75,1.25];
elec_spec = [0.1]; % Radius of circular electrodes

idx= [8,11,12,15,16,19,20,23,24,27, 28,31,32,3,4,7;
      9,10,13,14,17,18,21,22,25,26,29,30,1,2,5,6]';

fmdl = ng_mk_extruded_model(shape, elec_pos, elec_spec);
fmdl.electrode(idx) = fmdl.electrode(:);
ng_opt_file = []; % use the default
[fmdl.stimulation, fmdl.meas_select]=mk_stim_patterns(32,1,[0,5],[0,5],{'no_meas_current_next2'},1); % Skip 4 
fmdl.normalize = 0;

% Show the forward model to confirm
figure(1); clf;
% show_fem(fmdl, [0 1 0]);
axis off

% Load data
cd 'C:\Users\Mark\Documents\GraduateStudies\LAB\WeightedRestraint\data\Mali'; % insert your data directory here
files = filenames();
for f= 1:length(files)
    makeAbstractFigure(files(f, :), fmdl);
end % end for

function makeAbstractFigure(file, fmdl)
title= remove_underscores(file);
[vv, fs, auxdata, stim_x] = data_loader(file);
clip= round(5* fs);

% Filter and clip the data
vv= lowpass(vv', 5, fs)';
vv= vv(:, clip:(size(vv,2)-clip) );

% Visualize raw data
% show_fft(vv, fs);
% figure; plot(sum(vv));
% figure; plot(vv(1:100, :)');

% Find 5 inpiration/expiration pairs for each file
numBreaths= 5;
breaths= findBreaths((sum(vv)));
if length(breaths.insIdx) > length(breaths.expIdx)
    useBreaths= round( linspace(2, length(breaths.expIdx), numBreaths) );
else
    useBreaths= round( linspace(2, length(breaths.insIdx), numBreaths) );
end % end if
end_ex = breaths.expIdx(useBreaths);
end_in = breaths.insIdx(useBreaths);

% Experimental
% measStd= std(vf')';
% keepIdx= abs(measStd- mean(measStd)) < 2* mean(measStd);
% vf(~keepIdx,:)= repmat( mean(vf(keepIdx,:),1),  length(find(~keepIdx)), 1);
% figure; plot(vf(1:100, :)');

time = linspace(0,length(vv)/fs,length(vv));
figure('position',[0,0,1100,350]);
axes('Position', [.07,0.6,0.9,0.3]);
plot(time,sum(vv),'Color',[0.63,0.63,0.63])
hold on;
from= min(sum(vv));
to= max(sum(vv));
for i=1:5
    plot([end_ex(i)/fs, end_ex(i)/fs],[from, to],'Color',[0 0.4470 0.7410],'LineWidth',2);
    plot([end_in(i)/fs, end_in(i)/fs],[from, to],'Color',[0 0.4470 0.7410],'LineWidth',2);
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
vopt.zvec = linspace(0.75,1.25,6); % Target reconstruction planes
vopt.save_memory = 1;
opt.noise_figure = 2.0; % The higher this is the more likely it converges to the correct noise figure
[imdl, opt.distr] = GREIT3D_distribution(fmdl, vopt);
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
    view(225,45)
    axis off
    %annotation('textbox',[ishift-0.07,0.05,0.35,0.35],'String',labels(i),'EdgeColor','none')
end
cd 'C:\Users\Mark\Documents\GraduateStudies\LAB\WeightedRestraint\Figures';
print_convert(char(title+ "ventilation_images.png"));

end % end function
