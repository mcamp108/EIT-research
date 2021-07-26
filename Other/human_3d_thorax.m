function out = human_3d_thorax(in)
% -------------------------------------------------------------------------
% DESCRIPTION:
%
% -------------------------------------------------------------------------
% PARAMETERS:
% 
% -------------------------------------------------------------------------   
% RETURNS:
% 
% -------------------------------------------------------------------------   
% AUTHOR:
%   Mark Campbell
%   Carleton University
%   markacampbell@sce.carleton.ca
% -------------------------------------------------------------------------
% VERSION:
%   1.0.0
% -------------------------------------------------------------------------
% (C) 2019-2021 Mark Campbell. License: GPL version 2 or version 3
% -------------------------------------------------------------------------

fmdl= mk_library_model('adult_male_16el_lungs');
img = mk_image(fmdl, 1); % background conductivity
img.elem_data(fmdl.mat_idx{2}) = 0.3001; % lungs
img.elem_data(fmdl.mat_idx{3}) = 0.3002; % lungs
ROI = calc_slices(img,[inf,inf,0.5]);
llung_ROI = ~isnan(ROI) & (ROI==0.3001);
rlung_ROI = ~isnan(ROI) & (ROI==0.3002);
thorax_ROI= ~isnan(ROI); % include lungs, too

subplot(131); imagesc(thorax_ROI); axis image
subplot(132); imagesc(rlung_ROI);  axis image
subplot(133); imagesc(llung_ROI);  axis image
print_convert('GREIT_IBEX_01a.png');

[stim,msel] = mk_stim_patterns(16,1,[0,1],[0,1],{'no_meas_current'},1);
img.fwd_model.stimulation = stim;
img.fwd_model = mdl_normalize(img.fwd_model, 1);
opt.imgsz = [64 64];
opt.distr = 3;
opt.Nsim = 500;
opt.target_size = 0.03;
opt.target_offset = 0;
opt.noise_figure = .5; 
opt.square_pixels = 1;
imdl=mk_GREIT_model(img, 0.25, [], opt);
imdl.fwd_model.meas_select = msel;
end % end function


fname = 'DATA/STUDYNAME/SUBJECT_1/YYYYMMDD/Eit/Viasys/1001_c4.get';
vv= eidors_readdata(fname);
img= inv_solve(imdl, mean(vv,2), vv);

imgs= -calc_slices(img); % Negative to air is +
imgs(isnan(imgs(:)))= 0;

img.calc_colours.ref_level=0;
img.elem_data = img.elem_data(:,2:4:120);
img.show_slices.img_cols = 10;
clf; show_slices(img);

print_convert 'GREIT_IBEX_03a.jpg'


data.imageRate = 13;

data.patient.ROI.Inside =thorax_ROI*100; % to scale it up to 100
data.patient.ROI.RightLung =rlung_ROI*100;
data.patient.ROI.LeftLung =llung_ROI*100;
data.patient.ROI.Heart =zeros(size(imgs,1),size(imgs,2));

% put to dummy because they are missing
data.patient.halfChest = 'NaN';
data.patient.height = 'NaN';
data.patient.weight = 'NaN';
data.patient.gender = 'NaN';

data.measurement.Position.transversal = zeros (1, size(imgs, 3));
data.measurement.Position.longitudinal = zeros (1, size(imgs, 3));
data.measurement.ImageQuality = 100 * ones(1, size(imgs, 3));
data.measurement.ElectrodeQuality = zeros(size(imgs, 3), 32);
data.measurement.ZeroRef = imgs;

data.injctionPattern= 'NaN';
data.SensorBelt.NumEl= 'NaN';

data.measurement.CompositValue=squeeze(sum(sum(imgs,2),1));

save('file-for-IBEX.mat','data');