cd 'C:\EIDORS\eidors';
run 'startup.m';
cd 'C:\Users\Mark Campbell\Documents\GraduateStudies\LAB\Pulse Wave Velocity\Data';
%%
numElecs= 32;
numRings= 4;
avRingRad= 4.02;
elecPerRing= numElecs/numRings;
forearmLength= 15; % height
ringDist= 5;
skip4 = {32,1,[0,5],[0,5],{'meas_current'},1}; % now the stim pattern. currents injected from 0-5, measured from 0-5 (skip 4) - meas_current enables measurements from current-carrying electrodes
idx = reshape(1:32,numRings,[])'; %matrix of 32 elem numbered from 1:32, shaped with row length= L then transposed so numbering goes across rows
idx(2:2:end,:) = fliplr(idx(2:2:end,:)); % flips every second row to get proper square electrode arrangement
idx(2:end, :)= flipud(idx(2:end, :));

% ** Forward model **
fmdl= ng_mk_cyl_models([forearmLength+ ringDist, avRingRad, .5], [elecPerRing, linspace(0, forearmLength ,numRings)], [0.05]);
% fmdl= mk_forearm_fem([numElecs, numRings], [forearmLength, ringDist, 4.77, 4.54, 3.66, 3.1]);
fmdl.electrode(idx) = fmdl.electrode(:); % sets idx to fmdl.electrode in column order (models ordering of electrodes up and down the arm)
[fmdl.stimulation, fmdl.meas_select] = mk_stim_patterns(skip4{:}); 
img= mk_image(fmdl,1);
vh= fwd_solve(img);
show_fem(img, [0 1 0]);

% ** Inverse model **
vopt.imgsz = [32 32];
vopt.square_pixels = true;
vopt.zvec= linspace(0,forearmLength+ringDist,7); % try 6 next time. targetsize= 0.17 for 10 0.20 for 6
vopt.save_memory = 1;
[imdl_t,opt.distr] = GREIT3D_distribution(fmdl, vopt);

opt.noise_figure = 4;
% opt.target_size= 0.17; % this is the lowest i can get it for the target NF.
imdl= mk_GREIT_model(imdl_t, 0.2, [], opt);
imdl.fwd_model= fmdl;
%%
files= {'Seated_ArmLevelWithECG.eit'}; %'Seated_ArmLevelStraight.eit'};
% files= {'StandingArmLoweredNoEx.eit', 'StandingBentArmLoweredEx-18PU.eit', 'StandingArmRaisedEx-17PU.eit', 'StandingArmLoweredEx-30PU.eit', 'StandingBentOverArmLoweredNoEx.eit', 'StandingArmRaisedNoEx.eit', 'SeatedRaisedNoEx.eit', 'SeatedLevelNoEx.eit', 'SeatedLowered.eit'};

for i= 1: length(files)
    F1 = files{i};
    [dd,auxdata]= eidors_readdata(F1); 
    FR = 1e6./median(diff(auxdata.t_rel)); %framerate is median dif of time points/ 1000000 (convert to s)
    filename= F1(1:end-4);
    % ** FILTER WITH HEART RATE FREQUENCY **
    filtdd= dd;
    for meas= 1:1024
        dr= real(dd(meas,:));
        filtdd(meas, :)= bandpassfilt(6, dr, FR, 0.8, 1.3);
    end
    dd= filtdd;
 
    % ** INVERSE SOLVE ** (should throw away first and last 200 frames)
    firstFrame= 640;
    lastFrame= floor(length(dd));
    data1= mean(dd, 2);
    data2= real(dd(:,firstFrame:lastFrame));
    reconst= inv_solve(imdl, data1, data2);
    
    % ** DETREND **
%     reconst= detrend_3dmdl(reconst);
    
    % ** SET COLOUR LIMITS **
    clim= max(reconst.elem_data(:));
    cmin= min(reconst.elem_data(:));
    reconst.calc_colours.ref_level= 0;
    reconst.calc_colours.lim= clim;
    
    % ** FIND BLOOD FLOW PULSE PEAKS **
    ECGpmtEIT= csvread('ECGpmtEIT.csv');
    ECGdata= csvread('ECGdatamtEIT.csv');
        
    % ** ENSEMBLE AVERAGING **
    ensembl_opt.num= 1;
    ensembl_opt.min_frames= 38;
    ensembl_opt.offset= -5;
    reconst_temp = reconst;
    [reconst_temp.elem_data, peaks_used]= ensmblav(reconst_temp.elem_data, ECGpmtEIT, ensembl_opt);
    [ECG_ensemble, ecgpeaks_used]= ensmblav(ECGdata', ECGpmtEIT, ensembl_opt); figure; plot(ECG_ensemble);
    
    zpln= linspace(1,19,6);
    show_slices(reconst_temp, [inf, inf, zpln(1); inf, inf, zpln(2); inf, inf, zpln(3); inf, inf, zpln(4); inf, inf, zpln(5); inf, inf, zpln(6); inf, inf, zpln(7); inf, inf, zpln(8); inf, inf, zpln(9)],39);
    keyboard;
%%
reconst_temp2= reconst_temp;
angles= [270];
angle2= -45;
for angle1= angles
    % ** VISUALIZE BLOODFLOW **
%     vidfile = VideoWriter(filename+ "withECGref_Pulse_1Beat_3D_6ring_meanRef_"+ angle1+ angle2+ ".mp4",'MPEG-4');
    vidfile = VideoWriter(filename+ "scrap"+ angle1+ angle2+ ".mp4",'MPEG-4');
    vidfile.FrameRate= FR/4;
    reconst_temp2= reconst_temp;
    open(vidfile);
    figure;
    for frame= 1:size(reconst_temp.elem_data, 2)
        reconst_temp2.elem_data = reconst_temp.elem_data(:,frame);
        
        show_3d_slices(reconst_temp2, [0:20], [-3:-1], [-5:5]); view(angle1, angle2);
%         subplot(1, 2, 2); plot(ECG_ensemble); hold on; plot(ECG_ensemble(1:frame), 'red'); hold off;
%                 show_fem(reconst_temp2); view(angle1, angle2);
        colorbar; caxis([cmin clim]);
        writeVideo(vidfile,getframe(gcf));
    end
    close(vidfile);
end
%% Figure for paper
figure;

ha= tight_subplot(4,10,[0.01, 0.01], [0.01, 0.01], [0.01, 0.01]);
for frame= 1:size(reconst_temp.elem_data, 2)
    reconst_temp2.elem_data = reconst_temp.elem_data(:,frame);
    axes(ha(frame));
    show_3d_slices(reconst_temp2, [1:20], [-5:5], [-4:4]); view(250, -30);
end
colorbar; caxis([cmin clim]);
%%
    keyboard;
    
    % ** Visualize reconstruction **
    reconst_temp = reconst;
    filename= F1(1: end-4);    
    figureHandle = figure;
    vidfile = VideoWriter(filename+ "FEM_NF1_PointElec_detrended3D.mp4",'MPEG-4');
    vidfile.FrameRate= FR;
    open(vidfile);
    for frame= 1: size(reconst_slices1, 3)-1
        subplot(5,2,2); show_slices(reconst_slices5(:,:,frame)); title("15 cm (Distal)");
        subplot(5,2,4); show_slices(reconst_slices4(:,:,frame)); title("12.5 cm");
        subplot(5,2,6); show_slices(reconst_slices3(:,:,frame)); title("10 cm");
        subplot(5,2,8); show_slices(reconst_slices2(:,:,frame)); title("7.5 cm");
        subplot(5,2,10); show_slices(reconst_slices1(:,:,frame)); title("5 cm (Proximal)");
        axesHandles = findobj(get(figureHandle,'Children'), 'flat','Type','axes');
        axis(axesHandles,'square');
        reconst_temp.elem_data = reconst.elem_data(:,frame);
        subplot(5,2,[1,3,5]); show_fem(reconst_temp);
        colorbar; caxis([cmin clim]);
        writeVideo(vidfile,getframe(gcf));
    end    
    close(vidfile);
end

%%
figure; show_slices(reconst_slices1); figure; show_slices(reconst_slices2); figure; show_slices(reconst_slices3); figure; show_slices(reconst_slices4); figure; show_slices(reconst_slices5);
% Plotting data

% PLOT DATA FROM ARTERY
figure; show_slices(reconst_slices1); figure; show_slices(reconst_slices2); figure; show_slices(reconst_slices3); figure; show_slices(reconst_slices4); figure; show_slices(reconst_slices5);
% figure; show_slices(detrend_slices1); figure; show_slices(detrend_slices2); figure; show_slices(detrend_slices3); figure; show_slices(detrend_slices4); figure; show_slices(detrend_slices5);
watch_artery= [];
for row= 27
    for column= 25
        add= reconst_slices1(row, column, :);
        watch_artery= [watch_artery; add(:)'];
    end
end

hold on;
for i= 1:size(watch_artery, 1)
    plot(watch_artery(i,:));
end
hold off;
[pks, locs] = findpeaks(watch_artery);

% pulsemean= mean(watch, 1);
% detrend_pulsemean = detrend(pulsemean);
% trend = pulsemean - detrend_pulsemean;
% 
% plot(pulsemean);
% hold on;
% plot(detrend_pulsemean);
% hold off;


% PLOT CONVOLUTION
vec1= sum(data2(:, 1:floor(FR*2)), 1);
vec2= sum(data2, 1);
w = conv(vec2, vec1);
plot(w);

first= w(86:100);
second= w(86:3600);
figure;
v = conv(first, second);
plot(v);

figure;
autocorr(sum(real(dd(:,:)),1), 'NumLags', 1600);

% plot frequency spectrum of measurements
tax = linspace(0,FR,size(dd,2));
figure;
hold on;
for meas= 1:1024
    dr= real(dd(meas,:));
    semilogy(tax,abs(fft(dr)));
end
xlim([0,FR/2]);
hold off;

figure; plot(sum(real(dd(:,1:500)),1));
xlim([0,100]);
xticks([1:5:length(dd)]);
xticklabels([1:5:length(dd)]/FR);
xtickangle(90);
% % 
figure;
clf; subplot(221); plot(vh.meas); ylim(0.8*[-1,1]); % plot vh measurement data
xlim(100+[0,100]); 
subplot(223); plot(real(dd(:,150:800))); % plot real portion
xlim(100+[0,100]);
% 
subplot(222); 
plot(dd(:,:),'b*'); 
idx = abs(vh.meas) > 0.1; 
hold on; 
plot(dd(idx,:),'r*'); 
hold off; 

subplot(224); % plot frequency data
dr = real(dd(idx,:));
plot(sum(dr));
tax = linspace(0,FR,size(dd,2)); 
semilogy(tax,abs(fft(sum(dr)))); 
xlim([0,FR/2]);

dr= sum(real(dd(:,:)),1);
tax = linspace(0,FR,size(dd,2)); 
semilogy(tax,abs(fft(dr))); 
xlim([0,FR/2]);


function ensemble= get_ensemble_average(data, window_size)
    ensemble= data(:, 1: window_size);
    for i= 2:floor( length(data)/ window_size )
        start= (i-1)*window_size +1;
        stop= i*window_size;
        ensemble= ensemble+ data(:, start: stop);
    end
    ensemble= ensemble/i;
    
end

function detrended_slices= detrend_slices(slices, options)
% DETREND SLICE DATA

% 1. Find the mean value of each pixel over the time series. 
% 2. Detrend using a copy of the data then subtract original from detrended to get the trend of that pixel over time. 
% 3. Average all trends to get average slope
% 4. Add mean of original measurement to measurement and subtract (timepoint* slope) from each measurement.
% If asGroup is false, then each pixel is detrended individually and the mean of each pixel is preserved
asGroup= false;
if nargin >1
    if ~isfield(options, 'as_group')
        asGroup = true;
    end
end
detrended_slices= slices;
sz_data= size(slices);
rows= sz_data(1);
cols= sz_data(2);
ts= sz_data(3);
db= struct;
db.idx= [];
for i= 1:rows
    for j= 1:cols
        if ~isnan(slices(i, j, :))
            db.idx= [db.idx; i j];
        end
    end
end

valid= size(db.idx, 1);
db.means= zeros(valid, 1);
db.trends= zeros(valid, ts);
db.dts= zeros(valid, ts);
for i =1:valid
    pixel_idx= db.idx(i, :);
    temp= slices(pixel_idx(1), pixel_idx(2), :);
    temp= temp(:);
    if ~asGroup
        detrended_slices(pixel_idx(1), pixel_idx(2), :)= (detrend(temp) + mean(temp))';
    else
        db.dts(i, :)= temp;
        db.means(i)= mean(temp);
        trend= temp- detrend(temp);
        db.trends(i, :)= trend;
    end
end

if asGroup % detrend based on mean of signal trends
    mean_trend= mean(db.trends);
    for i =1:valid
        pixel_idx= db.idx(i, :); row= pixel_idx(1); col= pixel_idx(2);
        temp= slices(row, col, :);
        detrended_slices(row, col, :)= (temp(:)'- mean_trend) + db.means(i);
    end
end

end

function plot_sidebyside( list_of_slices, frames)
% ** TODO MAKE THIS MORE LIKE SHOW_SLICES
amalgum= [];
for i= 1:length(list_of_slices)
    temp= list_of_slices{i};
    amalgum= cat(3, amalgum, temp(:,:,frames));
end
show_slices(amalgum, [], length(frames));
end