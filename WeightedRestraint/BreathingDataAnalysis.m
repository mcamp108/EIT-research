run('myStartup.m');

%% MODEL
numElecs= 32;
numRings= 2;
elecPerRing= numElecs/numRings;
height= 1;
radius= 0.5;
inj= [0,5];
meas= [0,5];
amplitude= 3e-3;

fmdl= ng_mk_cyl_models([height, radius], [elecPerRing, 0.4, 0.6], [0.05, 0 ,0.025] );  %create fmdl

skip4 = {16, 2, inj , meas, {'no_meas_current_next2'}, amplitude}; % copied from Symon's display_3D. Do not make measurements on the two electrodes next to the current carrying
idx = reshape(1:32,numRings,[])'; %matrix of 32 elem numbered from 1:32, shaped with row length= L then transposed so numbering goes across rows
idx(2:2:end,:) = fliplr(idx(2:2:end,:)); % flips every second row to get proper square electrode arrangement
idx= flipud(idx);
fmdl.electrode(idx) = fmdl.electrode(:); % sets (was 16x2, now 32x1)
[fmdl.stimulation, fmdl.meas_select] = mk_stim_patterns(skip4{:}); % 'meas_sel' indicates which electrodes are used
vopt.imgsz = [32 32];
vopt.square_pixels = true;
vopt.zvec = linspace(0.3,0.7,4);
vopt.save_memory = 1;
opt.noise_figure = 1.0;
opt.keep_intermediate_results= true;
[imdl_t, opt.distr] = GREIT3D_distribution(fmdl, vopt);
imdl= mk_GREIT_model(imdl_t, 0.2, [], opt);

% DATA
sref= '2019_08_01_P1_standingReference.eit';
pref= '2019_08_01_P1_proneRef.eit';
wref= '2019_08_01_P1_proneRefWeight.eit';
wepos= '2019_08_01_P1_proneWeightExercisePos3.eit';
epos = '2019_08_01_P1_proneNoWeightExercisePos3.eit';
files= {sref, pref, wref, wepos, epos};

%%
msel= imdl.fwd_model.meas_select;
mm = find(msel);
impedanceFactor =  2.048 / (2^12 * 0.173 * 0.003) / 2^15; % = 0.9633 / 2^15;

for i= files
    [dd,auxdata]= eidors_readdata(i); 
    FR = 1e6./median(diff(auxdata.t_rel)); %framerate is median dif of time points/ 1000000 (convert to s)

    %% CLEAN DATA
    
    elec_impedance= mean(real(abs(auxdata.elec_impedance* impedanceFactor)), 2);
    [imdl_comp, vv_prime]= compensate_bad_elec(dd, elec_impedance, imdl);
    use_data= real(dd(mm, :));

    %% SOLVE
    imgr= inv_solve(imdl_comp, mean(use_data, 2), use_data); % data 1 == referene frame, data 2== other frames in time series.
end % end for
%% ANALYZE
show_slices(imgr, [inf,inf,0.5]);
show_fem(imgr);
figure; plot(sum(real(dd(1:500,:)),1)); % time slice control
xlim([0,100]);
xticks([1:2:length(dd)]);
xticklabels([1:5:length(dd)]/FR);
xtickangle(90);

figure;
clf; subplot(221); plot(vh.meas); ylim(0.8*[-1,1]); % plot vh measurement data
xlim(100+[0,100]); 
subplot(223); plot(real(dd(:,1:16))); % plot real portion
xlim(100+[0,100]);
timeSec= (auxdata.t_rel(end)- auxdata.t_rel(1)) /1e6; % time in seconds
xticks([1:length(dd)]);
xticklabels([1:length(dd)]/FR);
xtickangle(90);

subplot(222); 
plot(dd(:,1:16),'b*'); 
idx = abs(vh.meas) > 0.2; 
hold on; 
plot(dd(idx,1:10),'r*'); 
hold off; 

subplot(224); % plot frequency data
dr = real(dd(idx,:)); 
plot(sum(dr)) 
tax = linspace(0,FR,size(dd,2)); 
semilogy(tax,abs(fft(sum(dr)))) 
xlim([0,FR/2]); 
print_convert SeatedLevelNoEx.png
