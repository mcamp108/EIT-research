% This function was written by :
%                               Mark Campbell
%                               Carleton University
run 'myStartup.m';
% ** Model selection **
models= {"cylinder", "pighead", 'PigHeadMeshWithElectrodes.vol', "mri-ct"};
model_sel= models{4};
fsTrue= 47.68;
fs= 48;
if model_sel== "cylinder"
    fmdl= ng_mk_cyl_models([2,1,.15],[32,0.5],.02);
    img = mk_image(ng_mk_cyl_models([2,1,.15],[32,0.5],.02), 1);
elseif model_sel== "pighead"
    maxsz= 0.05; maxh= 0.01; imgsize= [32 32];
    [fmdl, imdl]= mk_pighead_fmdl(maxsz, maxh, imgsize);
elseif model_sel== "mri-ct"
    maxsz= 0.2; maxh= 2; imgsize= [64 64];
    [fmdl, imdl]= mk_pighead_fmdl(maxsz, maxh, imgsize);
end % end if

% Reference frames are synchronization spikes from EIT files
data_files= {'Airchill_27.11.15.eit',                  'EIT_11.2_Nativ_1_Schleuse.eit', 'EIT_11.2_Nativ_2_ZVK.eit', ...
             'EIT_16.12_Baseline_Hypertonic_1.eit',    'EIT_8.2_nativ_1.eit',           'EIT_8.2_nativ_2.eit', ...
             'EIT_9.2_Nativ_1.eit',                    'EIT_9.2._Nativ_2.eit',          'eit_basal_2_20ml_hypertonic.eit', ...
             'EIT_cerebral_Baseline_2_Hypertonic.eit', 'EIT_nativ_2_zvk_10.2.eit'};

ref_frames= [ [786, 0, 0, 1192, 0];                    [1227, 6371, -477, 524, 0];      [213, 5696, 524, 763, 5436]; ...
              [957, 4684, 763, 1192, 0];               [683, 6372, 2003, 2432, 4673];   [1, 4992, 939, 1654, 4277]; ...
              [248, 5241, 620, 1287, 3910];            [293, 5950, 572, 906, 4625];     [725, 3747, 0, 0, 0]; ...
              [1163, 9994, 0, 954, 0];                 [264, 6731, 477, 763, 4148] ];

% Frames of interest from preliminary analysis
fsel1= {1750:3500, 0, 0,         2450:2650, 2251:2475, 1800:3250, 700:1100,  1400:2050, 0, 0, 4975:5175};
fsel2= {0,         0, 2750:3150, 2100:2800, 2100:3500, 2050:3500, 2750:3150, 1900:2800, 0, 0, 3100:3400};

%%
% ** Mode selection **
modes= ["Show all slices", "Plot ROI", "Plot frames of interest", "fft", "plot sum data"];
do= modes(3);
for f= 1:length(data_files)
    % ** GET DATA BASED ON SYNCHRONIZATION SPIKES, STOP VENTILATION, SALINE
    % INJECTION, START VENTILATION **
    file= data_files{f};
    title_name= remove_underscores(file);
    vv= real(eidors_readdata(file));
    vf= lowpass(vv', 18, fsTrue)';
%     vf= bandpass(real(vv)', [0.5, 3], fstrue)';
    reference_frames= ref_frames(f, :);
    [firstSync, secondSync, stopVent, inject, startVent, firstFrame, lastFrame, vh]= get_reference_frames(vv, reference_frames, fs);
    imgr= inv_solve(imdl,vh,vf(:, firstFrame: lastFrame));
    % Set colour limits
    clim= max(imgr.elem_data(:));
    cmin= min(imgr.elem_data(:));
    imgr.calc_colours.ref_level= 0;
    imgr.calc_colours.lim= clim;
    
    meanFrame = real(mean(calc_slices(imgr), 3));
    minMF= min(meanFrame(:));
    maxMF= max(meanFrame(:));
    rangeMF= maxMF- minMF;
    meanFrame= meanFrame*(200/ rangeMF); % Put pixel value range in reasonable limits for mrfDenoise.
    meanFrameDenoised= mrfDeNoiseV3(meanFrame, 1);
    figure; imagesc(meanFrameDenoised);
    figure; imagesc(meanFrame);
    keyboard;

    
%     test= mk_stim_patterns(32,1,[0,4],[0,4])
    
    if do== "Show all slices"
        figure; plot(sum( vf(:,:) ,1) ); title(title_name); 
        hold on; plot([stopVent, stopVent], [0, 0.5]); plot([inject, inject], [0, 0.5]); plot([startVent, startVent], [0, 0.5]); legend('Data', 'Stop Ventilation' , 'Inject Saline', 'Start Ventilation'); hold off;
        imgr.show_slices.img_cols = 100;
        figure; show_slices(imgr); title(title_name);
    
    elseif do== "Plot ROI"
        acvOpt= struct;
        acvModes= ["show each", "show mean", "vs mean", "show ROC", "each mean vs all mean"];
        acvOpt.mode= acvModes(1);
        acvOpt.title= title_name;
        acvOpt.roiTitle= "Brain";
        acvOpt.pixels= (21:40);
        timePoints.stopVent= stopVent- firstFrame;
        timePoints.inject= inject-firstFrame;
        timePoints.startVent= startVent- firstFrame;
        timePoints.fsTrue= fsTrue;
%         timePoints.start= timePoints.inject;
%         timePoints.stop= timePoints.inject + 45*fs; % watch for 45 s after injection
        plotSegmentations(imgr, model_sel, timePoints, acvOpt);
        keyboard;
        close all;

    elseif do== "Plot frames of interest"
        if length(fsel1{f}) > 1
            imgs = calc_slices(imgr);
            figure;
            subplot(211)
            ctr = squeeze(imgs(48,round(linspace(18,32,5)),:))';
            xax = (1:length(ctr))/fsTrue;
            plot(xax,ctr);
            framesel= fsel1{f};
            framesel= round(linspace(framesel(1), framesel(length(framesel)), 40));
            if xax(framesel(end)+fs) < xax(end)
                xlim([xax(framesel(1)-fs), xax(framesel(end)+fs)]);
            end % end if
            yl = ylim;
            line([1;1]*xax(framesel),yl'*ones(size(framesel)),'Color',.7*[1,1,1]);
            line([1;1]*xax(inject), yl); % injection of saline
            set(gca,'YTickLabel',[]); ylim(yl);
            title(title_name);
            legend('left', 'center-left', 'center', 'center-right', 'right');
            subplot(212)
            imgr.get_img_data.frame_select = framesel;
            imgr.show_slices.img_cols = 10;
            show_slices(imgr);
            keyboard;
            print_convert(""+ title_name+ model_sel+ framesel(1)+ "-"+ framesel(length(framesel))+ ".jpg");
        end % end if
    
    elseif do== "fft"
        show_fft(vf, firstFrame, lastFrame, fsTrue);
        title(title_name+ " FFT");
        keyboard;
    
    elseif do== "plot sum data"
        dat= sum(vf, 1);
        xax = (1: size(dat, 2))/fsTrue;
        tsMin= min(dat);
        tsMax= max(dat);
        figure; hold on;
        plot(xax, dat);
        plot([stopVent/fsTrue, stopVent/fsTrue], [tsMin, tsMax]); 
        plot([inject/fsTrue, inject/fsTrue], [tsMin, tsMax]); 
        plot([startVent/fsTrue, startVent/fsTrue], [tsMin, tsMax]); 
        xlabel("Time (s)");
        ylabel("Voltage");
        title(title_name+ " 18 Hz Lowpass Total Boundary Voltage");
        legend('Total Boundary Voltage', 'Stop Ventilation' , 'Inject Saline', 'Start Ventilation');
        hold off;
        keyboard;
    end % end if
end % end for


function [firstSync, secondSync, stopVent, inject, startVent, firstFrame, lastFrame, vh]= get_reference_frames(vv, reference_frames, fs)
    firstSync= reference_frames(1);
    secondSync= reference_frames(2);
    if reference_frames(2)== 0
        lastFrame= size(vv,2)- 50;
    else
        secondSync= reference_frames(2); % could be 0 if NA
        lastFrame= secondSync- 50;
    end % end if
    if reference_frames(3)== 0
        stopVent= firstSync+ 200;
    elseif reference_frames(3)<= 0
        stopVent= firstSync+ reference_frames(3);
    else
        stopVent= firstSync+ reference_frames(3); % could be 0 if NA
    end % end if
    if reference_frames(4)== 0
        inject= firstSync+ 50; % only 1 file without inject ref
    else
        inject= firstSync+ reference_frames(4);
    end % end if
    if reference_frames(5)== 0
        if secondSync== 0
            startVent= lastFrame- 50;
        else
            startVent= secondSync- 50;
        end % end if
    else
        startVent= firstSync+ reference_frames(5);
    end % end if
    if stopVent< firstSync
        vh = mean(vv(:, fs:stopVent),2);
        firstFrame= firstSync+ 50;
    else
        vh = mean(vv(:, firstSync+50:stopVent),2);
        firstFrame= stopVent;
    end % end if
end % end function


function show_fft(data, start, stop, fs)
% start: The first frame to be included in the time series 
% stop:  The last frame to be included in the time series.
% fs:    The data framerate 
    figure;
    dr = real(data(:,start:stop));
    tax = linspace(0,fs,size(dr,2)); 
    semilogy(tax,abs(fft(sum(dr, 1)))); 
    xlim([0,fs/2]);
end % end function


function title_name= remove_underscores(file)
% Replace underscores with '-' and remove .eit suffix from file name
title_name= "";
for c= file(1: length(file)-4)
    if c== "_"
        title_name= title_name+ "-";
    else
        title_name= title_name+ c;
    end % end if
end % end for
end % end function
