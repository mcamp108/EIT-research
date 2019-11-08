%% MODEL
run('myStartup.m');
[fmdl, imdl]= mk_weighted_restraint_model();
%% DATA
cd 'C:\Users\Mark\Documents\GraduateStudies\LAB\WeightedRestraint\data\Mali Weighted Restraint';
clip= 50;
dirs= ls;

for i= 3: size(ls, 1)
    folder= dirs(i, :);
    while strcmp(folder(end), ' ')
        folder= folder(1:end-1);
    end % end while
    cd(folder);
    hdr= length(horzcat(folder, '_'))+ 1;
    files= ls;
    for f= 3: size(files, 1)
        file= files(f, :);
        while strcmp(file(end), ' ')
            file= file(1:end-1);
        end % end while
        
        if strcmp(file(hdr:hdr+8), 'proneRef.')
            pref= file;
        elseif strcmp(file(hdr:hdr+13), 'proneRefWeight')
            wref= file;
        elseif strcmp(file(hdr:hdr+16), 'standingReference')
            sref= file;
        elseif strcmp(file(hdr:hdr+18), 'proneWeightExercise')
            wepos= file;
        elseif strcmp(file(hdr:hdr+20), 'proneNoWeightExercise')
            epos= file;
        else
            disp("Unrecognized file name!")
        end % end if
    end % end for
    files= {sref, pref, wref, wepos, epos};
    
    for f= files
        f= f{1};
        [fmdl, imdl]= mk_weighted_restraint_model();
        
        % Inspect data quality
%         inspect_eit_elec_and_data(f, imdl);
%         keyboard;
        
        % Import data
        [dd,auxdata]= eidors_readdata(f); 
        FR = 1e6./median(diff(auxdata.t_rel)); %framerate is median dif of time points/ 1000000 (convert to s)
        
        % Clean data
        msel= imdl.fwd_model.meas_select;
        mm = find(msel);
        imdl_comp= compensate_bad_elec(f, imdl);
        use_data= real(dd(mm, :));
        
        % for removing cardiac component
        % find heart rate
        tax = linspace(0,FR,size(use_data,2)); 
        freq= abs(fft(sum(use_data)- mean(use_data, 1)));
        greq= movmedian(freq, 6);
        zone= find(((tax>0.9) + (tax<2))==2);
        heart_rate= min(tax(greq==(max(greq(zone)))));
        
        % apply bandpass +/- 20% of heart rate
        cardiac_data= bandpass(use_data', [heart_rate*0.8/(FR/2), heart_rate*1.2/(FR/2)])';
        
        % for lung component
        use_data= lowpass(use_data', 5, FR)';
        
        % trim filter edge artifacts
        use_data= use_data(:, clip:(size(use_data,2)-clip) );
        cardiac_data= cardiac_data(:, clip:(size(cardiac_data,2)-clip) );
        
        % Solve
        imgr= inv_solve(imdl_comp, mean(use_data, 2), use_data); % data 1 == referene frame, data 2== other frames in time series.
        imgc= inv_solve(imdl_comp, mean(cardiac_data, 2), cardiac_data);
        
        % BREATH BOUNDARIES 
        tbv= sum(use_data, 1);
        breaths= findBreaths(tbv);
        
        if length(breaths.insIdx) > length(breaths.expIdx)
            breaths.insIdx= breaths.insIdx(1: length(breaths.expIdx));
        else
            breaths.expIdx= breaths.expIdx(1: length(breaths.insIdx));
        end % end if
        end_in= breaths.insIdx;
        end_ex= breaths.expIdx;
        
        % calculate tidal volume
        tidal_volume= abs(tbv(end_ex)- tbv(end_in));
        
        % calculate images for each breath
        tidal_variation= zeros(32,32, length(end_ex));
        GI_mat= zeros(length(end_in), 11);
        CV_mat= GI_mat;
        for j= 1: length(end_in)
            vd= use_data(:,end_in(j));
            vh= use_data(:,end_ex(j));
            imgr= inv_solve(imdl_comp, vd, vh);
            
            % subtract cardiac component
            vcd= cardiac_data(:,end_in(j));
            vch= cardiac_data(:,end_ex(j));
            imgc= inv_solve(imdl_comp, vcd, vch);   
            tidal_variation(:,:,j)= calc_slices(imgr, [inf,inf,1]);
            
            
            % calculate coefficient of variation
            % calculate global inhomogeneity index
            breath= calc_slices(imgr, [inf,inf,1]);
            heart= calc_slices(imgc, [inf,inf,1]);
            heart_roi= breath> (max(heart)*0.5);
            slice_max= max(breath(:));
            thresholds= (5:5:95)* 0.01;
            for t= 1:length(thresholds)
                thresh= thresholds(t)* slice_max;
                lung_roi= breath> (slice_max* thresh);
                
                % make symmetrical about x axis
                lung_roi= lung_roi | fliplr(lung_roi);
                
%                 lung_roi= lung_roi- heart_roi;
                vals= breath((lung_roi));
                rm_nan_vals= vals(~isnan(vals));
                med_lung= median(rm_nan_vals);
                GI_mat(j, t)= sum(abs(rm_nan_vals - med_lung))/ sum(breath(~isnan(breath)));
                CV_mat(j, t)= std(rm_nan_vals)/ mean(rm_nan_vals);
            end % end for
            subplot(1, 2, 1); plot(GI_mat(j,:));
            hold on;
            subplot(1, 2, 2); plot(CV_mat(j,:));
        end % end for
        
    end % end for
        
        ins_to_ins= diff(end_in)./FR;
        exp_to_exp= diff(end_ex)./FR;
        ins_to_exp= abs((end_ex- end_in)./FR);
        
end % end for       
%         % functional ROI
%         slices= calc_slices(imgr, [inf,inf,1]);
%         slice_sd= std(slices, 0, 3);
%         slice_sd_vals= slice_sd(~isnan(slice_sd));
%         max_sd= max(slice_sd_vals);
%         roi= zeros(32, 32);
%         thresholds= (15:5:50)* 0.01;
%         for t= 1:length(thresholds)
%             thresh= thresholds(t);
%             roi= roi+ (slice_sd> (max_sd* thresh));
%         end % end for
%         imagesc(roi);
%         %%
%         ann_max= max(tbv);
%         ann_min= min(tbv);
%     subplot(3,2,[5 6])
%         plot(tbv); hold on; for i= 1:length(end_ex); plot([end_ex(i) end_ex(i)], [ann_min ann_max]); plot([end_in(i) end_in(i)], [ann_min ann_max]); end; hold off;
%     subplot(3,2,1)    
%         hold on; plot(ins_to_ins); hold off;
%     subplot(3,2,2)
%         hold on; plot(exp_to_exp); hold off;
%     subplot(3,2,3)
%         hold on; plot(ins_to_exp); hold off;
%     subplot(3,2,4)
%         hold on; plot(tidal_volume); hold off;
%     end % end for
%     
%     subplot(3,2,1);legend(["Standing Reference", "Prone Reference", "Weighted Prone Reference", "Weighted Posture Exercise", "Unweighted Posture Exercise"]);
%     subplot(3,2,2);legend(["Standing Reference", "Prone Reference", "Weighted Prone Reference", "Weighted Posture Exercise", "Unweighted Posture Exercise"]);
%     subplot(3,2,3);legend(["Standing Reference", "Prone Reference", "Weighted Prone Reference", "Weighted Posture Exercise", "Unweighted Posture Exercise"]);
%     subplot(3,2,4);legend(["Standing Reference", "Prone Reference", "Weighted Prone Reference", "Weighted Posture Exercise", "Unweighted Posture Exercise"]);


%%
figure; show_slices(imgr, [inf,inf,1]);
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

    cd ../;

