run 'myStartup.m';
% DATA
cd 'C:\Users\Mark\Documents\GraduateStudies\LAB\WeightedRestraint\data\Mali Weighted Restraint';
clip= 50;
dirs= ls;

for i= 4: size(ls, 1)
    folder= dirs(i, :);
    
    while strcmp(folder(end), ' ')
        folder= folder(1:end-1);
    end % end while
    
    if ~strcmp(folder, 'Other')
        cd(folder);
    else
        continue
    end % end for
    
    hdr= length(horzcat(folder, '_'))+ 1;
    files= ls;
    
    for f= 3: size(files, 1)
        file= files(f, :);
        while strcmp(file(end), ' ')
            file= file(1:end-1);
        end % end while
        if strcmp(file(end-3:end), '.eit')
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
        end % end if
    end % end for
    
    files= {sref, pref, wref, wepos, epos};
    
    for file= files
        f= file{1};
        [fmdl, imdl]= mk_weighted_restraint_model();
        
        % Import data
        [dd,auxdata]= eidors_readdata(f);
        dd= dd(:, 1: min(size(dd, 2), 14300));
        auxdata.elec_impedance= auxdata.elec_impedance(:, 1: size(dd, 2));
        auxdata.t_abs= auxdata.t_abs(:, 1: size(dd, 2));
        auxdata.t_rel= auxdata.t_rel(:, 1: size(dd, 2));
        FR = 1e6./median(diff(auxdata.t_rel)); %framerate is median dif of time points/ 1000000 (convert to s)

        % Inspect data quality
        inspect_eit_elec_and_data({dd, auxdata}, imdl); sgtitle(remove_underscores(f));
        keyboard;
        
        % Clean data
        msel= imdl.fwd_model.meas_select;
        mm = find(msel);
        imdl_comp= compensate_bad_elec(f, imdl);
        use_data= real(dd(mm, :));
        
        % for lung component
        use_data= lowpass(use_data', 5, FR)';
        
        % trim filter edge artifacts
        use_data= use_data(:, clip:(size(use_data,2)-clip) );
        
        % Solve
        imgr= inv_solve(imdl_comp, mean(use_data, 2), use_data); % data 1 == referene frame, data 2== other frames in time series.
        slices= calc_slices(imgr, [inf,inf,1]);
        save(horzcat(f(1:end-4), '_reconstr_data.mat'), 'slices');
    end % end for
    cd ../
end % end for
%         
%         % Breath boundaries
%         breaths= findBreaths(use_data);
%         end_in= breaths.insIdx;
%         end_ex= breaths.expIdx;
%         
%         % calculate total boundary voltage
%         tbv= sum(use_data, 1);
%         
%         % calculate images for each breath
% %         tidal_variation= zeros(32,32, length(end_ex));
%         
%         b= length(end_in);
%         % parameters
%         header= {   'ins_to_ins',               'exp_to_exp',               'ins_to_exp',...
%                     'del_z_lung_max',           'del_z_lung_min',           'left_lung_area',...
%                     'right_lung_area',          'left_lung_filling',        'right_lung_filling',...
%                     'left_lung_median_value',	'right_lung_median_value',  'global_inhomogeneity',...
%                     'center_of_ventilation'};
%                 
%         tidal_volume= abs(tbv(end_ex)- tbv(end_in));
%         ins_to_ins= diff(end_in)./FR;
%         exp_to_exp= diff(end_ex)./FR;
%         ins_to_exp= abs((end_ex- end_in)./FR);
%         del_z_lung_max= zeros(b, 1);
%         del_z_lung_min= zeros(b, 1);
%         left_lung_area= zeros(b, 1); % number of pixels for lung in breath
%         right_lung_area= zeros(b, 1);
%         left_lung_filling= zeros(b, 1); % total value of pixels for lung in breath
%         right_lung_filling= zeros(b, 1);
%         left_lung_median_value= zeros(b, 1); 
%         right_lung_median_value= zeros(b, 1);
%         global_inhomogeneity= zeros(b, 1);
%         center_of_ventilation= zeros(b, 1);
%         
%         for j= 1: b
%             % Get difference data for this breath
%             vd= use_data(:,end_in(j));
%             vh= use_data(:,end_ex(j));
%             imgr= inv_solve(imdl_comp, vd, vh);
%             lungs= calc_slices(imgr, [inf,inf,1]);
%             
%             % 1. maximum delta Z determined for lungs
%             del_z_lung_min(j)= min(lungs,[],'all');
%             del_z_lung_max(j)= max(lungs,[],'all')- del_z_lung_min(j);
%             
%             % Determine fROI cutoff and lung segmentation
%             thresh= 0.5;    
%             fROI_cutoff= del_z_lung_max(j)* thresh+ del_z_lung_min(j);
%             lung_roi= lungs> (fROI_cutoff); % lung segmentation
%             lung_roi_area_projection= sum(lung_roi, 1); % project segmentation to 1d vector
%             lung_val_projection= sum(lungs(lung_roi), 1); % ?
%             lung_boundaries= diff(find(lung_roi_area_projection)) > 1; 
%             % space in between lungs is col where sum is 0.
%             % find(lung_boundaries) returns the index of
%             % lung_roi_area_projection that is the right side of the left
%             % lung in the image, left side of right lung in the participant
%             
%             if sum(lung_boundaries)== 0 % Image is weird, simply divide it in half.
%                 mid_col= size(lung_roi, 2)/2;
%                 left_lung_area(j)= sum(lung_roi(:, mid_col+1:end), 'all'); % number of pixels in left lung
%                 right_lung_area(j)= sum(lung_roi(:, 1: mid_col), 'all'); % number of pixels in right lung
%                 lung_divide_idx= mid_col* size(lung_roi, 1) + 1;
%             else % If there is a clear boundary between lungs.
%                 left_lung_area(j)= sum(lung_roi_area_projection(find(lung_boundaries)+1: end)); % number of pixels in left lung
%                 right_lung_area(j)= sum(lung_roi_area_projection(1: find(lung_boundaries))); % number of pixels in right lung
%                 lung_divide_idx= find(lung_boundaries, 1)* size(lung_roi, 1) + 1;
%             end % end if
%             
%             % calculate lung filling and area
%             is_lung= find(lung_roi);
%             left_lung_idx= is_lung(is_lung>= lung_divide_idx);
%             right_lung_idx= is_lung(is_lung< lung_divide_idx);
% 
%             left_lung_filling(j)= sum(lungs(left_lung_idx)); % sum of voltage difference in left lung
%             right_lung_filling(j)= sum(lungs(right_lung_idx)); % sum of voltage difference in right lung
%             
%             left_lung_median_value(j)= median(lungs(left_lung_idx));
%             right_lung_median_value(j)= median(lungs(right_lung_idx));
%             vals= lungs([left_lung_idx; right_lung_idx]);
%             med_lung= median(vals);
%             
%             % calculate global inhomogeneity index
%             global_inhomogeneity(j)= sum(abs(vals - med_lung))/ sum(vals);
%             
%             % calculate coefficient of variation
%             center_of_ventilation(j)= std(vals)/ mean(vals);
%         end % end for
%         
%         out= [  [0;ins_to_ins'],                 [0;exp_to_exp'],               ins_to_exp',...
%                 del_z_lung_max,             del_z_lung_min,           left_lung_area,...
%                 right_lung_area,            left_lung_filling,        right_lung_filling,...
%                 left_lung_median_value,     right_lung_median_value,  global_inhomogeneity,...
%                 center_of_ventilation];
%             
%         save_name= horzcat(f(1:end-4), '_features.csv');
%         table_out = array2table(out, 'VariableNames', header);
%         writetable(table_out, save_name);    
%     end % end for
%     cd ../   
% end % end for
