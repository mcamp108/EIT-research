run 'myStartup.m';
cd 'C:\Users\Mark\Documents\GraduateStudies\LAB\WeightedRestraint\data\Mali Weighted Restraint';
clip= 50;
dirs= ls;

% model
[fmdl, imdl]= mk_weighted_restraint_model();

% lung roi
lung_roi= lung_segmentation();
mid_col= size(lung_roi, 2)/2;
lung_divide_idx= mid_col* size(lung_roi, 1) + 1;
is_lung= find(lung_roi);
left_lung_idx= is_lung(is_lung>= lung_divide_idx);
right_lung_idx= is_lung(is_lung< lung_divide_idx);
% frame numbers for start of ten 30-second intervals
ten_thrty_s_inc= linspace(1, round(270*47.6826),10);
n_frame_30_s= min(diff(ten_thrty_s_inc));

% parameters
header= {   'Time_s', 'ins_to_ins',               'exp_to_exp',               'ins_to_exp',...
                    'del_z_lung_max',           'del_z_lung_min',           'left_lung_filling',...          
                    'right_lung_filling',       'left_lung_median_value',	'right_lung_median_value',...  
                    'global_inhomogeneity',     'center_of_ventilation'};

figure('units','normalized','outerposition',[0 0 1 1]);

for i= 3: size(ls, 1)
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
    end % end for f
    
    files= {sref, pref, wref, wepos, epos};
    
    for j= 1:length(files)
        f= files{j};
        
        % Import data
        [dd,auxdata]= eidors_readdata(f);
        dd= dd(:, 1: min(size(dd, 2), 14300));
        auxdata.elec_impedance= auxdata.elec_impedance(:, 1: size(dd, 2));
        auxdata.t_abs= auxdata.t_abs(:, 1: size(dd, 2));
        auxdata.t_rel= auxdata.t_rel(:, 1: size(dd, 2));
        FR = 1e6./median(diff(auxdata.t_rel)); %framerate is median dif of time points/ 1000000 (convert to s)

%         % Inspect data quality
%         clf;
%         inspect_eit_elec_and_data({dd, auxdata}, imdl, 500); 
%         ax= gcf;
%         ylim_adjust= (0.03- (ax.Children(6).YLim(2)- ax.Children(6).YLim(1)))/2;
%         ax.Children(6).YLim(1)= ax.Children(6).YLim(1)- ylim_adjust;
%         ax.Children(6).YLim(2)= ax.Children(6).YLim(2)+ ylim_adjust;
        name= char(remove_underscores(f));
        suffix= ' clean breaths';
        sgtitle(horzcat(name, suffix));
        front= name(1:13);
        back= name(14:end);
%         cd 'C:\Users\Mark\Documents\GraduateStudies\LAB\WeightedRestraint\figures\data quality';
%         print_convert(horzcat(front, '_', num2str(i), '_', back, '.png'));
%         cd 'C:\Users\Mark\Documents\GraduateStudies\LAB\WeightedRestraint\data\Mali Weighted Restraint';
%         cd(folder);
%         keyboard;
        
        % Clean data
        msel= imdl.fwd_model.meas_select;
        mm = find(msel);
        if strcmp(f, '2019_08_07_P3_proneRef.eit')
            imdl_comp= compensate_bad_elec(f, imdl, 1000);
        else
            imdl_comp= imdl;
        end % end if
        use_data= real(dd(mm, :));
        
        % for lung component
        use_data= lowpass(use_data', 4, FR)';
        
        % trim filter edge artifacts
        use_data= use_data(:, clip:(size(use_data,2)-clip) );
        
        % Solve
        imgr= inv_solve(imdl_comp, mean(use_data, 2), use_data); % data 1 == referene frame, data 2== other frames in time series.
        if j==1
            % compile images from all conditions into a single image
            imgr_copy= imgr;
            % 5 conditions, max 10 images per condition.
            imgr_copy.elem_data= zeros(size(imgr.elem_data, 1), 50);
            imgr_copy.show_slices.img_cols= 10;
        end % end if
        slices= calc_slices(imgr, [inf,inf, 1]);
%         for plane= [0.85, 1, 1.15]
%             slices= calc_slices(imgr, [inf,inf, plane]);
%             save(horzcat(f(1:end-4), '_reconstr_data_zplane',num2str(plane), '.mat'), 'slices');
%         end % end for        
        
        % Breath boundaries
        tbv= sum(use_data, 1);
        breaths= findBreaths(tbv);
        end_in= breaths.insIdx;
        end_ex= breaths.expIdx;
        
        % apply breath selection
        b= length(end_in)- 1;
        tidal_volume= abs(tbv(end_ex(1:b))- tbv(end_in(1:b)));
        ins_to_ins= (diff(end_in)./FR)';
        exp_to_exp= (diff(end_ex)./FR)';
        ins_to_exp= (abs((end_ex(1:b)- end_in(1:b))./FR))';
        criteria= [ins_to_ins, exp_to_exp, ins_to_exp];
        reject_breath= zeros(length(ins_to_ins), 3);
        
        for k= 1:3
            testing= criteria(:, k);
            qntls= quantile(testing, 5);
            reject_breath(:,k)= (testing< qntls(1)) + (testing> qntls(end));
        end % end for k
        
        % calculate total boundary voltage
%         keep= find(sum(reject_breath, 2)==0);
        keep= find(reject_breath(:,3)==0);
        ins_to_ins= ins_to_ins(keep);
        exp_to_exp= exp_to_exp(keep);
        ins_to_exp= ins_to_exp(keep);
        tidal_volume= tidal_volume(keep);
        end_in= end_in(keep)';
        end_ex= end_ex(keep)';
        breath_start_stop= [end_in, end_ex];
        Time= min(breath_start_stop, [], 2)/FR;
        
        % plot breath boundaries
        plot(tbv, 'k');hold on;
        for k= 1:length(keep)
            x= sort([end_in(k), end_ex(k)]);
            x= x(1):x(2);
            xax= min(x):max(x);
            plot(xax, tbv(xax), 'm');
        end % end for k
        ylabel('Voltage');xlabel('Frame'); legend('total boundary voltage', 'retained data')
%         keyboard;
        cd 'C:\Users\Mark\Documents\GraduateStudies\LAB\WeightedRestraint\figures\data quality';
        print_convert(horzcat(front, '_', num2str(j), '_breaths_', back, '.png'));
        cd 'C:\Users\Mark\Documents\GraduateStudies\LAB\WeightedRestraint\data\Mali Weighted Restraint';
        cd(folder);     
        hold off;
        
        del_z_lung_max= zeros(k, 1);
        del_z_lung_min= zeros(k, 1);
        left_lung_area= zeros(k, 1); % number of pixels for lung in breath
        right_lung_area= zeros(k, 1);
        left_lung_filling= zeros(k, 1); % total value of pixels for lung in breath
        right_lung_filling= zeros(k, 1);
        left_lung_median_value= zeros(k, 1); 
        right_lung_median_value= zeros(k, 1);
        global_inhomogeneity= zeros(k, 1);
        center_of_ventilation= zeros(k, 1);
        
        % lung images
%         lung_imgs= slices(:,:,end_in)- slices(:,:,end_ex);
        lung_imgs= imgr.elem_data(:,end_in)- imgr.elem_data(:,end_ex);
        j_idx= (j-1)* 10 + 1;
        if (j<3) || (j==5) % conditions sref, pref or epos. Show a single image, the average of 30 seconds of all breaths            
%                 lung_img= mean(lung_imgs, 3);
            imgr_copy.elem_data(:,j_idx)= mean(lung_imgs, 2);
        else % conditions wref or wepos. Show 10 images, each the average of 30 seconds of breathing
            for m= 1:length(ten_thrty_s_inc)
                idx= find( sum([(breath_start_stop> ten_thrty_s_inc(m)) , (breath_start_stop<= ten_thrty_s_inc(m)+ n_frame_30_s)], 2)==4 );
                if isempty(idx)
%                     disp("Something's up with Jack.");
%                         lung_img(:,:,m)= zeros(32,32)./0;
                else
%                     lung_img(:,:,m)= mean(lung_imgs(:,:,idx), 3);
                    imgr_copy.elem_data(:, j_idx+ m-1)= mean(lung_imgs(:,idx), 2);                        
                end % end if
            end % end for m                
        end % end if j
        
%         clf; show_slices(lung_img); colorbar; title(horzcat(name,' reconstruction(s)'));
%         print_convert(horzcat(front, '_', num2str(j), '_breaths_', back, '.png'));
        
        for k= 1: length(keep)
            lungs= slices(:,:,end_in(k))- slices(:,:,end_ex(k));
            % 1. maximum delta Z determined for lungs
            del_z_lung_min(k)= min(lungs,[],'all');
            del_z_lung_max(k)= max(lungs,[],'all')- del_z_lung_min(k);            

            % calculate lung filling
            left_lung_filling(k)= sum(lungs(left_lung_idx)); % sum of voltage difference in left lung
            right_lung_filling(k)= sum(lungs(right_lung_idx)); % sum of voltage difference in right lung
            left_lung_median_value(k)= median(lungs(left_lung_idx));
            right_lung_median_value(k)= median(lungs(right_lung_idx));
            vals= lungs(is_lung);
            med_lung= median(vals);
            
            % calculate global inhomogeneity index
            global_inhomogeneity(k)= sum(abs(vals - med_lung))/ sum(vals);
            
            % calculate coefficient of variation
            center_of_ventilation(k)= std(vals)/ mean(vals);
        end % end for k
        
        out= [      Time, ins_to_ins,               exp_to_exp,               ins_to_exp,...
                    del_z_lung_max,           del_z_lung_min,           left_lung_filling,...          
                    right_lung_filling,       left_lung_median_value,	right_lung_median_value,...  
                    global_inhomogeneity,     center_of_ventilation];
            
        save_name= horzcat(f(1:end-4), '_features.csv');
        table_out = array2table(out, 'VariableNames', header);
%         writetable(table_out, save_name);    
    end % end for j
    show_slices(imgr_copy, [inf,inf, 1]);
    title(horzcat(front, 'sref, pref, wref, wepos, epos'));
    print_convert(horzcat(front, '_', 'allconditions_breaths', '.png'));
    cd ../   
end % end for i
