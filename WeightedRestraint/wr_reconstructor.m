run 'myStartup.m';
cd 'C:\Users\Mark\Documents\GraduateStudies\LAB\EIT-restraint\zzMC\data\Mali Weighted Restraint';
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

% number of images to produce forwreft and wepos
we_segments= 5;
sep= 3; % number of pixels between figure images
% starting frames of we_segments number of evenly spaced time increments
av_frm_s_inc= linspace(1, round(270*47.6826), we_segments);
% minimum number of frames in time increments
n_frame_s= min(diff(av_frm_s_inc));

figure('units','normalized','outerposition',[0 0 1 1]);
% parameters
header= {   'Time_s', 'ins_to_ins',               'exp_to_exp',               'ins_to_exp',...
                    'del_z_lung_max',           'del_z_lung_min',           'left_lung_filling',...          
                    'right_lung_filling',       'left_lung_median_value',	'right_lung_median_value',...  
                    'global_inhomogeneity',     'center_of_ventilation'};

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
    imgr_idx= 1;
    
    for j= 1:length(files)
        f= files{j};
%--------------------------------------------------------------------------        
        % Data pre-processing and breath selection
        [dd,auxdata]= eidors_readdata(f);
        dd= dd(:, 1: min(size(dd, 2), 14300));
        auxdata.elec_impedance= auxdata.elec_impedance(:, 1: size(dd, 2));
        auxdata.t_abs= auxdata.t_abs(:, 1: size(dd, 2));
        auxdata.t_rel= auxdata.t_rel(:, 1: size(dd, 2));
        FR = 1e6./median(diff(auxdata.t_rel)); %framerate is median dif of time points/ 1000000 (convert to s)
        name= char(remove_underscores(f));
        suffix= ' clean breaths';        
        front= name(1:13);
        back= name(14:end);
        
        % Clean data
        msel= imdl.fwd_model.meas_select;
        mm = find(msel);

        if strcmp(f, '2019_08_07_P3_proneRef.eit')
            imdl_comp= compensate_bad_elec(f, imdl, 1000);
        else
            imdl_comp= imdl;
        end % end if
        
        use_data= real(dd(mm, :));
        use_data= lowpass(use_data', 1, FR)';
        use_data= use_data(:, clip:(size(use_data,2)-clip) ); % trim filter edge artifacts
        
        % apply breath selection
        tbv= sum(use_data, 1);
        [end_in, end_ex, ins_to_exp, exp_to_exp]= select_breaths(tbv, FR);
        n_breaths= length(end_in);   
%--------------------------------------------------------------------------        
        % plot selected breath boundaries
        figure(1);clf;
        sgtitle(horzcat(name, suffix));
        plot(tbv, 'k');
        hold on
        for k= 1:n_breaths
            x= sort([end_in(k), end_ex(:,k)']);
            x= x(1):x(end);
            plot(x, tbv(x), 'm');
        end % end for k      
        ylabel('Voltage'); xlabel('Frame'); legend('total boundary voltage', 'retained data')
        hold off;
        starting_dir= cd;
        cd 'C:\Users\Mark\Documents\GraduateStudies\LAB\EIT-restraint\zzMC\figures\data quality';
        print_convert(horzcat(front, '_', num2str(j), '_breaths_', back, '.png'));
        cd(starting_dir);     
        
%--------------------------------------------------------------------------        
        % Image reconstruction
        imgr= inv_solve(imdl_comp, mean(use_data, 2), use_data); % solve
        
        if j==1
            imgr_copy= imgr; % compile images from all conditions into a single image
            imgr_copy.elem_data= zeros(size(imgr.elem_data, 1), 2*we_segments+3); % 5 conditions, max we_segments images for 2 conditions, 1 image for 3 conditions.
            imgr_copy.show_slices.img_cols= 2*we_segments+3;
            imgr_copy.show_slices.sep= sep;
        end % end if
 
        % lung images... inspiration minus mean of flanking expirations
        lung_imgs= imgr.elem_data(:,end_in)- (imgr.elem_data(:,end_ex(1,:)) + imgr.elem_data(:,end_ex(2,:)))./2;
        
        if (j<3) || (j==5) % conditions sref, pref or epos. Show a single image, the average of 30 seconds of all breaths            
            imgr_copy.elem_data(:,imgr_idx)= mean(lung_imgs, 2);
            imgr_idx= imgr_idx+ 1;
            plot_average_breath(tbv, end_in, end_ex);
            av_b_title= horzcat(name, ' average breath');
            sgtitle(av_b_title);
            print_convert(horzcat(av_b_title, '.png'));
        else % conditions wref or wepos. Show we_segments number of images
            
            for m= 1:length(av_frm_s_inc)
                % check breath begins after beginning of time period and at or before time period ends.
                begins_after= end_ex(1,:) > av_frm_s_inc(m);
                ends_before_at= end_ex(2,:) <= (av_frm_s_inc(m)+ n_frame_s);
                idx= find((begins_after + ends_before_at)==2);
                
                if isempty(idx)
                else
                    imgr_copy.elem_data(:, imgr_idx)= mean(lung_imgs(:,idx), 2);
                    plot_average_breath(tbv, end_in(idx), end_ex(:,idx));
                    av_b_title= horzcat(name, ' average breath ', num2str(m-1), ' to ', num2str(m), ' minutes');
                    sgtitle(av_b_title);
                    print_convert(horzcat(av_b_title, '.png'));                
                end % end if
                
                imgr_idx= imgr_idx+ 1;

            end % end for m  
            
        end % end if j
        
%--------------------------------------------------------------------------
%         % Statistics
%         flow volume loops
        % changes in tidal volume
        end_ex_av= mean(tbv(end_ex), 1);
        TV= tbv(end_in)- end_ex_av;
        % changes in FRC
        
        % changes in breathing frequency
        BF= 1/ exp_to_exp;
        % changes in center of ventilation
        
%         global_inhomogeneity= zeros(n_breaths, 1);
%         center_of_ventilation= zeros(n_breaths, 1);    
%         
%         for k= 1:n_breaths
%             lungs= slices(:,:,end_in(k))- slices(:,:,end_ex(k));
%             % 1. maximum delta Z determined for lungs
%             del_z_lung_min(k)= min(lungs,[],'all');
%             del_z_lung_max(k)= max(lungs,[],'all')- del_z_lung_min(k);            
% 
%             % calculate lung filling
%             left_lung_filling(k)= sum(lungs(left_lung_idx)); % sum of voltage difference in left lung
%             right_lung_filling(k)= sum(lungs(right_lung_idx)); % sum of voltage difference in right lung
%             left_lung_median_value(k)= median(lungs(left_lung_idx));
%             right_lung_median_value(k)= median(lungs(right_lung_idx));
%             vals= lungs(is_lung);
%             med_lung= median(vals);
%             
%             % calculate global inhomogeneity index
%             global_inhomogeneity(k)= sum(abs(vals - med_lung))/ sum(vals);
%             
%             % calculate coefficient of variation
%             center_of_ventilation(k)= std(vals)/ mean(vals);
%         end % end for k
%--------------------------------------------------------------------------
%         % Export
%         out= [      Time, ins_to_ins,               exp_to_exp,               ins_to_exp,...
%                     del_z_lung_max,           del_z_lung_min,           left_lung_filling,...          
%                     right_lung_filling,       left_lung_median_value,	right_lung_median_value,...  
%                     global_inhomogeneity,     center_of_ventilation];
%             
%         save_name= horzcat(f(1:end-4), '_features.csv');
%         table_out = array2table(out, 'VariableNames', header);
%--------------------------------------------------------------------------        
    end % end for j
    figure(1);clf;
    levels= [inf,inf,0.75; inf,inf,1; inf,inf,1.25];
    clim= max(imgr_copy.elem_data(:));
    imgr_copy.calc_colours.ref_level= 0;
    imgr_copy.calc_colours.lim= clim;
    final_img= show_slices(imgr_copy, levels);
    boundaries= [1,2,7,12] * (32+ sep) - (median(1:sep)-1);
    final_img(:, boundaries)= 2; % black
    image(final_img); colorbar;
    axis image
    axis off
    axis equal
    axis tight
    figg= gcf;
    colormap jet
    figg.Colormap(2,:)= [0 0 0];
    figg.Colormap(1,:)= [1 1 1];
    sgtitle(horzcat(front, ' sref, pref, wref, wepos, epos'));
    cd 'C:\Users\Mark\Documents\GraduateStudies\LAB\EIT-restraint\zzMC\figures\reconstructions';
    print_convert(horzcat(front, '_', '3_zplanes_all_conditions', '.png'));
    
    cd 'C:\Users\Mark\Documents\GraduateStudies\LAB\EIT-restraint\zzMC\data\Mali Weighted Restraint';
end % end for i

function [end_in, end_ex, ins_to_exp, exp_to_exp]= select_breaths(tbv, FR)

    % Breath boundaries
    breaths= find_breaths(tbv);
    end_in= breaths.ins_idx;
    end_ex= breaths.exp_idx;
    exp_to_exp= (diff(end_ex, 1))';
    ins_to_exp= ((end_ex(2,:)- end_in))';
    exp_to_ins= (abs(end_ex(1,:)- end_in))';
    boundaries= ins_to_exp+ exp_to_ins;
    end_ex_dif= abs(diff(tbv(end_ex), 1))';

    qntls= quantile(boundaries, 7);
    reject_breath= (boundaries> qntls(end)) + (boundaries< qntls(1));
    qntls= quantile(end_ex_dif, 7);
    reject_breath= reject_breath+ (end_ex_dif> qntls(6));
    qntls= quantile(exp_to_ins, 7);
    reject_breath= reject_breath+ (exp_to_ins> qntls(end));
    keep= find(reject_breath==0);
    
    exp_to_exp= exp_to_exp(keep)./FR;
    ins_to_exp= ins_to_exp(keep)./FR;
    end_in= end_in(keep);
    end_ex= end_ex(:,keep);

    % plot
    plot(tbv, 'k');
    hold on;
    plot(end_in, tbv(end_in), 'ob');
    plot(end_ex, tbv(end_ex), 'or');

end % end function


function plot_average_breath(tbv, end_in, end_ex)

center_ins= true;
n_breaths= size(end_ex, 2);
b_lens= abs(diff(end_ex));
b_centers= end_in- (end_ex(1,:)) + 1;
b_start= repmat(max(b_centers), 1, n_breaths)+ 1- b_centers;
av_breath= zeros(1, max(b_start+b_lens));
figure(2);clf;
hold on;

for i= 1:n_breaths
    b_len= b_lens(i);
    x1= end_ex(1,i);
    x2= end_ex(2,i);
    y= tbv(x1:x2);
    y= y-min(y);
    if center_ins
        xax= b_start(i):b_len+b_start(i);
    else
        xax= 1:b_len+1;
    end
    plot(xax, y, 'Color',[0 0 0]+0.75, 'linewidth', 2);
    av_breath(xax)= av_breath(xax)+ y;
end % end for

av_breath= av_breath./ n_breaths;
plot(av_breath, 'm', 'linewidth', 8);
axis off
axis tight
hold off;

end % end function
