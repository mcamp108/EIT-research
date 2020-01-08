run 'myStartup.m';
cd 'C:\Users\Mark\Documents\GraduateStudies\LAB\EIT-restraint\zzMC\data\Mali Weighted Restraint';
figure('units','normalized','outerposition',[0 0 1 1]);clf;
global glob_ref;
global stats;
global end_ex;
global end_in;
global all_seq_av_breath;
global imgr_idx_list;
global av_frm_s_inc; % starting frames of we_segments number of evenly spaced time increments
global n_frame_s; % minimum number of frames in time increments
global figure_img;
global cmap;

we_segments= 5; % number of images to produce forwreft and wepos
sep= 7; % number of pixels between figure images
av_frm_s_inc= linspace(1, round(270*47.6826), we_segments); % starting frames of we_segments number of evenly spaced time increments
n_frame_s= min(diff(av_frm_s_inc)); % minimum number of frames in time increments
imgr_idx_list= [1 2 3 8 13];

dirs= ls;

[fmdl, imdl]= mk_weighted_restraint_model(); % model

% configure colormap
calc_colours('cmap_type', 'blue_black_red');
cmap=colormap;
cmap= [[1, 1, 1]; flipud(cmap(2:end, :))];
cmap(2,:)= [0 0 0];
cmap(3,:)= [1 1 1]*0.5;

lung_roi= lung_segmentation(); % lung roi
mid_col= size(lung_roi, 2)/2;
lung_divide_idx= mid_col* size(lung_roi, 1) + 1;
is_lung= find(lung_roi);
left_lung_idx= is_lung(is_lung>= lung_divide_idx);
right_lung_idx= is_lung(is_lung< lung_divide_idx);
recordings= {'sref', 'pref', 'wref', 'wepos', 'epos'};

for i= 3: size(ls, 1) % for each participant
    stats=struct;
    all_seq_av_breath=struct;
    clf;
    folder= dirs(i, :);
    
    if contains(folder, '2019_')
        files= load_wr_data(folder);
    else
        continue
    end % end for
    
    for j= [2,3,4,5,1] % for each condition
        f= files{j};
        name= char(remove_underscores(f));
        suffix= ' clean breaths';        
        front= name(1:13);
        back= name(14:end);
        imgr_idx= imgr_idx_list(j);
%--------------------------------------------------------------------------        
        % Data pre-processing and breath selection
        [use_data, imdl_comp,FR]= wr_pp(f, imdl);
        tbv= sum(use_data, 1);
        [end_in, end_ex, Ti_by_Tt, BF]= select_breaths(tbv, FR);  
%--------------------------------------------------------------------------        
        % plot selected breath boundaries
%         plot_breath_boundaries(name, suffix, tbv, end_in, end_ex, front, back, j);
%--------------------------------------------------------------------------        
        glob_ref= mean(use_data(:, end_ex(1,1):end_ex(2,1)), 2); % use mean of first breath as reference
%         if j==2
%             glob_ref= mean(use_data, 2);

%         end % end if j
        % Image reconstruction
        imgr= inv_solve(imdl_comp, glob_ref, use_data); % solve
        
        if j==2
            figure_img= imgr; % compile images from all conditions into a single image
            figure_img.elem_data= zeros(size(imgr.elem_data, 1), 2*we_segments+3); % 5 conditions, max we_segments images for 2 conditions, 1 image for 3 conditions.
            figure_img.show_slices.img_cols= 2*we_segments+3;
            figure_img.show_slices.sep= sep;
            imgFRC= figure_img;
            imgTV= figure_img;
        end % end if

        % lung images... inspiration minus mean of flanking expirations        
        imgFRC.elem_data= (imgr.elem_data(:,end_ex(1,:)) + imgr.elem_data(:,end_ex(2,:)))./2;
        imgTV.elem_data= imgr.elem_data(:,end_in)- imgFRC.elem_data;
        
        calc_img_stats(imgTV, j, 'TV');
        calc_img_stats(imgTV, j, 'fig');
        calc_img_stats(imgFRC, j, 'FRC');
        calc_time_stats(BF, j, 'BF');
        calc_time_stats(Ti_by_Tt, j, 'FIT');
        calc_time_stats(tbv, j, 'av_breath');
    end % end for j
    
% --------------------------------------------------------------------------
    % Reconstructed images figure
    levels= [inf,inf,1.25; inf,inf,1; inf,inf,0.75];
    clim= max(figure_img.elem_data(:));
    figure_img.calc_colours.ref_level= 0;
    figure_img.calc_colours.clim= clim;
    final_img= calc_slices(figure_img, levels);
    final_img(isnan(final_img)) = 0;
    cntr_of_vnt= calc_cntr_of_vnt(final_img, lung_roi);
    final_img=calc_colours(calc_slices(figure_img, levels));
    cntr_of_vnt= reshape(cntr_of_vnt,2,size(final_img, 3))';

    for j =1:size(final_img, 3)
        cntr_x= cntr_of_vnt(j,1);
        cntr_y= cntr_of_vnt(j,2);
        % COV x coordinate and line
        final_img( (cntr_y- min(4,cntr_y-1)  : min(cntr_y+4,32)), cntr_x, j)= 3;
        % COV y coordinate and line
        final_img( cntr_y, :, j )= 3; 
    end % end for j

    final_img2= mk_mosaic(final_img(6:27,:,:), [sep, 3], [], 13);
    final_img2(isnan(final_img2)) = 0;
    boundaries= zeros(1,4); boundaries(1)= 32+ median(1:sep); boundaries(2:end)= boundaries(1)+ [1,6,11] * (32+ sep);
    final_img2(:, boundaries)= 2; % black

%--------------------------------------------------------------------------
    % Normalize stats to unweighted prone
    
    % FRC (mean end ex) as a difference from prone FRC baseline, translated
    % to relative volume by taking it as a ratio of prone TV (mean end ex -
    % mean end ins). 
    stats.FRCSD= stats.FRCSD/stats.FRC(2);
    stats.FRC= (stats.FRC - stats.FRC(2)) ./ stats.TV(2); 
    stats.TVSD= stats.TVSD./ stats.TV(2);
    stats.TV= stats.TV./ stats.TV(2);
    stats.BFSD= stats.BFSD./ stats.BF(2);
    stats.BF= stats.BF./ stats.BF(2);
%% --------------------------------------------------------------------------
    % Make Figure
    paper_figure('reconst', final_img2, i-2);
    paper_figure('breath', all_seq_av_breath.figcols, i-2);
    paper_figure('stat', stats, 1);
% --------------------------------------------------------------------------
    % Export    
    cd 'C:\Users\Mark\Documents\GraduateStudies\LAB\EIT-restraint\zzMC\figures\paper';
    saveas(gcf, horzcat(front, '.svg'));
    
    cd 'C:\Users\Mark\Documents\GraduateStudies\LAB\EIT-restraint\zzMC\data\Mali Weighted Restraint\features';
    header= {'BF','BFSD','BFN', 'FRC','FRCSD','FRCN', 'TV','TVSD','TVN'};
    out=[];
    for h =1:length(header)
        out= [out, stats.(header{h})];
    end
    save_name= horzcat('P',num2str(i-2), '_features.csv');
    table_out = array2table(out, 'VariableNames', header);
    writetable(table_out, save_name);
    
    cd 'C:\Users\Mark\Documents\GraduateStudies\LAB\EIT-restraint\zzMC\data\Mali Weighted Restraint';
end % end for i

% -------------------------------------------------------------------------
%                       S U B F U N C T I O N S
% -------------------------------------------------------------------------

function paper_figure(command, data, num)
global cmap;
B= 0;
W= 0.95;
H= 0.13; % height of each plot
L= 0.01; % left margin of reconstruction and stats plots
BL= 0.0043; % left margin of first breath plot
R= 0.98;
bpb= linspace(BL, W+0.012, 13+1);
spc= W/13;
shift=0.04;    

if strcmp(command, 'breath')
    mxsz= 0;
    datamx= 0;
    len_data= length(data);
    data_sz= zeros(len_data, 1);
    for i= 1:len_data
        data_sz(i)= size(data{i}, 2);
        mxsz= max( data_sz(i),  mxsz);
        datamx= max( max(data{i},[],'all') ,  datamx);
    end % end for len
    start= floor((mxsz- data_sz)./ 2) + 1;
    for num= 1:len_data
        ax= axes('Position',[bpb(num) 0.003+3.5*H spc H],'Box','on');
        s= start(num);
        e= start(num)+data_sz(num)-1;
        hold on;
        for j=1:size(data{num}, 1) - 1
            plot( (s:e), data{num}(j,:), 'Color',[0 0 0]+0.75, 'linewidth', 1 ); % grey
        end % end for
        plot( (s:e), data{num}(j+1,:), 'm', 'linewidth', 2 ); % magenta
        ax.XTickLabel=[];ax.YTickLabel=[];ax.XTick=[];ax.YTick=[];
        ax.XLim= ([1 mxsz]);
        ax.YLim= ([0 datamx]);
    end % end for num

elseif strcmp(command, 'stat')
    LW= 3; 
    ax1=axes('Position',[L  B+4.75*H+shift   W   2*H],'Box','on'); 
    plot(   data.FRC,       'linewidth', LW); hold on;
    plot(   data.TV,        'linewidth', LW);
    plot(   data.BF,        'linewidth', LW);
    ax1.XLim=[0.5 13.5];    ax1.XTick=1:13;     ax1.XLabel.FontSize= 20;
    ax1.XTickLabels={'SR 0-2 mins', 'PR 0-2 mins', 'PW 0-1 min', 'PW 1-2 mins', 'PW 2-3 mins', 'PW 3-4 mins', 'PW 4-5 mins', 'PWE 0-1 min', 'PWE 1-2 mins', 'PWE 2-3 mins', 'PWE 3-4 mins', 'PWE 4-5 mins', 'PP 0-2 mins'};
    ax1.YAxisLocation= 'right';
    xtickangle(30);
    legend('Normalized Delta FRC',  'Normalized V_T', 'Normalized F_B');

elseif strcmp(command, 'reconst')
    fg1=axes('Position',[L B+H-shift R H*3]);
    image(data);colorbar;
    axis image
    axis off
    axis equal
    axis tight
    fg1.Colormap= cmap;
    sgtitle(horzcat('Participant ', num2str(num)));
end % end if

end % end function

% -------------------------------------------------------------------------
function [end_in, end_ex, Ti_by_Tt, BF]= select_breaths(tbv, FR)

    % Breath boundaries
    breaths= find_breaths(tbv);
    end_in= breaths.ins_idx;
    end_ex= breaths.exp_idx;
    exp_to_exp= (diff(end_ex, 1))';
    Te= ((end_ex(2,:)- end_in))';
    Ti= (abs(end_ex(1,:)- end_in))';
    boundaries= Te+ Ti;
    end_ex_dif= abs(diff(tbv(end_ex), 1))';

    qntls= quantile(boundaries, 7);
    reject_breath= (boundaries> qntls(end)) + (boundaries< qntls(1));
    qntls= quantile(end_ex_dif, 7);
    reject_breath= reject_breath+ (end_ex_dif> qntls(6));
    qntls= quantile(Ti, 7);
    reject_breath= reject_breath+ (Ti> qntls(end));
    keep= find(reject_breath==0);
    
    end_in= end_in(keep);
    end_ex= end_ex(:,keep);
    Ti= Ti(keep);
    Te= Te(keep);
    exp_to_exp= exp_to_exp(keep);
    
    Ti_by_Tt= Ti ./ (Ti + Te);
    BF= 60./ (exp_to_exp./ FR); % instantaneous breaths per minute

end % end function

% -------------------------------------------------------------------------
function all_seq_av_breath= calc_average_breath(tbv, end_in, end_ex)

n_breaths= size(end_ex, 2);
b_lens= abs(diff(end_ex)); % lengths of breaths
b_centers= end_in- (end_ex(1,:)) + 1; % index relative to 1 of where inspiration peak occurs
b_start= repmat(max(b_centers), 1, n_breaths)+ 1- b_centers;
all_seq_av_breath= zeros(n_breaths+1, max(b_start+b_lens));

for i= 1:n_breaths
    b_len= b_lens(i);
    x1= end_ex(1,i);
    x2= end_ex(2,i);
    y= tbv(x1:x2);
    y= y-min(y);
    xax= b_start(i):b_len+b_start(i);
    all_seq_av_breath(i,xax)= y;
end % end for

av_breath= mean(all_seq_av_breath, 1);
all_seq_av_breath(i+1,:)= av_breath;

end % end function

% -------------------------------------------------------------------------
function cntr_of_vnt= calc_cntr_of_vnt(img, lung_roi)
% img is a 4D array. Output will be 2 x 13 x 3
% remove image components not related to respiration
img_= -(img.*lung_roi); % invert pixel values
row_pix= sum(lung_roi, 2);
row_pix(row_pix==0)=1; % avoid NaN in division
col_pix= sum(lung_roi, 1);
col_pix(col_pix==0)=1; % avoid NaN in division
% 2 coordinates per image, images/row, number of rows
n_imgs= size(img_, 3);
n_rows= size(img_, 4);
cntr_of_vnt= zeros(2, n_imgs, n_rows); 

for h=1:n_rows
    for i= 1:n_imgs
        this_img= img_(:,:,i,h);
        cntr_col= sum(this_img, 1) ./ col_pix; % take as weighted average
        cntr_row= sum(this_img, 2) ./ row_pix; % same
        cntr_of_vnt(1,i,h) = round(sum(cntr_col.* (1:32), 'all') ./ sum(cntr_col));     % X coor
        cntr_of_vnt(2,i,h) = round(sum(cntr_row.* (1:32)', 'all') ./ sum(cntr_row));    % Y coor
        if isnan(cntr_of_vnt(1,i,h)) || isnan(cntr_of_vnt(2,i,h))
            keyboard;
        end % end if
    end % end for i
end % end for h
cntr_of_vnt=abs(cntr_of_vnt);
end % end function

% -------------------------------------------------------------------------
function files= load_wr_data(folder)
    
while strcmp(folder(end), ' ')
    folder= folder(1:end-1);
end % end while

hdr= length(horzcat(folder, '_'))+ 1;
cd(folder);
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
cd ../.
end % end function

% -------------------------------------------------------------------------
function plot_breath_boundaries(name, suffix, tbv, end_in, end_ex, front, back, j)

n_breaths= length(end_in);
fg1= figure(1);clf;
sgtitle(horzcat(name, suffix));
plot(tbv, 'k');
hold on
for k= 1:n_breaths
    x= sort([end_in(k), end_ex(:,k)']);
    x= x(1):x(end);
    plot(x, tbv(x), 'm');
end % end for k
plot(end_in, tbv(end_in), 'ob');
plot(end_ex, tbv(end_ex), 'or');
ylabel('Voltage'); xlabel('Frame'); legend('total boundary voltage', 'retained data');

hold off;
starting_dir= cd;
cd 'C:\Users\Mark\Documents\GraduateStudies\LAB\EIT-restraint\zzMC\figures\data quality\breath selection';
saveas(fg1, horzcat(front, '_', num2str(j), '_breaths_', back, '.svg'));
cd(starting_dir);

end % end function

% -------------------------------------------------------------------------
function [use_data, imdl_comp,FR]= wr_pp(f, imdl)
    clip= 50;    
    [dd,auxdata]= eidors_readdata(f);
    dd= dd(:, 1: min(size(dd, 2), 14300));
    auxdata.elec_impedance= auxdata.elec_impedance(:, 1: size(dd, 2));
    auxdata.t_abs= auxdata.t_abs(:, 1: size(dd, 2));
    auxdata.t_rel= auxdata.t_rel(:, 1: size(dd, 2));
    FR = 1e6./median(diff(auxdata.t_rel)); %framerate is median dif of time points/ 1000000 (convert to s)

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
end % end function

% -------------------------------------------------------------------------
function calc_time_stats(data, j, param)
global stats;
global end_in;
global end_ex;
global imgr_idx_list;
global av_frm_s_inc;
global n_frame_s;
global all_seq_av_breath;

param_list= {'BF', 'FIT', 'av_breath'};
if ~contains(param_list, param)
    disp('Unrecognized statistical parameter. Accepted parameters are: BF, FIT, av_breath');
    return
end % end if

if strcmp(param, 'av_breath')
    do_breaths=true;
else
    do_breaths=false;
    paramSD= horzcat(param, 'SD');
    paramN= horzcat(param, 'N');
end

stat_idx= imgr_idx_list(j);

if j<3 || j==5
    if do_breaths
        all_seq_av_breath.figcols{stat_idx}= calc_average_breath(data, end_in, end_ex);
    else
        stats.(param)(stat_idx,1)=  mean(data);
        stats.(paramSD)(stat_idx,1)=  std(data);
        stats.(paramN)(stat_idx,1)=  length(data);
    end
elseif j==3 || j==4
    for m= 1:length(av_frm_s_inc)
        begins_after= end_ex(1,:) > av_frm_s_inc(m);
        ends_before_at= end_ex(2,:) <= (av_frm_s_inc(m)+ n_frame_s);
        img_idx= find((begins_after + ends_before_at)==2);
        if ~isempty(img_idx)
            if do_breaths
                all_seq_av_breath.figcols{stat_idx}= calc_average_breath(data, end_in(img_idx), end_ex(:,img_idx));
            else
                stats.(param)(stat_idx,1)=  mean(data(img_idx));
                stats.(paramSD)(stat_idx,1)=  std(data(img_idx));
                stats.(paramN)(stat_idx,1)=  length(data(img_idx));
            end
        end % end if
        stat_idx= stat_idx + 1;
    end % end for m
else
    disp("Unrecognized experiment phase.")
end % end if

end

% -------------------------------------------------------------------------
function calc_img_stats(img, j, param)
global stats;
global end_ex;
global imgr_idx_list;
global av_frm_s_inc;
global n_frame_s;
global figure_img;
param_list= {'FRC', 'TV', 'fig'};
if ~contains(param_list, param)
    disp('Unrecognized statistical parameter. Accepted parameters are: FRC, TV, fig');
    return
end % end if

if strcmp(param, 'FRC') || strcmp(param, 'TV')
    img_slices= calc_slices(img, [inf inf 1]);
    img_slices(isnan(img_slices))= 0;
    img_slices= img_slices .* lung_segmentation();
    img_slices= -img_slices;
    img_sums= squeeze(sum(img_slices, [1 2]));
    paramSD= horzcat(param, 'SD');
    paramN= horzcat(param, 'N');
end % end if

stat_idx= imgr_idx_list(j);

if j<3 || j==5
    if strcmp(param, 'fig')
        figure_img.elem_data(:,stat_idx)= mean(img.elem_data, 2);
    else
        stats.(param)(stat_idx,1)= mean(img_sums);    
        stats.(paramSD)(stat_idx,1)= std(img_sums); 
        stats.(paramN)(stat_idx,1)= length(img_sums);
    end % end if
    
elseif j==3 || j==4
    for m= 1:length(av_frm_s_inc)
        begins_after= end_ex(1,:) > av_frm_s_inc(m);
        ends_before_at= end_ex(2,:) <= (av_frm_s_inc(m)+ n_frame_s);
        img_idx= find((begins_after + ends_before_at)==2);
        if ~isempty(img_idx)
            if strcmp(param, 'fig')
                figure_img.elem_data(:,stat_idx)= mean(img.elem_data(:,img_idx), 2);
            else
                stats.(param)(stat_idx,1)= mean(img_sums(img_idx));
                stats.(paramSD)(stat_idx,1)= std(img_sums(img_idx)); 
                stats.(paramN)(stat_idx,1)= length(img_sums(img_idx));
            end % end if
        end % end if
        stat_idx= stat_idx + 1;
    end % end for m
else
    disp("Unrecognized experimental phase.")
end % end if

end % end function