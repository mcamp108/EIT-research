run 'myStartup.m';
dataDir = 'E:\University\Masters\EIT-restraint\zzMC\data\phys_meas\EIT\';
saveDir = 'E:\University\Masters\EIT-restraint\zzMC\figures\paper\';
featDir = 'E:\University\Masters\EIT-restraint\zzMC\features\';
feat_save_name = sprintf('%sper_observation_features_selfRef2.csv', featDir);
% figure('units','normalized','outerposition',[0 0 1 1]); clf;
global stats;
global end_ex;
global end_in;
global all_seq_av_breath;
global imgr_idx_list;
global av_frm_s_inc; % starting frames of we_segments number of evenly spaced time increments
global n_frame_s; % minimum number of frames in time increments
global figure_img;
global cmap;
global we_segments;
global fs
referenceStyle = 'self';
fs = 47.6826;
we_segments     = 5; % number of simages to produce for wref and wepos
sep             = 7; % number of pixels between figure images
timeSeriesLength= 300; % shortest length of time series under investigation. This will set the time per segment.
av_frm_s_inc    = linspace(1, round(timeSeriesLength * fs), we_segments + 1); % starting frames of we_segments number of evenly spaced time increments
av_frm_s_inc    = av_frm_s_inc(1 : we_segments);
n_frame_s       = min(diff(av_frm_s_inc)); % minimum number of frames in time increments
[fmdl, imdl]    = mk_weighted_restraint_model(); % get model
cmap = config_cmap(); % configure colormap
lung_roi        = lung_segmentation(); % lung roi
mid_col         = size(lung_roi, 2)/2;
lung_divide_idx = mid_col* size(lung_roi, 1) + 1;
is_lung         = find(lung_roi);
left_lung_idx   = is_lung(is_lung>= lung_divide_idx);
right_lung_idx  = is_lung(is_lung< lung_divide_idx);
bsln            = 2; % prone unweighted default is recording #2. This is used as the reference for all recordings.
dirs            = ls(dataDir); 
dirs            = dirs(3:end,:);
fileNums        = [2, 1, 3, 4, 5];
imgr_idx_list   = [1 2 3 3 + we_segments (3 + we_segments * 2)];
dataFn = {'seq_1','seq_2','seq_3','seq_4','seq_5'};
fig = figure(1);
fig.Units = 'normalized';
fig.OuterPosition = [0 0 1 1];
clf(); 
hold on;
fig.PaperOrientation = 'landscape';

% stats output
TABLE = table();

for i = 1: size(dirs, 1) % for each participant
    if i == 14 || i == 16 % these subjs were noisy
        continue
    end
    stats = struct;
    all_seq_av_breath = struct;
    clf;
    files = load_wr_data( sprintf('%s%s', dataDir, regexprep(dirs(i, :), ' +', '')) );
    [D, imdl_comp] = wr_pp(files, imdl, timeSeriesLength);
    % number of columns in figure
    NCOLS = 2 * we_segments + length(fileNums) - 2;
    for j = fileNums % for each condition
        try
            f = files{j}; 
        catch e
            continue
        end
        if isempty(f)
           continue 
        end
        name = char(remove_underscores(f));
        suffix = ' clean breaths';        
        front = name(1:14); % indices of date and participant number
        back = name(15:end); % indices of experimental phase in file name
        imgr_idx = imgr_idx_list(j);
%--------------------------------------------------------------------------        
        % Data pre-processing and breath selection
        try
            use_data = D.(dataFn{j}).useData;
        catch e
            continue
        end
        FR = D.(dataFn{j}).fs;
        tbv = sum(use_data, 1); % total boundary voltage        
        
        % reconstruct with average frame to find breath landmarks
        imgr = inv_solve(imdl_comp, mean(use_data,2), use_data);
        img_slices = calc_slices(imgr, [inf inf 1]);
        img_slices(isnan(img_slices)) = 0;
        img_slices = img_slices .* lung_segmentation();
        lungZ = squeeze( sum( img_slices, [1 2]) )';
        [end_in, end_ex, Ti_by_Tt, BF] = select_breaths(-lungZ, FR);
%--------------------------------------------------------------------------
        % Plot selected breath boundaries
%         plot_breath_boundaries(name, suffix, -lungZ, end_in, end_ex, front, back, j);
%         saveDir = 'C:\Users\Mark\Documents\GraduateStudies\LAB\EIT-restraint\zzMC\figures\data quality\breath selection\';
%         saveas(gcf, sprintf('%s%s-%.0f_fig.svg', saveDir, front, j));
%         clf();
%--------------------------------------------------------------------------        
        % Image reconstruction
        if strcmp(referenceStyle, 'self')
            glob_ref = mean(use_data(:, end_ex(1, 1): end_ex(2, 1)), 2); % use mean of first breath as reference
        else
            if j==2
                glob_ref = mean(use_data, 2);
            end
        end
        imgr = inv_solve(imdl_comp, glob_ref, use_data); % solve
        
        img_slices = calc_slices(imgr, [inf inf 1]);
        img_slices(isnan(img_slices)) = 0;
        img_slices = img_slices .* lung_segmentation();
%         lungZ = squeeze( sum( img_slices, [1 2]) )';
%         [end_in, end_ex, Ti_by_Tt, BF] = select_breaths(-lungZ, FR);
        
        if j==2
            figure_img = imgr; % compile images from all conditions into a single image
            figure_img.elem_data = zeros(size(imgr.elem_data, 1), NCOLS); % 5 conditions, max we_segments images for 2 conditions, 1 image for 3 conditions.
            figure_img.show_slices.img_cols = NCOLS;
            figure_img.show_slices.sep = sep;
            imgFRC  = figure_img;
            imgTV   = figure_img;
        end % end if
        
        % FRC is average EELI
        imgFRC.elem_data = (imgr.elem_data(:, end_ex(1, :)) + imgr.elem_data(:, end_ex(2, :))) ./ 2;
        % TV is EILI - FRC
        imgTV.elem_data = imgr.elem_data(:, end_in) - imgFRC.elem_data;
        % calculate parameters, store in stats global variable
        calc_img_stats(imgTV, j, 'TV');
        calc_img_stats(imgTV, j, 'fig'); % Add average TV to appropriate column in figure_img
        calc_img_stats(imgFRC, j, 'FRC');
        calc_time_stats(end_in ./ fs, j, 'time');
        calc_time_stats(BF, j, 'BF');
        calc_time_stats(tbv, j, 'av_breath'); % plots
    end % end for j
    
% --------------------------------------------------------------------------
    % Reconstructed images figure
%     cInFig  = 2:12;
    cInFig = 1:13;
    levels = [inf, inf, 1.25; inf, inf, 1; inf, inf, 0.75];
    clim = max(figure_img.elem_data(:, cInFig), [], 'all');
    figure_img.calc_colours.ref_level = 0;
    figure_img.calc_colours.clim = clim;
    final_img = calc_slices(figure_img, levels);
    final_img(isnan(final_img)) = 0;
    
    cntr_of_vnt = calc_cntr_of_vnt(final_img, lung_roi);
    final_img = calc_colours(calc_slices(figure_img, levels));
    cntr_of_vnt = reshape(cntr_of_vnt,2,size(final_img, 3))';

    for j = 1: size(final_img, 3)
        cntr_x = cntr_of_vnt(j, 1);
        cntr_y = cntr_of_vnt(j, 2);
        % COV x coordinate and line
        final_img( (cntr_y - min(4, cntr_y - 1) : min(cntr_y + 4, 32)), cntr_x, j) = 3;
        % COV y coordinate and line
        final_img( cntr_y, :, j ) = 3; 
    end % end for j
    
    final_img2 = mk_mosaic(final_img(6: 27, :, [cInFig, cInFig + NCOLS, cInFig + NCOLS * 2]), [sep, 3], [], length(cInFig));
    final_img2(isnan(final_img2)) = 0;
    if length(cInFig) == 13 
        boundaries = zeros(1, length(files) - 1);
        boundaries(1) = 32 + median(1: sep);
        boundaries(2: end) = boundaries(1) + [1, 1 + we_segments, 1 + we_segments * 2] * (32 + sep);
    elseif length(cInFig)== 11
        boundaries = zeros(1, 2);
        boundaries(1) = 32 + median(1: sep);
        boundaries(2) = boundaries(1) + (we_segments) * (32 + sep);
    end
    % colour background black
    black = find(sum(cmap == [0, 0, 0], 2) == 3);
    final_img2(:, boundaries) = black; % black

%--------------------------------------------------------------------------
    % Normalize stats to unweighted prone
    
    % FRC (mean end ex) as a difference from standing FRC baseline,
    % translated to relative volume by taking it as a ratio of prone TV
    % (mean end ex - mean end ins).
    
    stats.FRCTREL   = (stats.FRCT - stats.FRCT(bsln)) ./ stats.TV(bsln);    
    stats.FRCMREL   = (stats.FRCM - stats.FRCM(bsln)) ./ stats.TV(bsln);
    stats.FRCBREL   = (stats.FRCB - stats.FRCB(bsln)) ./ stats.TV(bsln);
    stats.TVREL     = stats.TV./ stats.TV(bsln);    
    stats.BFREL     = stats.BF./ stats.BF(bsln);    
    stats.MINVNT    = stats.TV .* stats.BF;    
    stats.MINVNTREL = stats.MINVNT ./ stats.MINVNT(bsln);    
    
    for idx = 1: NCOLS
        stats.FRCTRELind{idx} = (stats.FRCTind{idx} - stats.FRCT(bsln)) ./ stats.TV(bsln);
        stats.FRCMRELind{idx} = (stats.FRCMind{idx} - stats.FRCM(bsln)) ./ stats.TV(bsln);
        stats.FRCBRELind{idx} = (stats.FRCBind{idx} - stats.FRCB(bsln)) ./ stats.TV(bsln);
        stats.TVRELind{idx}   = stats.TVind{idx} ./ stats.TV(bsln);
        stats.BFRELind{idx}   = stats.BFind{idx} ./ stats.BF(bsln);
        stats.MINVNTind{idx}  = stats.TVind{idx} .* stats.BFind{idx};
        stats.MINVNTRELind{idx} = stats.MINVNTind{idx} ./ stats.MINVNT(bsln);
        
    end

% --------------------------------------------------------------------------
    % Make Figure
%     paper_figure('reconst', final_img2, i);
%     paper_figure('breath', all_seq_av_breath.figcols, cInFig);
%     paper_figure('stat', stats, cInFig);
% --------------------------------------------------------------------------
    % Export  
%     printPDF(sprintf('%s%s_fig.pdf', saveDir, front));
    header  = fieldnames(stats);
    header1 = {};
    maxLength = 0;
    count1 = 0;
    
    phases = {'U','R','W1','W2','W3','W4','W5','X1','X2','X3','X4','X5','P'};
    % export per-observation parameters
    name = strsplit(front);
    name = name{end};
    
    outParams = {'FRCMRELind', 'MINVNTRELind', 'TVRELind', 'BFRELind'};
    subjTable = mk_subj_table(name, phases, outParams);
    TABLE = [TABLE; subjTable];
end % end for i

writetable(TABLE, feat_save_name);

% -------------------------------------------------------------------------
%                       S U B F U N C T I O N S
% -------------------------------------------------------------------------

function outTable = mk_subj_table(name, phases, outParams)
    global stats
    [age, height, weight, bmi, pos] = get_ahwbp(name);
    outTable = table();
    for i = 1 : length(phases)
        skip = false;
        observations = [];
        phase = phases{i};
        for j = 1 : length(outParams)
            thisParam = outParams{j};
            if isempty(stats.(thisParam){i})
                skip = true;
                break
            else
                observations(:, j) = stats.(thisParam){i};
            end
        end
        if skip
            continue
        end
        
        nObs = size(observations, 1);
        
        nameObs = repmat({name}, nObs, 1);
        ageObs = repmat(age, nObs, 1);
        heightObs = repmat(height, nObs, 1);
        weightObs = repmat(weight, nObs, 1);
        bmiObs = repmat(bmi, nObs, 1);
        posObs = repmat(pos, nObs, 1);
        phaseObs = repmat({phase}, nObs, 1);
        timeObs = stats.timeind{i}';
        T1 = table(nameObs, ageObs, heightObs, weightObs, bmiObs, posObs, phaseObs, timeObs,...
            'VariableNames', {'subj', 'age', 'height', 'weight', 'bmi', 'posture', 'phase', 'time'});
        T2 = array2table(observations, 'VariableNames', outParams);
        thisTable = [T1, T2];
        outTable = [outTable; thisTable];
    end
end

function [age, height, weight, bmi, pos] = get_ahwbp(name)
    ages = [22    41    22    20    38    28    27    21    21    21    19    21    21    21    21    22    44    21    42];
    heights = [170   187   168   180   174   175   187   178   169   172   189   176   172   170   175   183   178   175   169];
    weights = [63.9000   97.3000   65.5000   83.7000   80.1000   84.0000   89.9000   68.2000   60.0000   74.1000   81.4000   68.6000   73.6000   92.3000   64.5000  110.0000  105.5000   89.5000   61.8000];
    bmis = [22.1000   27.8000   23.2000   25.8000   26.5000   27.4000   25.7000   21.5000   21.0000   25.0000   22.8000   22.2000   24.9000   31.9000   21.1000   32.8000   33.3000   29.2000   21.6000];
    postures = [3     3     3     2     2     2     1     1     2     1     1     1     3     3     2     2     3     2     3];
    
    ID = str2double(strrep(name, 'P', ''));
    
    age = ages(ID);
    height = heights(ID);
    weight = weights(ID);
    bmi = bmis(ID);
    pos = postures(ID);
end


function paper_figure(command, data, num)
global cmap;
global we_segments;
B   = 0;
W   = 0.95;
H   = 0.13;     % height of each plot
L   = 0.01;     % left margin of reconstruction and stats plots
BL  = 0.0043;   % left margin of first breath plot
R   = 0.98;
bpb = linspace(BL, W+0.012, length(num)+1);
spc = W/length(num);
shift = 0.04;

if strcmp(command, 'breath')
    mxsz    = 0;
    datamx  = 0;
    datamn  = 0;
    len_data= length(num);
    data_sz = zeros(len_data, 1);
    
    for i = 1:len_data
        idx = num(i);
        try
            if ~isempty(data{idx}) % changed
                data_sz(i) = size(data{idx}, 2);
                mxsz    = max( data_sz(i),  mxsz);
                datamx  = max( max(data{idx},[],'all') ,  datamx);
                datamn  = min( min(data{idx},[],'all') ,  datamn);
            end % end if
        catch
            continue
        end
    end % end for len
    
    start = floor((mxsz- data_sz)./ 2) + 1;
    
    for i = 1:len_data
        idx = num(i);
%     for num = 1:len_data
        ax  = axes('Position',[bpb(i) 0.003+3.5*H spc H],'Box','on');
        s   = start(i);
        e   = start(i) + data_sz(i) - 1;
        hold on;
        
        for j = 1:size(data{idx}, 1) - 1
            plot( (s:e), data{idx}(j,:), 'Color',[0 0 0]+0.75, 'linewidth', 1 ); % grey
        end % end for
        
        plot( (s:e), data{idx}(j+1,:), 'm', 'linewidth', 2 ); % magenta
        ax.XTickLabel=[];ax.YTickLabel=[];ax.XTick=[];ax.YTick=[];
        ax.XLim= ([1 mxsz]);
        ax.YLim= ([datamn datamx]);
    end % end for num

elseif strcmp(command, 'stat')
    LW= 3; 
    xlbls = cell(1,length(data.FRCTREL));
    xlbls{1} = 'U'; xlbls{2} = 'R';
    idx = 3;
    ax1 = axes('Position',[L  B+4.75*H+shift   W   2*H],'Box','on');
    hold on;
    for i = 1:we_segments
        xlbls{idx} = sprintf('W_%i',i);
        xlbls{idx+we_segments} = sprintf('X_%i',i);
        idx=idx+1;
    end
    
    if length(data.FRCTREL) > (2*we_segments+2)
       xlbls{end} = 'P';
       IDX = {[1],[2],[3:idx-1],[idx:length(data.FRCTREL)-1], [length(data.FRCTREL)]};
    else
       IDX = {[1],[2],[3:idx-1],[idx:length(data.FRCTREL)]};
    end
    if length(num)==11
        IDX = {2,[3:idx-1],[idx:length(data.FRCTREL)-1]};
    end
    xlbls = {xlbls{num}};
    for i=1:length(IDX)
        if i == length(IDX)
            e1 = plot( IDX{i}, data.FRCTREL(IDX{i}),  '-o',   'linewidth',LW, 'Color','#A2142F'); 
            e2 = plot( IDX{i}, data.FRCMREL(IDX{i}),  '-o',   'linewidth',LW, 'Color','#D95319');
            e3 = plot( IDX{i}, data.FRCBREL(IDX{i}),  '-o',   'linewidth',LW, 'Color','#0072BD');
            e4 = plot( IDX{i}, data.TVREL(IDX{i}),    '--^',  'linewidth',LW, 'Color','#7E2F8E');
            e5 = plot( IDX{i}, data.BFREL(IDX{i}),    ':s',   'linewidth',LW, 'Color','#77AC30');            
            legEntries = [e1,e2,e3,e4,e5];
            legend(legEntries, '\Delta FRC (top)', '\Delta FRC (middle)', '\Delta FRC (bottom)', 'V_T', '{\it f}_R', 'Location', 'northoutside', 'NumColumns', 5);
        else
            plot( IDX{i}, data.FRCTREL(IDX{i}),  '-o',   'linewidth',LW, 'Color','#A2142F'); 
            plot( IDX{i}, data.FRCMREL(IDX{i}),  '-o',   'linewidth',LW, 'Color','#D95319');
            plot( IDX{i}, data.FRCBREL(IDX{i}),  '-o',   'linewidth',LW, 'Color','#0072BD');
            plot( IDX{i}, data.TVREL(IDX{i}),    '--^',  'linewidth',LW, 'Color','#7E2F8E');
            plot( IDX{i}, data.BFREL(IDX{i}),    ':s',   'linewidth',LW, 'Color','#77AC30');            
        end
    end
    
    ax1.XLim=[num(1)-.5 num(end)+.5];    ax1.XTick=num;
    grid on;
    set(gca,'XTickLabel', xlbls, 'fontsize', 20);
    yrule = ax1.YAxis;
    yrule.FontSize = 10;
    ax1.YAxisLocation= 'right';
    
elseif strcmp(command, 'reconst')
    fg1=axes('Position',[L B+H-shift R H*3]);
    image(data);colorbar;
    axis image
    axis off
    axis equal
    axis tight
    fg1.Colormap= cmap;
%     sgtitle(horzcat('Participant ', num2str(num)));
end % end if

end % end function

% -------------------------------------------------------------------------
function [end_in, end_ex, Ti_by_Tt, BF] = select_breaths(tbv, FR)

    % Breath boundaries
    breaths= find_breaths(tbv);
    end_in= breaths.ins_idx;
    end_ex= breaths.exp_idx;
    exp_to_exp= (diff(end_ex, 1))';
    Te= ((end_ex(2,:)- end_in))';
    Ti= (abs(end_ex(1,:)- end_in))';
    
    reject1 = abs( (end_in - end_ex(1,:)) ./ (end_in - end_ex(2,:)) );
    reject2 = abs( (end_in - end_ex(2,:)) ./ (end_in - end_ex(1,:)) );
    reject3 = (tbv(end_in) - tbv(min(end_ex,[],1))) ./ abs(diff(tbv(end_ex),1));
    
    keep = ((reject1 < 2.5) + (reject2 < 2) + (reject3 > 3)) == 3;
    
    end_in = end_in(keep);
    end_ex = end_ex(:, keep);
    Ti = Ti(keep);
    Te = Te(keep);
    exp_to_exp = exp_to_exp(keep);
    
    Ti_by_Tt= Ti ./ (Ti + Te);
    BF = 60./ (exp_to_exp./ FR); % instantaneous breaths per minute

end % end function

% -------------------------------------------------------------------------

function all_seq_av_breath= calc_average_breath(tbv, end_in, end_ex)

n_breaths   = size(end_ex, 2);
b_lens      = abs(diff(end_ex)); % lengths of breaths
b_centers   = end_in- (end_ex(1,:)) + 1; % index relative to 1 of where inspiration peak occurs
b_start     = repmat(max(b_centers), 1, n_breaths)+ 1- b_centers;
all_seq_av_breath= zeros(n_breaths+1, max(b_start+b_lens));

for i= 1:n_breaths
    b_len   = b_lens(i);
    x1      = end_ex(1,i);
    x2      = end_ex(2,i);
    y       = tbv(x1:x2);
    y       = y-min(y);
    xax     = b_start(i):b_len+b_start(i);
    all_seq_av_breath(i,xax)= y;
end % end for

av_breath   = mean(all_seq_av_breath, 1);
all_seq_av_breath(i+1,:)= av_breath;

end % end function

% -------------------------------------------------------------------------

function cntr_of_vnt= calc_cntr_of_vnt(img, lung_roi)
% img is a 4D array. Output will be 2 x 13 x 3
% remove image components not related to respiration
img_= -(img.*lung_roi); % invert pixel values
% 2 coordinates per image, images/row, number of rows
n_imgs= size(img_, 3);
n_rows= size(img_, 4);
cntr_of_vnt= zeros(2, n_imgs, n_rows); 

for h=1:n_rows
    for i= 1:n_imgs
        this_img= img_(:,:,i,h);
        this_img= (this_img./(sum(this_img,'all'))) *100;
        
        col_sums= sum(this_img, 1);
        c= cumsum(col_sums);
        cn= find(c <=50,1,'last');
        if isempty(cn)
            cntr_of_vnt(1,i,h) = 1;
            cntr_of_vnt(2,i,h) = 1;
            continue
        end % end if
        k= (50- sum(col_sums(1:cn))) / cn;
        cntr_of_vnt(1,i,h)= round( ((cn+k+0.5) / (32+1)) *32); % X coor
        
        row_sums= sum(this_img, 2);
        r= cumsum(row_sums);
        rn= find(r <=50,1,'last');
        k= (50- sum(row_sums(1:rn))) / rn;
        cntr_of_vnt(2,i,h)= round( ((rn+k+0.5) / (32+1)) *32); % Y coor
        
        if isnan(cntr_of_vnt(1,i,h)) || isnan(cntr_of_vnt(2,i,h))
            keyboard;
        end % end if
    end % end for i
end % end for h
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

param_list= {'BF', 'FIT', 'av_breath', 'time'};
if ~contains(param_list, param)
    disp('Unrecognized statistical parameter. Accepted parameters are: BF, FIT, av_breath');
    return
end % end if

if strcmp(param, 'av_breath')
    do_breaths = true;
else
    do_breaths  = false;
    paramSD     = horzcat(param, 'SD');
    paramN      = horzcat(param, 'N');
    paramI      = horzcat(param, 'ind');
    if ~isfield(stats, param)
        stats.(param) = nan(max(imgr_idx_list), 1);
    end
    if ~isfield(stats, paramSD)
        stats.(paramSD) = nan(max(imgr_idx_list), 1);
    end
    if ~isfield(stats, paramN)
        stats.(paramN) = nan(max(imgr_idx_list), 1);
    end
    if ~isfield(stats, paramI)
        stats.(paramI) = cell(max(imgr_idx_list), 1);
    end
end

stat_idx= imgr_idx_list(j);

if j<3 || j==5
    if do_breaths
        all_seq_av_breath.figcols{stat_idx} = calc_average_breath(data, end_in, end_ex);
    else
        stats.(param)(stat_idx,1)   = mean(data);
        stats.(paramSD)(stat_idx,1) = std(data);
        stats.(paramN)(stat_idx,1)  = length(data);
        stats.(paramI){stat_idx}    = data;
    end
elseif j==3 || j==4
    for m= 1:length(av_frm_s_inc)
        begins_after    = end_ex(1,:) > av_frm_s_inc(m);
        ends_before_at  = end_ex(2,:) <= (av_frm_s_inc(m)+ n_frame_s);
        img_idx         = find((begins_after + ends_before_at)==2);
        if ~isempty(img_idx)
            if do_breaths
                all_seq_av_breath.figcols{stat_idx}= calc_average_breath(data, end_in(img_idx), end_ex(:,img_idx));
            else
                stats.(param)(stat_idx,1)   = mean(data(img_idx));
                stats.(paramSD)(stat_idx,1) = std(data(img_idx));
                stats.(paramN)(stat_idx,1)  = length(data(img_idx));
                stats.(paramI){stat_idx}    = data(img_idx);
            end
        end % end if
        stat_idx = stat_idx + 1;
    end % end for m
else
    disp("Unrecognized experiment phase.")
end % end if

end

% -------------------------------------------------------------------------
function calc_img_stats(img, j, param)

param_list= {'FRC', 'TV', 'fig'};
if ~contains(param_list, param)
    disp('Unrecognized statistical parameter. Accepted parameters are: FRC, TV, fig');
    return
end % end if

if strcmp(param, 'FRC')
    param_list  = {'FRCT', 'FRCM', 'FRCB'};
    levels      = [inf,inf,1.25; inf,inf,1; inf,inf,0.75];
elseif strcmp(param, 'TV')
    param_list  = {param};
    levels      = [inf,inf,1];
else
    do_img_stat_calc(param, img, [], j);
    return;
end % end if

for i=1:length(param_list)
    param       = param_list{i};
    img_slices  = calc_slices(img, levels(i,:));
    img_slices(isnan(img_slices))= 0;
    img_slices  = img_slices .* lung_segmentation();
    img_slices  = -img_slices;
    img_sums    = squeeze(sum(img_slices, [1 2]));
    do_img_stat_calc(param, img, img_sums, j);
end

end % end function

function do_img_stat_calc(param, img, img_sums, j)

global stats;
global end_ex;
global imgr_idx_list;
global av_frm_s_inc;
global n_frame_s;
global figure_img;

if ~strcmp(param, 'fig') && ~strcmp(param, 'av_breath')
    paramSD = horzcat(param, 'SD');
    paramN  = horzcat(param, 'N');
    paramI  = horzcat(param, 'ind');
    if ~isfield(stats, param)
        stats.(param) = nan(max(imgr_idx_list), 1);
    end
    if ~isfield(stats, paramSD)
        stats.(paramSD) = nan(max(imgr_idx_list), 1);
    end
    if ~isfield(stats, paramN)
        stats.(paramN) = nan(max(imgr_idx_list), 1);
    end
    if ~isfield(stats, paramI)
        stats.(paramI) = cell(max(imgr_idx_list), 1);
    end
end % end if
stat_idx = imgr_idx_list(j);

if j<3 || j==5
    if strcmp(param, 'fig')
        figure_img.elem_data(:,stat_idx) = mean(img.elem_data, 2);
    else
        stats.(param)(stat_idx,1)   = mean(img_sums);    
        stats.(paramSD)(stat_idx,1) = std(img_sums); 
        stats.(paramN)(stat_idx,1)  = length(img_sums);
        stats.(paramI){stat_idx}    = img_sums;
    end % end if
    
elseif j==3 || j==4
    for m = 1:length(av_frm_s_inc)
        begins_after = end_ex(1,:) > av_frm_s_inc(m);
        ends_before_at = end_ex(2,:) <= (av_frm_s_inc(m)+ n_frame_s);
        img_idx = find((begins_after + ends_before_at)==2);
        if ~isempty(img_idx)
            if strcmp(param, 'fig')
                figure_img.elem_data(:,stat_idx) = mean(img.elem_data(:,img_idx), 2);
            else
                stats.(param)(stat_idx,1)   = mean(img_sums(img_idx));
                stats.(paramSD)(stat_idx,1) = std(img_sums(img_idx)); 
                stats.(paramN)(stat_idx,1)  = length(img_sums(img_idx));
                stats.(paramI){stat_idx}    = img_sums(img_idx);
            end % end if
        end % end if
        stat_idx = stat_idx + 1;
    end % end for m
else
    disp("Unrecognized experimental phase.")
end % end if

end


function cmap = config_cmap()
global cmap;
calc_colours('cmap_type', 'blue_black_red');
colormap(calc_colours('colourmap'));
cmap=colormap*2;
cmap(cmap>1)=1;
black= find(sum(cmap==[0,0,0],2)==3);

white_grad= linspace(0,1,black-1)';
cmap(black:end,1)=white_grad;
cmap(black:end,2)=white_grad;

white_grad= linspace(0,1,black)';
white_grad=flipud(white_grad);
cmap(1:black,2)=white_grad;
cmap(1:black,3)=white_grad;

cmap= [[1, 1, 1]; flipud(cmap(2:end, :))];

end % end function

function plot_breath_boundaries(name, suffix, tbv, end_in, end_ex, front, back, j)

n_breaths= length(end_in);
fg1 = figure(1);clf;
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
end % end function
