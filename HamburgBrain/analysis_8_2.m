% -------------------------------------------------------------------------
% Description: This script is used for analyzing and visualizing swine
% stroke data.
% -------------------------------------------------------------------------
% Author:
%   Mark Campbell
%   Carleton University
%   markacampbell@cmail.carleton.ca
% -------------------------------------------------------------------------
global seg;
run 'myStartup.m';
% pigs= {'8-2','9-2','10-2'};
pigs= {'10-2'};
for q = 1:length(pigs)
    close all
    bigFig();
    bigFig();
    pig = pigs{q};
    [fmdl, imdl] = get_pig_mdl(pig);
    % Load data
    ref = 'self';
    % ref = 'baseline';
    
    D = load_HamburgBrain_data(pig, ref);
    seg = get_brain_segmentation_for_pig(pig);
    fn = fieldnames(D);
    suffix = date;

    switch pig
        case '8-2'   
            opt.period = [2,2,2,1];
            sel = [2,4,3];
            temp = D.seq1.eit.apn; D.seq1.eit.apn = 1868; 
            temp2 = D.seq2.eit.apn; D.seq2.eit.apn = 1068;
        case '9-2'     
            opt.period = [0,1,1,1];
%             sel = [1,3,4];
            sel = [1,3];
        case '10-2'    
            opt.period = [1,1,1,1,1,1];
%             sel = [1,3,5];
            sel = [2,4,6];
        case '11-2'    
            opt.period = [1,2,0,2,1,1];
            sel = [1,3,5];
            temp = D.seq3.eit.apn; D.seq3.eit.apn = 843;
        case '12-2'    
            opt.period = [2,0,0,0,2,0];
            sel = [1,3,5];
            temp = D.seq3.eit.apn; D.seq3.eit.apn = 939; % 1223;
    end % end switch
    
% -------------------------------------------------------------------------
    % 1. Injection Images Figure
    cd('C:\Users\Mark\Documents\GraduateStudies\LAB\HamburgBrain\Figures\paper\inj');

    switch pig
%         case '8-2';     ts = [5,5,5,5];       te = [20,20,20,20]; % use seq2 as ref
%         case '9-2';     ts = [5,5,5,5];       te = [20,20,20,20];
%         case '10-2';    ts = [5,5,5,5,5,5];   te = [20,20,20,20,20,20];
%         case '11-2';    ts = [5,5,5,5,5,5];   te = [20,20,20,20,20,20];
%         case '12-2';    ts = [5,5,5,5,5,5];   te = [20,20,20,20,20,20];
        case '8-2';     ts = [1,1,1,1];       te = [20,20,20,20]; % use seq2 as ref
        case '9-2';     ts = [1,1,1,1];       te = [20,20,20,20];
        case '10-2';    ts = [1,1,1,1,1,1];   te = [20,20,20,20,20,20];
        case '11-2';    ts = [1,1,1,1,1,1];   te = [20,20,20,20,20,20];
        case '12-2';    ts = [1,1,1,1,1,1];   te = [20,20,20,20,20,20];
    end % end switch

    % options
    nFrames = length(ts(1):te(1));
    % make figure
    show_inj_fig(D, ts, te, nFrames, sel);
    % save figures
    figure(1); pause(0.5);
    printPDF( sprintf('%sInjectionFigure', char(D.(fn{1}).pig)) );
    figure(2);
    printPDF( sprintf('%sInjectionFigurePlots', char(D.(fn{1}).pig)) );   
    pause(1);

% -------------------------------------------------------------------------
%     % 2. Total change over cardiac cycle
%     clf();
%     cd('C:\Users\Mark\Documents\GraduateStudies\LAB\HamburgBrain\Figures\paper\tcoc');
%     hold on;
%     % otptions
%     name = horzcat(char(pig), ' ',ref,' reference. Total change over average cardiac cycle');
%     title( name );
%     % make figure
%     delta_heatmap(D, sel, opt);
%     colorbar();
%     % title(horzcat(num2str(sel), ' - ',D.(fn{sel}).name));
%     colormap jet;
%     fig = gcf;
%     fig.Colormap(1,:) = [1 1 1] * 0.75;
%     printPDF( horzcat(char(pig),'TCOCC') );
%     pause(1);

% -------------------------------------------------------------------------
    % 3. Compare average CC for all sequences with adjusted clim
    figure(1); clf();
    figure(2); clf();
    cd('C:\Users\Mark\Documents\GraduateStudies\LAB\HamburgBrain\Figures\paper\ensemble');
    % options
    opt.sidelen = inf;
    opt.align = 2;
    % opt.sidelen = 5;
    switch pig
        case '8-2'
            TITLE = '8-2 Ensemble Average of Pulsatile Signal. baseline (top), after perfusion reduction (middle), 30 minutes after perfusion reduction (bottom)';
        case '9-2'
%             TITLE = '9-2 Ensemble Average of Pulsatile Signal. baseline (top), after perfusion reduction (middle), 60 minutes after perfusion reduction (bottom)';
            TITLE = '9-2 Ensemble Average of Pulsatile Signal. baseline (top), after perfusion reduction (bottom)';
        case '10-2'
            TITLE = '10-2 CVC Ensemble Average of Pulsatile Signal. baseline (top), after stroke induction (middle), 6 hours after stroke induction (bottom)';
        case '11-2'
            TITLE = '11-2 CVC Ensemble Average of Pulsatile Signal. baseline (top), after stroke induction (middle), 4 hours after stroke induction (bottom)';
        case '12-2'
            TITLE = '12-2 CVC Ensemble Average of Pulsatile Signal. baseline (top), after stroke induction (middle), 3.5 hours after stroke induction (bottom)';
    end % end switch
    % make figures
    compare_pre_inj_imgs(D, sel, opt);
    % save figures
    figure(1); 
    title(TITLE); colorbar(); pause(0.5);
    printPDF(horzcat(char(pig), 'ensemble'));
    figure(2);
    title(TITLE);
    printPDF(horzcat(char(pig), 'ensemblePlots'));
% -------------------------------------------------------------------------
    % revert landmarks
    switch pig
        case '8-2';     D.seq.eit.apn = temp; D.seq2.eit.apn = temp2;
        case '11-2';    D.seq3.eit.apn = temp;
        case '12-2';    D.seq3.eit.apn = temp;
    end % end switch
end % end for

% ======================================================================= %

function temp = adjust_for_poi(seq, start, stop)

    shift = start - 1;
    useLndmrk = ( (seq.eit.vals(1,:) >= start) + (seq.eit.vals(2,:) <= stop) ) == 2;
    temp = seq;
    temp.imgr.elem_data = seq.imgr.elem_data(:, start:stop);
    temp.eit.peaks = seq.eit.peaks( useLndmrk ) - shift;
    temp.eit.vals = seq.eit.vals( :, useLndmrk ) - shift ;

end

% ======================================================================= %

function show_inj_imgs(seq, opt)
    SECONDS = 15;
    if nargin == 1
        opt = struct;
    end % end if
    
    opt.show = true;
    
    fs = round(seq.eit.fs);
    % take 2 seconds after inj to be sure synchronization flush has been
    % excluded.
    start = seq.eit.inj + 2 * fs;
    % take first 10 seconds after inj
    if isfield(opt, 'frames')
        start = start + frames - 1;
    else
        stop = start + SECONDS * fs - 1;
    end % end if
    
    temp = adjust_for_poi(seq, start, stop);
    get_mean_cc(temp, opt);
    title(horzcat(seq.name, ' Post-injection conductivity from individual and ensemble-averaged (bottom) Cardiac Cycles'));
    
end % end function

% ======================================================================= %

function show_pre_inj_img(seq, opt)
    
    if seq.eit.apn < 1
        start = 1;
    else
%         start = seq.eit.apn;
        start = 1;
    end % end if
    
    stop = seq.eit.inj;
    temp = adjust_for_poi(seq, start, stop);
    opt.show = true;
    get_mean_cc( temp, opt );
    title(horzcat(seq.name, ' Pre-injection conductivity from individual and ensemble-averaged (bottom) Cardiac Cycles'));

end % end function

% ======================================================================= %

function meanCycle = get_mean_cc(seq, opt)
    
    if nargin == 1
       opt.align = 1;
       opt.av = 'all';
    end % end if
    
    if isfield(opt, 'align')
        align = opt.align;
    else
        align = 2;
    end % end if
    
    if isfield(opt, 'av')
        av = opt.av;
    else
        av = 'all';
    end % end if
    
    if ~isfield(opt, 'show')
        opt.show = false;
    end % end if
    
    if isfield(opt, 'sidelen')
        SIDELEN = opt.sidelen;
    else
        SIDELEN = inf;
    end % end if
    
    alignDiast1 = false;
    alignSyst = false;
    alignDiast2 = false;
    
    if align == 1
        alignDiast1 = true;
    elseif align == 2
        alignSyst = true;
    elseif align == 3
        alignDiast2 = true;
    end

    use_syst = seq.eit.peaks;
    use_diast = seq.eit.vals;
    nCycles = length(use_syst);
    
    % set up matrix to view each CC as one row in show_sliecs.
    ccLen = diff(use_diast) + 1;
    z_diast1 =  seq.imgr.elem_data(:, use_diast(2,:));
    z_diast2 =  seq.imgr.elem_data(:, use_diast(1,:));
    z_diast = (z_diast1 + z_diast2) ./ 2;
    ls = use_syst - use_diast(1,:); % how many frames from left side systole occurs
    
    if alignSyst % show 5 frames before and 5 frames after systole
        lpad = abs( ls - max(ls) );
        systFrm = lpad(1) + ls(1);
    elseif alignDiast1
        lpad = zeros(nCycles, 1);
    elseif alignDiast2
        lpad = abs( ccLen - max(ccLen) );
    end % end if
    
    longestCC = max(lpad + ccLen);
    SIDELEN1 = max(1, systFrm - SIDELEN);
    SIDELEN2 = min(systFrm + SIDELEN, longestCC);
    showRng = SIDELEN1 : SIDELEN2;
    
    z_frames = zeros( size(seq.imgr.elem_data, 1), longestCC, nCycles);
    meanCycle = z_frames(:, 1:longestCC);
    nContributed = zeros(1, longestCC);
    % create full image set for each CC that is the difference between each
    % frame and the mean diastole for that CC.
    for i = 1:nCycles
        startIdx = lpad(i) + 1;
        endIdx = startIdx + ccLen(i) - 1;
        imgs = seq.imgr.elem_data( :, use_diast(1,i): use_diast(2,i) ) - z_diast(:,i);
        z_frames(:, startIdx:endIdx , i) = imgs;
        meanCycle(:, startIdx: endIdx) = meanCycle(:, startIdx: endIdx ) + imgs;
        nContributed( startIdx: endIdx) = nContributed( startIdx: endIdx) + 1;
    end % end for
    
    longestCC = length(showRng);
    meanCycle = meanCycle(:, showRng) ./ nContributed(showRng);
    z_frames = z_frames(:, showRng, :);
    z_frames = reshape( z_frames, size(z_frames,1), longestCC * nCycles );
    meanCycle( isnan(meanCycle) ) = 0;
    
    if opt.show
        if strcmp(opt.av, 'all')
            z_frames = [z_frames, meanCycle];
        elseif strcmp(opt.av, 'mean')
            z_frames = meanCycle;
        end % end if
        puls_img = seq.imgr;
        puls_img.elem_data = z_frames;
        puls_img.calc_colours.clim = max(meanCycle, [], 'all');
        puls_img.show_slices.img_cols = longestCC;
        show_slices(puls_img);
    end % end if
end % end function

% ======================================================================= %

function compare_post_vnt_imgs(D, sel, opt)

    fn = fieldnames(D);
    
    if nargin == 1
        sel = 1:length(fn);
    elseif isempty(sel)
        sel = 1:length(fn);
    end % end if
    
    MC = struct;
    maxCcLen = 0;
    CcLens = zeros(1, length(sel));
    for j = 1:length(sel)
        i = sel(j);
        seq = D.(fn{i});
        start = seq.eit.vnt;        
        stop = size(seq.imgr.elem_data, 2);
        temp = adjust_for_poi(seq, start, stop);
        MC.(fn{i}) = get_mean_cc( temp, opt );
        CcLens(j) = size(MC.(fn{i}), 2);
        maxCcLen = max( maxCcLen, CcLens(j) );
    end % end for
    
    fn = fieldnames(MC);
    nElems = size(MC.(fn{j}), 1);
    addFrames = abs(CcLens - maxCcLen);
    outImg = D.(fn{j}).imgr;
    outImg.elem_data = [];
    outImg.show_slices.img_cols = maxCcLen;
    
    for i = 1:length(sel)
        if addFrames(i) > 0
            paddedAvCc = [ zeros(nElems, addFrames(i)), MC.(fn{i}) ];
        else
            paddedAvCc = MC.(fn{i});
        end % end if
        outImg.elem_data = [outImg.elem_data, paddedAvCc];
    end % end for
    
    outImg.calc_colours.clim = max(outImg.elem_data, [], 'all');
    show_slices(outImg);
    
end % end function

% ======================================================================= %

function [start,stop] = define_period(seq, period)
    injOpt = [];
    if length(period) > 1
        injOpt = period(2:end);
        period = period(1);
    end % end if

    switch period
        case 0 % start to apn
            start = 1;
            stop = max(1, seq.eit.apn);
        case 1 % start to inj
            start = 1;
            stop = seq.eit.inj;
        case 2 % apn to inj
            if seq.eit.apn <= 0
                start = 1;
            else
                start = seq.eit.apn;
            end % end if
            stop = seq.eit.inj;
        case 3 % inj + 10 seconds
            if length(injOpt) == 1
                startExtend = 0;
                stopExtend = injOpt(1);
            else
                startExtend = injOpt(1);
                stopExtend = injOpt(2);
            end % end if
            start = seq.eit.inj + round( seq.eit.fs * startExtend );
            stop = seq.eit.inj + round( seq.eit.fs * stopExtend );

        case 4 % inj to vnt
            start = seq.eit.inj;
            stop = seq.eit.vnt;
        case 5 % vnt to end
            start = seq.eit.vnt;
            stop = size(seq.imgr.elem_data, 2);

    end % end switch
end % end function

% ======================================================================= %

function delta_heatmap(D, sel, opt)

    fn = fieldnames(D);
    
    if nargin == 1
        sel = 1:length(fn);
        opt = struct;
        opt.period = ones(length(fn),1);
    elseif isempty(sel)
        sel = 1:length(fn);
    end % end if
    
    opt.sidelen = inf;
    heatmapFrames = zeros( size(D.seq1.imgr.elem_data,1), length(sel) );
    
    for j = 1:length(sel)
        i = sel(j);
        seq = D.(fn{i});
        [start,stop] = define_period(seq, opt.period(i));
        temp = adjust_for_poi(seq, start, stop);
        meanCC = get_mean_cc( temp, opt );
        heatmapFrames(:,j) = max(meanCC, [], 2) - min(meanCC, [], 2);
        % normalize
        heatmapFrames(:,j) = ( heatmapFrames(:,j) - min(heatmapFrames(:,j)) ) ./ ( max(heatmapFrames(:,j)) - min(heatmapFrames(:,j)) );
    end % end for
    
    clim = max(heatmapFrames, [], 'all'); % adjust clim
    temp = D.(fn{1}).imgr;
    temp.calc_colours.clim = clim;
    temp.elem_data = heatmapFrames;
    temp.show_slices.img_cols = length(sel);
    bigFig();
	show_slices(temp);
%     assert(min(img(img~=1)) >= 128, 'Min error');
%     img = img - 127;
%     img(img < 0) = 0;
%     image(img*2);

end % end function

% ======================================================================= %

function compare_pre_inj_imgs(D, sel, opt)
    IMGCOLS = 10;
    fn = fieldnames(D);
    period = opt.period;
    if nargin == 1
        sel = 1:length(fn);
    elseif isempty(sel)
        sel = 1:length(fn);
    end % end if
    
    MC = struct;
    maxCcLen = 0;
    CcLens = zeros(1, length(sel));
    for j = 1:length(sel)
        i = sel(j);
        seq = D.(fn{i});
        [start,stop] = define_period(seq, period(i));        
        temp = adjust_for_poi(seq, start, stop);
        MC.(fn{i}) = get_mean_cc( temp, opt );
        CcLens(j) = size(MC.(fn{i}), 2);
        maxCcLen = max( maxCcLen, CcLens(j) );
    end % end for
    
    fn = fieldnames(MC);
    nElems = size(MC.(fn{j}), 1);
    if mod(maxCcLen,IMGCOLS) <= 5
        maxCcLen = maxCcLen - mod(maxCcLen,IMGCOLS);
    else
        maxCcLen = maxCcLen + IMGCOLS - mod(maxCcLen,IMGCOLS);
    end % end if
    
    addFrames = maxCcLen - CcLens;
    outImg = D.(fn{j}).imgr;
    outImg.elem_data = [];
    outImg.show_slices.img_cols = IMGCOLS;
    
    for i = 1:length(sel)
        if addFrames(i) > 0
%             paddedAvCc = [ nan(nElems, addFrames(i)), MC.(fn{i}), nan(nElems, mod(maxCcLen,IMGCOLS)) ];
            paddedAvCc = [ MC.(fn{i}), nan(nElems, addFrames(i)) ];
        else
            paddedAvCc = MC.(fn{i});
            paddedAvCc = paddedAvCc(:,1:maxCcLen);
        end % end if
        outImg.elem_data = [outImg.elem_data, paddedAvCc];
    end % end for
    
    show_slice_and_brain_seg_z(outImg, fn);
end % end function

% ======================================================================= %

function cmap = confg_cmap()

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

end % end funcion

% ======================================================================= %

function bigFig()
    figure('units','normalized','outerposition',[0 0 1 1]);
end % end function

% ======================================================================= %

function show_inj(seq, ts, te)
    [start,stop] = define_period(seq, [3, ts, te]);
    temp = seq.imgr;
    temp.elem_data = seq.imgr.elem_data(:, start:stop);
    temp.calc_colours.clim = max(temp.elem_data, [] ,'all');
    temp.show_slices.img_cols = round(seq.eit.fs);
    show_slices(temp);
end % end function

% ======================================================================= %

function show_inj_fig(D, ts, te, nFrames, sel)
    IMGCOLS = 10;
    fn = fieldnames(D);
    fn2 = cell(length(sel),1);
    ts = ts(sel);
    te = te(sel);    
    outImg = D.seq1.imgr;
    injFrames = zeros( size(outImg.elem_data, 1), nFrames, length(sel) );

    for i=1:length(sel)
        seq = D.(fn{sel(i)});
        [start,stop] = define_period(seq, [3, ts(i), te(i)]);
        injFrameIdx = round( linspace(start, stop, nFrames+1) );        
        injFrameIdx = injFrameIdx(1:nFrames); % each frame is the start of a second this way.
        injFrames(:,:,i) = seq.imgr.elem_data(:, injFrameIdx);
        fn2{i} = fn{sel(i)};
        
        % want to increase time resolution of line plots for brain
        % conductivity
        injFrameIdx2 = start:stop;
        if i == 1
            injFrames2 = zeros( size(outImg.elem_data, 1), length(injFrameIdx2), length(sel) );
        end % end if
        injFrames2(:,:,i) = seq.imgr.elem_data(:, injFrameIdx2);
    end % end for
    
    outImg.elem_data = reshape(injFrames, size(outImg.elem_data, 1), nFrames * i);
    outImg.show_slices.img_cols = IMGCOLS;
    
    % new stuff with higher temporal res for slices in inj
    plotimg = D.seq1.imgr;
    plotimg.elem_data = reshape(injFrames2, size(plotimg.elem_data, 1), length(injFrameIdx2) * i); % num elems x (n seq x length each series)
    plotimg.show_slices.img_cols = length(injFrameIdx2);  
    
    show_slice_and_brain_seg_z(outImg, fn2, plotimg);
    figure(1);
    title( sprintf('Subject %s - Reconstructed images from %i - %i seconds after saline bolus injection', char(D.(fn{i}).pig), ts(i), te(i)) );
    figure(2);    
    title( sprintf('Subject %s - Reconstructed images from %i - %i seconds after saline bolus injection', char(D.(fn{i}).pig), ts(i), te(i)) );
end % end function

% ======================================================================= %

function show_slice_and_brain_seg_z(img, fn, plotimg)
    global seg;
    
    % show slices
    figure(1);
    binSeg = (seg > 0)*5; % binary segmentation (all) 5x multiplier for brain in show slices.
    binSeg(binSeg == 0) = 1;
    slc = calc_slices(img).* binSeg;
    img.calc_colours.clim = max(slc, [] ,'all');
    rowPerSeq = size(img.elem_data, 2) / length(fn) / img.show_slices.img_cols;
    fivexSlc = calc_colours(slc, img);
    imagesc(mk_mosaic(fivexSlc,0,[],img.show_slices.img_cols)); axis equal; axis image
%     show_slices(img); axis on
    ax = gca;
    ax.YTick = 32;
    for i=2:length(fn)
        ax.YTick = [ax.YTick, ax.YTick(i-1) + rowPerSeq*64];
    end % end for
    ax.YTickLabel = cell(length(fn), 1);
    for i = 1:length(fn)
        ax.YTickLabel{i} = sprintf('Sequence %i', i);
    end % end for
    ax.XTickLabel = {};
    
    % show line plots of brain conductivity
    if nargin == 3 % if we want the higher temporal res plots
        img = plotimg;
        nFrames = img.show_slices.img_cols;
    else
        nFrames = img.show_slices.img_cols * rowPerSeq;
    end
    
    img.calc_colours.clim = max(img.elem_data, [] ,'all');
    slc = calc_slices(img);
    seqcount = 0;
    ymax = -inf;
    ymin = inf;
    figure(2);
    w = 10;      % number of subplot columns
    z = w-1;    % number of subplot columns for line plot
    for i = 1:length(fn)
        seqcount    = seqcount + 1;
        stop        = nFrames * i;
        start       = stop + 1 - nFrames;
        substop     = i * w;
        substart    = substop - z;
        subrng      = substart: substop-1;

        % plot brain segmentation conductivity
        subplot(length(fn), w, substop);
        show_seg(); % show brains segmentation as legend        
        subplot(length(fn), w, subrng);
        plot_brain_seg(slc(:, :, start: stop));
        fg = gcf;
        fg.Children(3).XGrid = 'on';
        if fg.Children(3).YLim(1) < ymin
            ymin = fg.Children(3).YLim(1);
        end % end if
        if fg.Children(3).YLim(2) > ymax
            ymax = fg.Children(3).YLim(2);
        end % end if
        
    end % end for

    for i=3:3:length(fn)*3
        fg.Children(i).YLim = [ymin*1.05 ymax*1.05];
    end

end % end function

% ======================================================================= %

function show_seg()
    global seg;
    trimseg = seg;
    rowsum = sum(trimseg, 2);
    colsum = sum(trimseg, 1);
    trimseg = trimseg(rowsum>0, colsum>0);
    imagesc(trimseg);
    fg = gcf();
    fg.Colormap = [ 1 1 1;...
                    0.8500 0.3250 0.0980;...
                    0.9290 0.6940 0.1250;...
                    0.4940 0.1840 0.5560;...
                    0.4660 0.6740 0.1880;...
                    0.3010 0.7450 0.9330;...
                    0.6350 0.0780 0.1840];
    axis equal; 
    axis off;
end % end function


function plot_brain_seg(slc)
    global seg;
    binSeg = (seg > 0)*1; % binary segmentation (all)
    slc(isnan(slc)) = 0;
    nPix = sum(binSeg,'all');
    sumBrain = squeeze( sum( slc .* binSeg, [1,2] ) );
    % plot global signal first, dark blue colour.
    plot(sumBrain./nPix, 'linewidth', 2);
    legend();
    xlim([0.5, size(slc,3)+0.5]);
    hold on;
    % plot each of the 6 regions next.
    for j = 1:max(seg,[],'all')
        mask = (seg == j) *1;
        nPix = sum(mask,'all');
        plot( squeeze( sum( slc .* mask, [1,2] ) )./nPix, 'linewidth', 2);  
        xlim([0.5, size(slc,3)+0.5]);
    end % end for
    
    ax = gca;
    if size(slc,3) > 75
        ax.XTick = linspace(1,size(slc,3),20);
        ax.XTickLabel = 1:20;
    else
        ax.XTick = 1:size(slc,3);
        
    end % end if
    
    % set axis labels
    fg = gcf;
    leg = fg.Children(1);
    leg.String = {'Whole Brain Average'};
    ax = axes(fg,'visible','off');
    ax.Title.Visible='on';
    ax.XLabel.Visible='on';
    ax.YLabel.Visible='on';
    ax.XGrid = 'on';
    ylabel(ax,'Average conductivity per pixel in segmentation', 'FontSize', 20);
    if size(slc,3) > 75
        xlabel('Time (s)', 'FontSize', 20);
    else
        xlabel('Frame of cardiac cycle', 'FontSize', 20);
    end % end if    

end % end function


function printPDF(filename)
    h = gcf;
    set(h, 'PaperUnits','centimeters');
    set(h, 'Units','centimeters');
    pos=get(h,'Position');
    set(h, 'PaperSize', [pos(3) pos(4)]);
    set(h, 'PaperPositionMode', 'manual');
    set(h, 'PaperPosition',[0 0 pos(3) pos(4)]);
    print('-dpdf',filename);
end % end function
