% -------------------------------------------------------------------------
% Description: This script is used for analyzing and visualizing swine
% stroke data.
% -------------------------------------------------------------------------
% Author:
%   Mark Campbell
%   Carleton University
%   markacampbell@cmail.carleton.ca
% NOTES:
%   NORMALIZE:
%       true:
%           normalize colours across sequences
%       false:
%           normalize colours within each sequence
%   ref:
%       self: 
%           use mean frame between injection and ventilation landmarks of
%           this sequence to reconstruct.
%       baseline: 
%           use mean frame between injection and ventilation landmarks of
%           baseline sequenceto reconstruct.
% -------------------------------------------------------------------------
global seg;
run 'myStartup.m';
pigs = {'8-2', '9-2', '10-2', '11-2', '12-2'};
pigs = {'12-2'};
do = 2.5;
modelDir = 'E:/University/Masters/HamburgBrain/Models';
figDir = 'E:/University/Masters/HamburgBrain/Figures/NF0_5';
for q = 1: length(pigs)
    close all
    bigFig();
    bigFig();
    pig = pigs{q};
    [fmdl, imdl] = get_pig_mdl(pig, modelDir);
    % Load data
%     ref = 'self';
    ref = 'baseline';
    
    D = load_HamburgBrain_data(pig, ref, modelDir);
    seg = get_brain_segmentation_for_pig(pig);
    fn = fieldnames(D);
    suffix = date;

    switch pig
        case '8-2'   
            opt.period = [2, 2, 2, 1];
            sel = [2, 4, 3];
            % baseline 1
            D.seq1.eit.apn = 1868;
            % baseline 2
            D.seq2.eit.apn = 1068;
            D.seq2.eit.inj = 1500;
            D.seq2.eit.vnt = 2475;
            % post-infarct
            D.seq3.eit.vnt = 2592;
            % infarct
            D.seq4.eit.inj = 2124;
            D.seq4.eit.vnt = 2916;
            
        case '9-2'     
            opt.period = [0, 1, 1, 1];
            sel = [1, 3];
            % baseline
            D.seq1.eit.inj = 1274;
            D.seq1.eit.vnt = 3149;
            % infarct
            D.seq3.eit.vnt = 2222;
            
        case '10-2'    
            opt.period = [1,1,1,1,1,1];
%             sel = [1,3,5];
            sel = [2, 4, 6];
            % baseline
            D.seq2.eit.inj = 744;
            D.seq2.eit.vnt = 1694;
            % infarct
            D.seq4.eit.vnt = 1755;
            
            D.seq6.eit.vnt = 1762;
            
        case '11-2'    
            opt.period = [1,2,0,2,1,1];
            sel = [1, 3, 5];
            % baseline
            D.seq1.eit.vnt = 1627;
            % infarct
            D.seq3.eit.apn = 843;
            D.seq3.eit.vnt = 3319;
            % post-infarct
            D.seq5.eit.vnt = 2151;
            
        case '12-2'    
            opt.period = [2,0,0,0,2,0];
            sel = [1, 3, 5];
            
            D.seq1.eit.inj = 1045;
            D.seq1.eit.vnt = 1682;
            % infarct
            D.seq3.eit.apn = 939; % 1223;
            D.seq3.eit.inj = 2795;
            D.seq3.eit.vnt = 4521;
            % post-infarct
            D.seq5.eit.inj = 2120; % 2704 at flat
            D.seq5.eit.vnt = 4601;
                        
    end % end switch
    clf();
%     seq = D.seq2;
%     seq = D.seq3;
%     seq = D.seq4;
%     plot(sum(seq.eit.fdata,1));
%     hold on;
%     xline(seq.eit.inj);
%     xline(seq.eit.vnt);
%     keyboard;
% -------------------------------------------------------------------------

    % 1. Injection Images Figure
    if do == 1
%         cd('E:\University\Masters\HamburgBrain\Docs\paper\figures\NF2\inj');
        switch pig
            case '8-2';     ts = [1,1,1,1];       te = [20,20,20,20]; % use seq2 as ref
            case '9-2';     ts = [1,1,1,1];       te = [20,20,20,20];
            case '10-2';    ts = [1,1,1,1,1,1];   te = [20,20,20,20,20,20];
            case '11-2';    ts = [1,1,1,1,1,1];   te = [20,20,20,20,20,20];
            case '12-2';    ts = [1,1,1,1,1,1];   te = [20,20,20,20,20,20];
        end % end switch

        % options
        opt.normalize = true;
        nFrames = length(ts(1):te(1));
        % make figure
        show_inj_fig(D, ts, te, nFrames, sel, opt);
        % save figures
        figure(1); pause(0.5);
        printPDF( sprintf('%s/%sInjectionFigure', figDir, D.(fn{1}).pig) );
        figure(2);
        printPDF( sprintf('%s/%sInjectionFigurePlots', figDir, D.(fn{1}).pig) );
        pause(1);
    end
% -------------------------------------------------------------------------
%     % 2. Total change over cardiac cycle
    if do == 2
        clf();
        hold on;
        % otptions
        name = horzcat(char(pig), ' ',ref,' reference. Total change over average cardiac cycle');
        title( name );
        % make figure
        delta_heatmap(D, sel, opt);
        colorbar();
        % title(horzcat(num2str(sel), ' - ',D.(fn{sel}).name));
        colormap jet;
        fig = gcf;
        fig.Colormap(1,:) = [1 1 1] * 0.75;
        printPDF( sprintf('%s/tcoc/%sTCOCC', figDir, pig) );
        pause(1);
    end
    
% -------------------------------------------------------------------------

    % 2.5. Compare average CC for all sequences during saline injection
    % Use full period between saline injection and resuming ventilation.
    % This is the timeframe where there are the most cardiac cycles over
    % which to average!
    if do == 2.5
        switch pig
            case '8-2';     ts = [1,1,1,1];       te = [20,20,20,20];
            case '9-2';     ts = [1,1,1,1];       te = [20,20,20,20];
            case '10-2';    ts = [1,1,1,1,1,1];   te = [20,20,20,20,20,20];
            case '11-2';    ts = [1,1,1,1,1,1];   te = [20,20,20,20,20,20];
            case '12-2';    ts = [1,1,1,1,1,1];   te = [20,20,20,20,20,20];
        end % end switch
        
        cd('E:\University\Masters\HamburgBrain\Figures\saline_ensemble');
        % options
        opt = struct;
        opt.period = 4;
%         opt.sidelen = inf;
        opt.sidelen = 15;
        opt.align = 2;
        opt.normalize = true;
%         opt.normalize = false;
        
        % make figures
        cc_ensemble_during_saline(D, sel, opt);
        
        % save figures
        figure(1); 
%         title(TITLE); 
        colorbar(); pause(0.5);
        ttl = sprintf('%s/saline_ensemble/%s_injectionEnsemble_%sRef_normColour=%s', figDir, pig, ref, num2str(opt.normalize));
        printPDF(ttl);
        
%         figure(2);
%         title(TITLE);
%         printPDF(horzcat(char(pig), 'injection_ensemblePlots_baselineRef_normColour'));
    end
    
% -------------------------------------------------------------------------

    % 3. Compare average CC for all sequences with adjusted clim
    if do == 3
        figure(1); clf();
        figure(2); clf();
        % options
        opt.sidelen = inf;
        % opt.sidelen = 10;
        opt.align = 2;

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
        opt.normalize = true;
        compare_pre_inj_imgs(D, sel, opt);
        % save figures
        figure(1); 
%         title(TITLE); 
        colorbar(); pause(0.5);
        printPDF(sprintf('%s/ensemble/%s_ensemble', figDir, pig));
        figure(2);
%         title(TITLE);
        printPDF(sprintf('%s/ensemble/%s_ensemblePlots', figDir, pig));
    end
% -------------------------------------------------------------------------
% 4. Systole Diastole difference image
    if do == 4
        img1 = mk_diast_syts_dif_fig(D.seq1);
        img2 = mk_diast_syts_dif_fig(D.seq2);
        img3 = mk_diast_syts_dif_fig(D.seq3);
        imagesc([img1; img2; img3]); axis equal;
    end
% -------------------------------------------------------------------------
% 5. Assess localization accuracy
    if do == 5
        NORMALIZE = true;
        outImg = compare_pre_inj_imgs(D, sel, opt);
        nFrames = size(outImg.elem_data, 2) / length(sel);
        if NORMALIZE
            tempImg = outImg;
            slc     = [];
            start   = 1;
            stop    = nFrames;
            for i = 1: length(sel)
                tempImg.elem_data = outImg.elem_data(:, start:stop);
                if i == 1
                    slc = calc_colours(calc_slices(tempImg));
                else
                    slc(:,:,start:stop) = calc_colours(calc_slices(tempImg));
                end % end if
                start   = start + nFrames;
                stop    = stop + nFrames;
            end % end for
        else
            slc = calc_slices(outImg);
        end % end if
        segX = sum(seg, 1);
        segY = sum(seg, 2);
        
        brainTop    = find(segY ~= 0, 1, 'first');
        brainBot    = find(segY ~= 0, 1, 'last');
        brainLeft   = find(segX ~= 0, 1, 'first');
        brainRight  = find(segX ~= 0, 1, 'last');
        
        for r = 1:size(slc, 3)
            img1 = slc(:,:,r);
            Dxy = localization_err( img1(brainTop:brainBot,brainLeft:brainRight), cbv_img(pig) );
            disp(Dxy);
        end
    end % end if do 4

% -------------------------------------------------------------------------
%     % revert landmarks
%     switch pig
%         case '8-2';     D.seq.eit.apn = temp; D.seq2.eit.apn = temp2;
%         case '11-2';    D.seq3.eit.apn = temp;
%         case '12-2';    D.seq3.eit.apn = temp;
%     end % end switch
end % end for

% ======================================================================= %

function perfImgs = mk_diast_syts_dif_fig(seq)
    global seg;
    binSeg = (seg > 0);
    apn = seq.eit.apn;
    inj = seq.eit.inj;
    
    thisImgr = seq.imgr;
    
    systole = seq.eit.peaks';
    diastole = seq.eit.vals';
    triad = systole;
    dIdx = sum( ((triad >= apn) + (triad <= inj)), 2) == 2;
    triad = triad(dIdx, :);
    
    perfImgElems = zeros(size(thisImgr.elem_data, 1), length(triad) - 1);
    for i=1:length(triad) - 1
        temp = thisImgr;
        temp.elem_data = thisImgr.elem_data(:, triad(i) : triad(i+1));
        elemRange = range(temp.elem_data);
        diast = find(elemRange == min(elemRange));
        syst = find(elemRange == max(elemRange));
        perfImgElems(:, i) = temp.elem_data(:, syst) - temp.elem_data(:, diast);
    end
    thisImgr.elem_data = mean(perfImgElems, 2);
    thisImgr.calc_colours.clim = max(thisImgr.elem_data(:));
    perfImgs = show_slices(thisImgr);
    
    
    
%     systole = seq.eit.peaks;
%     sIdx = ((systole >= apn) + (systole <= inj)) == 2;
%     systole = systole(sIdx);
%     systImgr = seq.imgr;
%     systImgr.elem_data = seq.imgr.elem_data(:, systole);
%     systImgs = calc_slices(systImgr);
%     
%     diastole = seq.eit.vals(1,:);
%     dIdx = ((diastole >= apn) + (diastole <= inj)) == 2;
%     diastole = diastole(dIdx);
%     diastImgr = seq.imgr;
%     diastImgr.elem_data = seq.imgr.elem_data(:, diastole);
%     diastImgs = calc_slices(diastImgr);
%     
%     tmp = diastImgr;
%     tmp.elem_data = systImgr.elem_data - diastImgr.elem_data;
%     perfImgs = calc_slices(tmp);
%     perfImgs = systImgs - diastImgs;

%     outImgs = zeros(32, 64, size(perfImgs, 3));
%     mask = find(binSeg);
%     for i = 1:size(perfImgs, 3)
% %         tmp = perfImgs(:,:, i);
%         tmp = trim_img(perfImgs(:,:, i));
%         mask = find(~isnan(tmp));
%         tmp(mask) = normalize(tmp(mask), 'range');
%         outImgs(:,:, i) = tmp;
%     end
%     img = squeeze(mean(perfImgs, 3));
%     img = squeeze(mean(outImgs, 3));
%     img1 = max(abs(outImgs), [], 3);
%     MM = max(outImgs, [], 3);
%     mm = min(outImgs, [], 3);
%     
%     tmp = abs(MM) > abs(mm);
%     img2 = squeeze(std(outImgs, [], 3));
    
%     img2 = squeeze(mean(outImgs, 3));
%     img = [img1, img2];
%     img = squeeze(std(outImgs, [], 3));
%     img = trim_img(img .* binSeg);
%     img = trim_img(img);
end

% ======================================================================= %

function newImg = trim_img(img)
    global seg;
    NROWSINIMG = 32;
    binSeg = (seg > 0);
    segBoundaries = sum(binSeg, 2) == 0;
    topRowSeg = find(segBoundaries == 0, 1, 'first');
    botRowSeg = find(segBoundaries == 0, 1, 'last');
    toDivide = NROWSINIMG - (botRowSeg - topRowSeg + 1);
    topIdx = topRowSeg - floor(toDivide / 2); % if toDivide is odd, will give more rows to bottom (as in MRI images)
    botIdx = botRowSeg + ceil(toDivide / 2);
    nRows = size(img, 1);
    if topIdx < 1
        botIdx = botIdx + (1 - topIdx); % give these rows to bottom if too high on image
        topIdx = 1;
    elseif botIdx > nRows
        topIdx = topIdx + (nRows - botIdx);
        botIdx = nRows;
    end % end if
    newImg = img(topIdx: botIdx, :, :);
end

% ======================================================================= %

function temp = adjust_for_poi(seq, start, stop)
    % Create a copy of seq with adjusted cardiac cycle landmarks for selected start and stop data
    % analysis window. The elem_data within temp will also be trimmed to
    % encompass only the region defined by start and stop
    shift = start - 1;
    useLndmrk = ( (seq.eit.vals(1, :) >= start) + (seq.eit.vals(2, :) <= stop) ) == 2;
    temp = seq;
    temp.imgr.elem_data = seq.imgr.elem_data(:, start: stop);
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

function cc_ensemble_during_saline(D, sel, opt)
    IMGCOLS = 10;
    NUMCYCLES = 25;
    fn = fieldnames(D);
    period = opt.period;
    if nargin == 1
        sel = 1: length(fn);
    elseif isempty(sel)
        sel = 1: length(fn);
    end % end if    
    MC = struct;
    SC = struct;
    maxCcLen = 0;
    CcLens  = zeros(1, length(sel));
    for j = 1:length(sel)
        i = sel(j);
        seq = D.(fn{i});
        
        % define injection event and period of clean data after
        [start, stop] = define_period(seq, period);
%         start = start + round(seq.eit.fs * 10);
        
        % get sequence object with elem_data and CC landmarks within
        % injection period
        temp = adjust_for_poi(seq, start, stop);
        
        % Choose at most last 45 peaks and valleys
        peakCount = length(temp.eit.peaks);
        if peakCount > NUMCYCLES
%             startIdx = peakCount - NUMCYCLES + 1;
            startIdx = peakCount - NUMCYCLES + 1;
            endIdx = startIdx + NUMCYCLES - 1;
            temp.eit.peaks = temp.eit.peaks(startIdx : endIdx);
            temp.eit.vals = temp.eit.vals(:, startIdx : endIdx);
            
            temp.perf.peaks = temp.perf.peaks(startIdx : endIdx);
            temp.perf.vals = temp.perf.vals(:, startIdx : endIdx);
        end
        
        % Get imgr object whose elem_data is the average of CCs within
        % injection frame. CLIM is adjusted to max of elem_data.
        [MC.(fn{i}), SC.(fn{i})] = get_mean_cc( temp, opt );
        CcLens(j) = size(MC.(fn{i}), 2);
        maxCcLen = max( maxCcLen, CcLens(j) );
    end % end for
    
    fn = fieldnames(MC);
    nElems = size(MC.(fn{j}), 1);
    addFrames = abs(CcLens - maxCcLen);
    outImg = D.(fn{j}).imgr;
    outImg.elem_data = [];
    outImg.show_slices.img_cols = IMGCOLS;
    
    for i = 1:length(sel)
        if addFrames(i) > 0
            paddedAvCc = [ zeros(nElems, addFrames(i)), MC.(fn{i}) ];
        else
            paddedAvCc = MC.(fn{i});
        end % end if
        outImg.elem_data = [outImg.elem_data, paddedAvCc];
        
        temp = D.(fn{i}).imgr;
        temp.elem_data = SC.(fn{i});
        temp.calc_colours.clim = max(temp.elem_data, [], 'all');
        figure(2 + i);
        show_slices(temp);
        title(sprintf('Standard deviation image: %s', fn{i}));
        printPDF(sprintf('%s %s ensemble SD', D.(fn{i}).pig, fn{i}));
    end % end for
    
    outImg.calc_colours.clim = max(outImg.elem_data, [], 'all');
    
    show_slice_and_brain_seg_z(outImg, fn, [], opt);
        
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

function [meanCycle, stdEnsemble] = get_mean_cc(seq, opt)
    
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
        systFrm = lpad(1) + ls(1); % location of systolic frame
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
    meanCycle = z_frames(:, 1: longestCC);
    nContributed = zeros(1, longestCC);
    % create full image set for each CC that is the difference between each
    % frame and the mean diastole for that CC.
    for i = 1:nCycles
        startIdx = lpad(i) + 1;
        endIdx = startIdx + ccLen(i) - 1;
        imgs = seq.imgr.elem_data( :, use_diast(1, i): use_diast(2, i) ) - z_diast(:, i); % elem_data has mean of flanking diast subtracted from cc sequence
        z_frames(:, startIdx: endIdx , i) = imgs;
        meanCycle(:, startIdx: endIdx) = meanCycle(:, startIdx: endIdx ) + imgs;
        nContributed( startIdx: endIdx) = nContributed( startIdx: endIdx) + 1;
    end % end for
    
    longestCC = length(showRng);
    stdEnsemble = std(meanCycle(:, showRng), [], 2);
    meanCycle = meanCycle(:, showRng) ./ nContributed(showRng);
    z_frames = z_frames(:, showRng, :);
    z_frames = reshape( z_frames, size(z_frames, 1), longestCC * nCycles );
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
    fprintf('Number of cycles %s: %.0f\n',seq.name, nCycles)
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
    if length(period) > 1
        injOpt = period(2: end);
        period = period(1);
    else
        injOpt = [10];
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
        [start, stop] = define_period(seq, opt.period(i));
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

function outImg = compare_pre_inj_imgs(D, sel, opt)
    IMGCOLS = 10;
    fn      = fieldnames(D);
    period  = opt.period;
    if nargin == 1
        sel = 1:length(fn);
    elseif isempty(sel)
        sel = 1:length(fn);
    end % end if    
    MC      = struct;
    maxCcLen = 0;
    CcLens  = zeros(1, length(sel));
    for j = 1:length(sel)
        i = sel(j);
        seq = D.(fn{i});
        [start, stop] = define_period(seq, period(i));        
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
    
    if nargout == 0
        show_slice_and_brain_seg_z(outImg, fn, [], opt);
    end % end if
        
end % end function

% ======================================================================= %

function cmap = confg_cmap()

    calc_colours('cmap_type', 'blue_black_red');
    colormap(calc_colours('colourmap'));
    cmap = colormap * 2;
    cmap(cmap > 1) = 1;
    black = find(sum(cmap == [0,0,0], 2) == 3);

    white_grad = linspace(0,1,black-1)';
    cmap(black:end,1) = white_grad;
    cmap(black:end,2) = white_grad;

    white_grad = linspace(0, 1, black)';
    white_grad = flipud(white_grad);
    cmap(1:black, 2) = white_grad;
    cmap(1:black, 3) = white_grad;

end % end funcion

% ======================================================================= %

function bigFig()
    figure('units','normalized','outerposition',[0 0 1 1]);
end % end function

% ======================================================================= %

function show_inj(seq, ts, te)
    [start, stop] = define_period(seq, [3, ts, te]);
    temp = seq.imgr;
    temp.elem_data = seq.imgr.elem_data(:, start:stop);
    temp.calc_colours.clim = max(temp.elem_data, [] ,'all');
    temp.show_slices.img_cols = round(seq.eit.fs);
    show_slices(temp);
end % end function

% ======================================================================= %

function show_inj_fig(D, ts, te, nFrames, sel, opt)    
    IMGCOLS = 10;
    fn = fieldnames(D);
    fn2 = cell(length(sel),1);
    ts = ts(sel);
    te = te(sel);    
    outImg = D.seq1.imgr;
    injFrames = zeros( size(outImg.elem_data, 1), nFrames, length(sel) );

    for i = 1: length(sel)
        seq = D.(fn{sel(i)});
        [start, stop] = define_period(seq, [3, ts(i), te(i)]);
        
        injFrameIdx = round( linspace(start, stop, nFrames + 1) );        
        injFrameIdx = injFrameIdx(1: nFrames); % each frame is the start of a second this way.
        fn2{i} = fn{sel(i)};
        % want to increase time resolution of line plots for brain
        % conductivity
        injFrameIdx2 = start: stop; % data to use for plotting
        if i == 1
            injFrames2 = zeros( size(outImg.elem_data, 1), length(injFrameIdx2), length(sel) );
        end % end if
        tmp = seq.imgr.elem_data(:, injFrameIdx);
        tmp = tmp - tmp(:, 1);
        injFrames(:,:,i) = tmp;
        
        tmp = seq.imgr.elem_data(:, injFrameIdx2);
        tmp = tmp - tmp(:, 1);
        injFrames2(:,:,i) = tmp;
    end % end for
    
    outImg.elem_data = reshape(injFrames, size(outImg.elem_data, 1), nFrames * i);
    outImg.show_slices.img_cols = IMGCOLS;
    
    % new stuff with higher temporal res for slices in inj
    plotimg = D.seq1.imgr;
    plotimg.elem_data = reshape(injFrames2, size(plotimg.elem_data, 1), length(injFrameIdx2) * i); % num elems x (n seq x length each series)
    plotimg.show_slices.img_cols = length(injFrameIdx2);  
    
    show_slice_and_brain_seg_z(outImg, fn2, plotimg, opt);
    figure(1);
%     title( sprintf('Subject %s - Reconstructed images from %i - %i seconds after saline bolus injection', char(D.(fn{i}).pig), ts(i), te(i)) );
    figure(2);    
%     title( sprintf('Subject %s - Reconstructed images from %i - %i seconds after saline bolus injection', char(D.(fn{i}).pig), ts(i), te(i)) );
end % end function

% ======================================================================= %

function show_slice_and_brain_seg_z(img, fn, plotimg, opt)
    global seg;
    if isfield(opt, 'normalize')
        NORMALIZE = opt.normalize;
    else
        NORMALIZE = true;
    end % end if
    
    NROWSINIMG = 32;
    % show slices
    figure(1);
    binSeg = (seg > 0) * 5; % binary segmentation (all) 5x multiplier for brain in show slices.
    binSeg(binSeg == 0) = 1; % multiply all other image elements by 1.
    
    % calculate number of rows in image for each sequence. Trim end of
    % sequences if needed so that they fit in an integer number of rows.
    framePerImg = size(img.elem_data, 2) / length(fn);
    rowPerSeq = framePerImg / img.show_slices.img_cols;
    if ~isinteger(rowPerSeq)
        rowPerSeq = floor(rowPerSeq);
        newFrameCount = rowPerSeq * img.show_slices.img_cols;
        newElemData = [];
        start = 1;
        for i = 1: length(fn)
            stop = start + newFrameCount - 1;
            seq_data = img.elem_data(:, start: stop);
            newElemData = [newElemData, seq_data];
            start = start + framePerImg;
        end % end for
        img.elem_data = newElemData;
    end % end if
    slices = calc_slices(img) .* binSeg;    
    
    % calculate normalized colours across all sequences if true.
    if NORMALIZE
        clr_slices = calc_colours(slices);
    else
        clr_slices = [];
        start = 1;
        stop = rowPerSeq * img.show_slices.img_cols;
        for i = 1: length(fn)    
            thisSlice = slices(:, :, start: stop);
            clr_slices(:, :, start: stop) = calc_colours(thisSlice);
            start = start + rowPerSeq * img.show_slices.img_cols;
            stop = stop + rowPerSeq * img.show_slices.img_cols;
        end
    end
    
    % trim EIT images so they are roughly the same height as the perfusion
    % MR images.
    fivexSlc = trim_img(clr_slices);
    
    % check to see if a full row of images is blank and remove if so.
    rowPerSeq2 = repmat(rowPerSeq, length(fn), 1); 
    elimRows = squeeze( sum(fivexSlc, [1,2]) / ( size(fivexSlc,1) * size(fivexSlc,2) * 128) )';
    elimRows = reshape(elimRows, img.show_slices.img_cols, rowPerSeq*length(fn));
    elimRows = sum(elimRows, 1) == img.show_slices.img_cols;
    if sum(elimRows) >= 1
        elimRowsIdxS = elimRows*1 .* (1: 10: 10 * rowPerSeq * length(fn));
        elimRowsIdxS = elimRowsIdxS(elimRowsIdxS >0);
        elimRowsIdxE = elimRowsIdxS + img.show_slices.img_cols - 1;
        for i = 1:length(elimRowsIdxS)
            fivexSlc(:, :, elimRowsIdxS(i): elimRowsIdxE(i)) = nan;
        end % end for i
        fivexSlc = fivexSlc(~isnan(fivexSlc)); % remove those rows
        fivexSlc = reshape(fivexSlc, NROWSINIMG, nRows, []); % reshape for show
        rng = [1,rowPerSeq];
        idx = find(elimRows); % get row numbers of all blank rows
        for i = 1:length(fn)
            rowPerSeq2(i) = rowPerSeq2(i) - sum( (idx >= rng(1)) + (idx <= rng(2)) == 2 ); % eliminate rows from rowsPerSeg to allow proper axis labelling.
            rng = rng + rowPerSeq;
        end % end for i
    end % end if
    
    imagesc(mk_mosaic(fivexSlc , 0, [], img.show_slices.img_cols)); 
    axis equal; axis image;
    
    % set ticklabels to denote start of each sequence/recording
    ax = gca;
    ax.YTick = NROWSINIMG / 2;
    for i = 2: length(fn)
        ax.YTick = [ax.YTick, ax.YTick(i - 1) + rowPerSeq2(i-1) * NROWSINIMG];
    end % end for
    ax.YTickLabel = cell(length(fn), 1);
    for i = 1:length(fn)
%         ax.YTickLabel{i} = sprintf('Sequence %i', i);
        ax.YTickLabel{i} = '';
    end % end for
    ax.XTickLabel = {};
    
    % show line plots of brain conductivity
    if isstruct(plotimg) % if we want the higher temporal res plots
        img = plotimg;
        nFrames = img.show_slices.img_cols;
    else
        nFrames = img.show_slices.img_cols * rowPerSeq;
    end
    
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
        plot_brain_seg(clr_slices(:, :, start: stop));
        fg = gcf;
        fg.Children(3).XGrid = 'on';
        if fg.Children(3).YLim(1) < ymin
            ymin = fg.Children(3).YLim(1);
        end % end if
        if fg.Children(3).YLim(2) > ymax
            ymax = fg.Children(3).YLim(2);
        end % end if
        
    end % end for

    for i= 3: 3: length(fn) * 3
        fg.Children(i).YLim = [ymin * .95 ymax * 1.05];
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

% ======================================================================= %

function plot_brain_seg(slc)
    global seg;
    binSeg = (seg > 0)*1; % binary segmentation (all)
    slc(isnan(slc)) = 0;
    nPix = sum(binSeg,'all');
    
    for i = 1:size(slc, 3)
        tmp = slc(:, :, i);
        tmp(binSeg == 0) = 0;
%         tmp(binSeg == 1) = normalize(tmp(binSeg == 1), 'range');
        slc(:, :, i) = tmp;
    end
    
    sumBrain = squeeze( sum( slc .* binSeg, [1,2] ) );
    sumBrain(sumBrain == nPix) = nan;
    % plot global signal first, dark blue colour.
    tmp = sumBrain ./ nPix;
    plot(tmp, 'linewidth', 2);
    legend();
    xlim([0.5, size(slc,3)+0.5]);
    hold on;
    % plot each of the 6 regions next.
    for j = 1:max(seg, [], 'all')
        mask = (seg == j) * 1;
        nPix = sum(mask, 'all');
        dd = squeeze( sum( slc .* mask, [1,2] ) );
        dd(dd == nPix)  = nan;
        plot( dd./nPix, 'linewidth', 2);  
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

% ======================================================================= %

function printPDF(filename)
    h = gcf;
    set(h, 'PaperUnits', 'centimeters');
    set(h, 'Units', 'centimeters');
    pos = get(h, 'Position');
    set(h, 'PaperSize', [pos(3) pos(4)]);
    set(h, 'PaperPositionMode', 'manual');
    set(h, 'PaperPosition', [0 0 pos(3) pos(4)]);
    print('-dpdf', filename);
end % end function

% ======================================================================= %
