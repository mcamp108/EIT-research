% -------------------------------------------------------------------------
% Description: This script is used for analyzing and visualizing swine
% stroke data.
% -------------------------------------------------------------------------
% Author:
%   Mark Campbell
%   Carleton University
%   markacampbell@cmail.carleton.ca
% -------------------------------------------------------------------------

run 'myStartup.m';
% bigFig();
maxsz= 0.2; maxh= 2; imgsize= [64 64]; 
% pigs= ["8.2","9.2","10.2","11.2","12.2"];
pig= "12.2";
[fmdl, imdl]= mk_pighead_fmdl(maxsz, maxh, imgsize, pig);
% Load data
ref = 'self';
% ref = 'baseline';
D= load_HamburgBrain_data(pig, ref);
fn= fieldnames(D);
cd(horzcat('C:\Users\Mark\Documents\GraduateStudies\LAB\HamburgBrain\Figures\', char(pig)));
suffix= date;

% configure colormap
% cmap = confg_cmap();

% 1. SHOW EIT AND PERFUSION DATA WITH PERFUSION ANNOTATIONS
cd(horzcat('C:\Users\Mark\Documents\GraduateStudies\LAB\HamburgBrain\Figures\', char(pig)));
%%
cd seqData;
for i= 1:numel(fn)
    plot_seq_data(D.(fn{i}));
    saveas( gcf, horzcat(D.(fn{i}).name, suffix, '.svg') );
end % end for
close all
cd ../
%% 2. Injection Images Figure
cd('C:\Users\Mark\Documents\GraduateStudies\LAB\HamburgBrain\Figures\paper');
nFrames = 15;
switch pig
    case "8.2";     ts = [0,0,0,0];             te = [10,10,10,10]; % use seq2 as ref
    case "9.2";     ts = [10,10,10,10];         te = [25,25,25,25];
    case "10.2";    ts = [10,10,10,10,10,10];   te = [25,25,25,25,25,25];
    case "11.2";    ts = [10,10,10,10,10,10];   te = [25,25,25,25,25,25];
%     case "12.2";    ts = [10,10,10,10,10,10];   te = [25,25,25,25,25,25];
    case "12.2";    ts = [1,1,1,1,1,1];   te = [25,25,25,25,25,25];
end % end switch

bigFig();
show_inj_fig(D, ts, te, nFrames);
saveas( gcf, sprintf('z %s injection figure simple.svg', char(D.(fn{1}).pig)) );

%% 3. Total change over cardiac cycle
cd('C:\Users\Mark\Documents\GraduateStudies\LAB\HamburgBrain\Figures\paper');
switch pig
    case "8.2";     opt.period = [2,2,1];        
        temp = D.seq1.eit.apn; D.seq1.eit.apn = 1868; 
        temp2 = D.seq2.eit.apn; D.seq2.eit.apn = 1068;
        sel = [2,3,4];
    case "9.2";     opt.period = [0,1,1,1];
        sel = [1,2,3];
    case "10.2";    opt.period = [1,1,1,1,1,1];
        sel = [1,3,5];
    case "11.2";    opt.period = [1,2,0,2,1,1];
        sel = [1,3,6];
        temp = D.seq3.eit.apn; D.seq3.eit.apn = 843;
    case "12.2";    opt.period = [2,0,0,0,2,0];
        temp = D.seq3.eit.apn; D.seq3.eit.apn = 939; % 1223; 
        sel = [2,3,5];
end % end switch

bigFig();
hold on;
name = horzcat(char(pig), ' ',ref,' reference. Total change over average cardiac cycle');
sgtitle( name );
% sel = 1:length(fn);
delta_heatmap(D, sel, opt);
colorbar();
% title(horzcat(num2str(sel), ' - ',D.(fn{sel}).name));
colormap jet;
fig = gcf;
fig.Colormap(1,:) = [1 1 1] * 0.9;
saveas( gcf, horzcat('z ', char(pig),' TCOCC period ',num2str(opt.period),' sequences ',num2str(sel), ' simple.svg') );
switch pig
    case "8.2";     D.seq.eit.apn = temp; D.seq2.eit.apn = temp2;
    case "11.2";    D.seq3.eit.apn = temp;
    case "12.2";    D.seq3.eit.apn = temp;
end % end switch



%% 4. Compare average CC for all sequences with adjusted clim
cd('C:\Users\Mark\Documents\GraduateStudies\LAB\HamburgBrain\Figures\paper');
opt.sidelen = 5;
switch pig
    case "8.2";     opt.period = [2,2,2,1];
    case "9.2";     opt.period = [0,1,1,1];
    case "10.2";    opt.period = [1,1,1,1,1,1];
    case "11.2";    opt.period = [1,2,0,2,1,1];
    case "12.2";    opt.period = [2,0,0,0,2,0];
end % end switch
switch pig
    case "8.2"
        temp = D.seq1.eit.apn; D.seq1.eit.apn = 1868; 
        temp2 = D.seq2.eit.apn; D.seq2.eit.apn = 1068;
        sel = [2,3,4];
        TITLE = '8.2 Ensemble Average of Pulsatile Signal. baseline (top), after embolism (middle), 30 minutes after embolism (bottom)';
    case "9.2"
        sel = [1,2,3];
        TITLE = '9.2 Ensemble Average of Pulsatile Signal. baseline (top), baseline 2 (middle), after embolism (bottom)';
    case "10.2"
        opt.select = 1;
        sel = [1,5];
%         TITLE = '10.2 Schleuse Ensemble Average of Pulsatile Signal. baseline (top), after stroke induction (middle), 6 hours after stroke induction (bottom)';
        TITLE = '10.2 Schleuse Ensemble Average of Pulsatile Signal. baseline (top), 6 hours after stroke induction (bottom)';
    case "11.2"
        temp = D.seq3.eit.apn; D.seq3.eit.apn = 843;
        sel = [1,3,6];
        TITLE = '11.2 Schleuse Ensemble Average of Pulsatile Signal. baseline (top), after stroke induction (middle), 4 hours after stroke induction (bottom)';
    case "12.2"
        temp = D.seq3.eit.apn; D.seq3.eit.apn = 939; % 1223; 
        sel = [2,3,5];
        TITLE = '12.2 Schleuse Ensemble Average of Pulsatile Signal. baseline (top), after stroke induction (middle), 3.5 hours after stroke induction (bottom)';
end % end switch

saveName = horzcat('z simple ', char(pig), ' ensemble. seq', num2str(sel), ' ', ref, '.svg');
opt.align = 2;
bigFig();
compare_pre_inj_imgs(D, sel, opt);

title(TITLE); colorbar();
saveas( gcf, saveName );

switch pig
    case "8.2";     D.seq.eit.apn = temp; D.seq2.eit.apn = temp2;
    case "11.2";    D.seq3.eit.apn = temp;
    case "12.2";    D.seq3.eit.apn = temp;
end % end switch





%% 2. Ensemble average perfusion images
cd(horzcat('C:\Users\Mark\Documents\GraduateStudies\LAB\HamburgBrain\Figures\', char(pig)));
opt.av = 'mean';
opt.align = 2;
% opt.sidelen
for i = 1:numel(fn)
    bigFig();
    show_pre_inj_img(D.(fn{i}), opt);
end % end for


%% 3.5
bigFig();
compare_post_vnt_imgs(D, [5], opt);
%% 4. Compare images after bolus injection
cd(horzcat('C:\Users\Mark\Documents\GraduateStudies\LAB\HamburgBrain\Figures\', char(pig)));
show_inj_imgs(D.seq1, opt);
%% 5. Show mean cycle from pre and post injection time periods
cd(horzcat('C:\Users\Mark\Documents\GraduateStudies\LAB\HamburgBrain\Figures\', char(pig)));
opt.av = 'mean';
for i = 1:numel(fn)
    bigFig();
    compare_pre_post_inj(D.(fn{i}), opt);
    saveas( gcf, horzcat(D.(fn{i}).name, '_PrePostInj_', suffix, '.svg') );
end % end for
%% 6. Visualize different timeframes after injection
cd(horzcat('C:\Users\Mark\Documents\GraduateStudies\LAB\HamburgBrain\Figures\', char(pig)));
ts = 1;
te = 10;
bigFig();
for i = 1:numel(fn)
    show_inj(D.(fn{i}),ts,te);
    title( horzcat(num2str(i), '_',D.(fn{i}).name, '1-10 seconds Post-Injection. 1 row = 1 second') );
%     saveas( gcf, horzcat(num2str(i), '_',D.(fn{i}).name,'1-10secondsPostInj.svg') );
end % end for
%% 8. Injection Images Exploration
switch pig
    case "8.2";     ts = [0,0,0,0];         te = [10,10,10,10];
    case "9.2";     ts = [10,10,10,10];         te = [25,25,25,25];
    case "10.2";    ts = [10,10,10,10,10,10];   te = [25,25,25,25,25,25];
    case "11.2";    ts = [10,10,10,10,10,10];   te = [25,25,25,25,25,25];
    case "12.2";    ts = [10,10,10,10,10,10];   te = [25,25,25,25,25,25];
end % end switch
for i=1:length(fn)
    bigFig();
    show_inj(D.(fn{i}), ts(i), te(i));
    colorbar();
    title( sprintf('%i - %s %i - %i seconds Post-Injection. 1 row = 1 second', i, D.(fn{i}).name, ts(i), te(i)) );
end % end for
%% VIDEO OF RECONSTRUCTED IMAGE AND BRAIN SEGMENTATION

cd(horzcat('C:\Users\Mark\Documents\GraduateStudies\LAB\HamburgBrain\Figures\', char(pig)));
for i= [1,3,5]
    start= D.(fn{i}).eit.inj;
    stop= D.(fn{i}).eit.inj+ 1000;
    mk_vid(D.(fn{i}), start, stop, suffix);
end % end for

%% EIT, PREFUSION, AND RECONSTRUCTED IMAGE ENSEMBLES IN ONE FIGURE
cd(horzcat('C:\Users\Mark\Documents\GraduateStudies\LAB\HamburgBrain\Figures\', char(pig)));
opt.pv= 3;
opt.usefData= 1;
opt.plotLM= 2;
opt.ensemble= 'each';
cd 'ensemble/rm CI greater than 400/each';
for i= 1:numel(fn)
    for j= 0:4
        opt.section= j;
        show_ensemble(D.(fn{i}), opt);
        figureTitle= horzcat(num2str(j), ' ', D.(fn{i}).name, ' - ensemble', suffix);
        print_convert(horzcat(figureTitle, '.png'));
    end % end for
    close all
end % end for
cd ../
%% TIME TO PEAK
opt.sel= [0 0 0 1];
opt.plotLM= 1;
ttp1= [];
for i= [1, 2, 3, 4]
%     opt.start= 50;
%     opt.stop= D.(fn{i}).eit.inj;
    opt.start= D.(fn{i}).eit.inj;
    opt.stop= D.(fn{i}).eit.inj+ 1000;
    ttp1= [ttp1, get_ttp(D.(fn{i}), 2, opt)];
%     title("time to peak for " + D.(fn{i}).name);
end % end for
figure; imagesc(ttp1); axis 'equal'; colorbar;

%% SHOW MEAN FRAME
cd meanFrame
for i= 1:numel(fn)
    meanFrame = mean(calc_slices(D.(fn{i}).imgr), 3);
    figure('units','normalized','outerposition',[0 0 1 1]);
    imagesc(meanFrame);
    title(horzcat(D.(fn{i}).name, suffix));
    colorbar;
    axis equal
    
    [imageDenoised, betheE]= mrfDeNoiseV3(meanFrame, 3);
    figure('units','normalized','outerposition',[0 0 1 1]);
    imagesc(imageDenoised);
    
%     print_convert(char(D.(fn{i}).name+ suffix+ ".png"));
end % end for
close all
cd ../


%% Plot selected pixels of brain segmentation over time
acvOpt.pixels= (1:100);
acvOpt.plotLM= 1;
result= plot_segmentations(D.seq1, acvOpt);

mtx= zeros(1, size(result, 1));
for r= 1:size(result, 1)
    row= result(r, :);
    mtx(r)= find(row== max(row));
end % end for
imagesc(mtx); % loks kinda cool

%%
opt.start= D.seq2.eit.inj;
opt.stop= D.seq2.eit.inj+ 600;
ttp1= get_ttp(D.seq2, 1, opt);

opt.start= D.seq6.eit.inj;
opt.stop= D.seq6.eit.inj+ 600;
ttp2= get_ttp(D.seq6, 1, opt);
imagesc([ttp1;ttp2])

figure; imagesc(ttp1- ttp2);
imageOrig= ttp1- ttp2;
out_img= mrfDeNoiseV2(ttp1, 5);
figure; image(out_img);

for i= 1: numel(fn)
    ensemble= get_ensembles(D.(fn{i}), 2, opt);
    figure('units','normalized','outerposition',[0 0 1 1]);
    sgtitle(D.(fn{i}).name);
    
    subplot(2, 2, [1,3]);
    show_slices(ensemble.imgr_ensemble);
    title("reconstructions averaged over cardiac cycles");
    
    subplot(2, 2, 2);
    xax= (1: size(ensemble.eit_ensemble, 2))/ D.(fn{i}).eit.fs;
    plot(xax, mean(ensemble.eit_ensemble, 1));
    title("Average EIT signal over all cardiac cycles");
    xlabel("Time (s)");
    
    subplot(2, 2, 4);
    xax= (1: size(ensemble.perf_ensemble, 2))/ D.(fn{i}).perf.tickrate;
    plot(xax, mean(ensemble.perf_ensemble, 1));
    title("Average perfusion signal over all cardiac cycles");
    xlabel("Time (s)");
    print_convert(char(D.(fn{i}).name+ " ensemble simplified model.png"));
end % end for
cd ../

%%


function simulate_ischemic_zone(pig)

img = mk_image(fmdl, 0.41); % Background conductivity is scalp
img.elem_data([fmdl.mat_idx{1}]) = 0.41;    %   1: scalp        0.41
img.elem_data([fmdl.mat_idx{2}]) = 0.016;   %   2: skull        0.016
img.elem_data([fmdl.mat_idx{3}]) = 0.47;    %   3: grey matter  0.47
img.elem_data([fmdl.mat_idx{4}]) = 0.0001;  %   4: air          0.0001
img.fwd_solve.get_all_meas = 1;
img.fwd_model.stimulation = imdl.fwd_model.stimulation;
img.fwd_model.normalize_measurements = 0;
vh = fwd_solve(img);

figure();
img2 = img;
img2.elem_data([fmdl.mat_idx{1}]) = 0.41;    %   1: scalp        0.41
img2.elem_data([fmdl.mat_idx{2}]) = 0.016;   %   2: skull        0.016
img2.elem_data([fmdl.mat_idx{3}]) = 0.47;    %   3: grey matter  0.47
img2.elem_data([fmdl.mat_idx{4}]) = 0.1;  %   4: air          0.0001
vi = fwd_solve(img2);
imgr = inv_solve(imdl, vh, vi);
show_fem(imgr);

nodes = fmdl.elems(fmdl.mat_idx{1}, :);
coors = fmdl.nodes(nodes,:);
end % end function



function compare_pre_post_inj(seq, opt)

    ax1 = subplot(2,1,1);
        show_pre_inj_img(seq, opt);
    ax2 = subplot(2,1,2);
        show_inj_imgs(seq, opt);
    
    for i = [ax1, ax2]
        i.Position(1) = 0.05;
        i.Position(3) = 0.9;
    end % end for
    
    ax2.Position(2) = 0.35;
    
end % end function

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
    title(horzcat(seq.name, " Post-injection conductivity from individual and ensemble-averaged (bottom) Cardiac Cycles"));
    
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
    title(horzcat(seq.name, " Pre-injection conductivity from individual and ensemble-averaged (bottom) Cardiac Cycles"));

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
        case 0 % start to inj
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
        [start,stop] = define_period(seq, opt.period(j));
        temp = adjust_for_poi(seq, start, stop);
        meanCC = get_mean_cc( temp, opt );
        heatmapFrames(:,j) = max(meanCC, [], 2) - min(meanCC, [], 2);
    end % end for
    
    clim = max(heatmapFrames, [], 'all'); % adjust clim
    temp = D.(fn{1}).imgr;
    temp.calc_colours.clim = clim;
    temp.elem_data = heatmapFrames;
    temp.show_slices.img_cols = 3;
    img = show_slices(temp);
    assert(min(img(img~=1)) >= 128, "Min error");
    img = img - 127;
    img(img < 0) = 0;
    image(img*2);

end % end function

% ======================================================================= %

function compare_pre_inj_imgs(D, sel, opt)

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
        [start,stop] = define_period(seq, period(j));        
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

function show_inj_fig(D, ts, te, nFrames)
    
    fn = fieldnames(D);
    switch D.seq1.pig
        case "8.2"; ts = ts(1:3); te = te(1:3); fn = {'seq1','seq3','seq4'};
    end
    
    assert( length(fn) == length(ts), 'There must be one time slice pair per sequence!');
    assert( length(fn) == length(te), 'There must be one time slice pair per sequence!');
    temp = D.seq1.imgr;
    injFrames = zeros( size(temp.elem_data, 1), nFrames, length(fn) );

    for i=1:length(fn)
        seq = D.(fn{i});
        [start,stop] = define_period(seq, [3, ts(i), te(i)]);
        injFrameIdx = round( linspace(start, stop, nFrames+1) );
        injFrameIdx = injFrameIdx(1:nFrames); % each frame is the start of a second this way.
        injFrames(:,:,i) = seq.imgr.elem_data(:, injFrameIdx);
    end % end for

    temp.elem_data = reshape(injFrames, size(temp.elem_data, 1), nFrames * i);
    temp.show_slices.img_cols = nFrames;
    temp.calc_colours.clim = max(temp.elem_data, [] ,'all');
    show_slices(temp);

    colorbar();
    title( sprintf('Subject %s - Reconstructed images from %i - %i seconds after saline bolus injection', char(D.(fn{i}).pig), ts(i), te(i)-1) );
    ax = gca;
    ax.Visible = 'on';
    % Y labels
    seqYTicksRef = linspace(ax.YLim(1), ax.YLim(2), i*2+1);
    ax.YTick = seqYTicksRef(2:2:length(seqYTicksRef));
    for j = 1:length(ax.YTick)
        ax.YTickLabel{j} = sprintf('Sequence %i', j);
    end % end for
    % X labels
    seqXTicksRef = linspace(ax.XLim(1), ax.XLim(2), nFrames*2+1);
    ax.XTick = seqXTicksRef(2:2:length(seqXTicksRef));
    for j = 1:length(ax.XTick)
        ax.XTickLabel{j} = sprintf('%is', ts(1)+j-1);
    end % end for
    ax.TickLength = [0,0];
    ax.FontSize = 15;

end % end function

% ======================================================================= %




