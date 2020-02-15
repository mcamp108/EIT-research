% -------------------------------------------------------------------------
% Description: This script is used for analyzing and visualizing swine stroke data.
% -------------------------------------------------------------------------
% Author:
%   Mark Campbell
%   Carleton University
%   markacampbell@cmail.carleton.ca
% -------------------------------------------------------------------------

% global cmap;

run 'myStartup.m';
maxsz= 0.2; maxh= 2; imgsize= [64 64]; pig= "9.2";
[fmdl, imdl]= mk_pighead_fmdl(maxsz, maxh, imgsize, pig);
% Load data
D= load_HamburgBrain_data(pig);
fn= fieldnames(D);
cd(horzcat('C:\Users\Mark\Documents\GraduateStudies\LAB\HamburgBrain\Figures\', char(pig)));
suffix= ' 20200215';

% configure colormap
% cmap = confg_cmap();

%% 1. SHOW EIT AND PERFUSION DATA WITH PERFUSION ANNOTATIONS
cd(horzcat('C:\Users\Mark\Documents\GraduateStudies\LAB\HamburgBrain\Figures\', char(pig)));
cd seqData;
for i= 1:numel(fn)
    plot_seq_data(D.(fn{i}));
    saveas( gcf, horzcat(D.(fn{i}).name, suffix, '.svg') );
end % end for
close all
cd ../

%% 2. Ensemble average perfusion images
for i = 1:numel(fn)
    show_pre_inj_img(D.(fn{i}));
end

%% 3. Compare average CC for all sequences with adjusted clim
compare_pre_inj_imgs(D, [2,4,6]);










%% VIDEO OF RECONSTRUCTED IMAGE AND BRAIN SEGMENTATION

cd(horzcat('C:\Users\Mark\Documents\GraduateStudies\LAB\HamburgBrain\Figures\', char(pig)));
for i= 4:numel(fn)
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
% noise analysis

% % cross-reference noisy to stim pattern
% % write function that based on measurements removed, return stimulating and
% % measuring pairs involved in that measurement
% 
% vv= real(D.seq1.eit.fdata);
% % remove noisy measurement pairs
% meas_std= std(vv,{},2);
% lim= mean(meas_std)+ std(meas_std); % remove measurements more than 1 sd of mean
% noisy= std(vv,{},2) > lim;
% vv(noisy,:)= 0;
% noisy= reshape(noisy, 32, 32);
% bad_meas= struct;
% meas_pattern= imdl.fwd_model.stimulation.meas_pattern;
% for i= 1:size(noisy, 2) % for each measurement combination
%     was_noisy= noisy(:, i);
%     % the stimulating electrodes for this measurement column were
%     bad_meas(i).stim_elec= find(imdl.fwd_model.stimulation(i).stim_pattern);
%     % the noisy electrodes for these stimulating electrodes were
%     meas_select= reshape(imdl.fwd_model.meas_select, 32, 32); % 32 x 32 logical mtx
%     this_meas_sel= meas_select(:, i)* 1;
%     meas_pattern= imdl.fwd_model.stimulation(i).meas_pattern; % 32 x 32 sparse mtx
%     search_rows= find(was_noisy);
%     % drop dim from 32 to 29 while preserving idx
%     for r= search_rows
%         if this_meas_sel(r)~= 0
%             this_meas_sel(r)= 2;
%         end % end if
%     end % end for
%     new_search_rows= find(this_meas_sel(this_meas_sel>0) == 2);
%     meas_elecs= [];
%     if ~isempty(new_search_rows)
%         for row= new_search_rows'
%             if row > 29
%                 keyboard;
%             else
%                 look_in= meas_pattern(row, :);    
%                 meas_elecs= [meas_elecs; find(look_in)];
%             end % end if
%         end % end for
%     end % end if
%     bad_meas(i).meas_elecs= meas_elecs;
% end % end for
% 
% figure; imagesc(meas_select); axis 'equal'
% figure; imagesc(meas_pattern); axis 'equal'
% %%
% % Stuff to show Andy
% 
% % show raw data
% dd= real(D.seq1.eit.data);
% 
% inspect_eit_elec_and_data(D.seq1, imdl);
% 
% for be=1:32; kk=meas_icov_rm_elecs(imdl, be); ee = find(diag(kk)~=1); plot(dd(ee,'k')'); title(sprintf('bad=%d',be)); pause; end
% 
% % load data
% dd= real(D.seq1.eit.data);
% % plot each electrode and look for worst ones.
% for be=1:32; kk=meas_icov_rm_elecs(imdl, be); ee = find(diag(kk)~=1); plot(dd(ee,:)','k'); title(sprintf('bad=%d',be)); pause; end
% plot(sum(dd(notee,:)))
% % look at weird ones, find which channel they belong to
% channel= find( abs( df(:,2205) - 0.0012848096189744)<1e-10 );
% % find which ellec this belongs to
% for be=1:32; kk=meas_icov_rm_elecs(imdl, be); disp(full([be, kk(channel,channel)])); end
% % plot with this electrode removed
% kk=meas_icov_rm_elecs(imdl, [2,3,7,9, 13, 14, 19, 20, 22, 24, 28, 29, 30]); ee = find(diag(kk)~=1);plot(1:4638,dd','k', 1:4638,dd(ee,:)','r')
% notee=1:544; notee(ee)=[];
% 
% ddm= dd- mean(dd, 2);
% vv= real(D.seq1.eit.data)';
% vvdt= detrend(vv);
% xmax= max(ddm(:));
% xmin= min(ddm(:));
% for i= 1:1024; plot(ddm(i,:)), ylim([xmin xmax]); title(i);pause(1);end
% rsqaured= zeros(1024, 1);
% for i= 1:1024; mdl = fitlm(1:size(ddm, 2),ddm(i,:)); rsqaured(i)=mdl.Rsquared.Adjusted;end
% xax= 1: size(D.seq1.eit.elec_impedance, 2);
% xax1= 1:size(dd, 2);
% plot(xax1, ddm');
% for i= 1:32; plot(xax, real(D.seq2.eit.elec_impedance(i,:))'); title(i); pause(); end
% for i= 1:32; plot(xax, hilbert(real(D.seq2.eit.elec_impedance(i,:))')); title(i); pause(); end
% 
% d= abs(D.seq2.eit.elec_impedance');
% dm= movmean(d, 5);
% ddt= detrend(dm);
% drng= range(ddt);
% plot(xax, ddt);
% bar(drng);

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

function show_pre_inj_img(seq)
    
    allignSyst = false;
    allignDiast1 = false;
    allignDiast2 = true;
    if seq.eit.apn < 1
        start = 1;
    else
        start = seq.eit.apn;
    end % end if
    stop = seq.eit.inj;
    use_lndmrk_idx = find( (seq.eit.peaks < stop) + (seq.eit.peaks >= start) == 2);
    use_syst = seq.eit.peaks(use_lndmrk_idx);
    use_diast = seq.eit.vals(:, use_lndmrk_idx);
    nCycles = length(use_syst);
    
    % set up matrix to view each CC as one row in show_sliecs.
    ccLen = diff(use_diast) + 1;
    z_diast1 =  seq.imgr.elem_data(:, use_diast(2,:));
    z_diast2 =  seq.imgr.elem_data(:, use_diast(1,:));
    z_diast = (z_diast1 + z_diast2) ./ 2;
    ls = use_syst - use_diast(1,:);
    rs = use_diast(2,:) - use_syst;
    if allignSyst
        lpad = abs( ls - max(ls) );
    elseif allignDiast1
        lpad = zeros(nCycles, 1);
    elseif allignDiast2
        lpad = abs(ccLen - max(ccLen));
    end % end if
    longestCC = max(lpad + ccLen);
    z_frames = zeros( size(seq.imgr.elem_data, 1), longestCC * nCycles);
    meanCycle = z_frames(:, 1:longestCC);
    nContributed = zeros(1,longestCC);
    
    % create full image set for each CC that is the difference between each
    % frame and the mean diastole for that CC.
    for i = 1:nCycles
        startIdx = (i - 1) * longestCC + 1 + lpad(i);
        endIdx = startIdx + ccLen(i) - 1;
        imgs = seq.imgr.elem_data( :, use_diast(1,i): use_diast(2,i) ) - z_diast(:,i);
        z_frames(:, startIdx:endIdx) = imgs;
        mcIdxS = 1 + lpad(i);
        mcIdxE = lpad(i) + ccLen(i);
        meanCycle(:, mcIdxS: mcIdxE) = meanCycle(:, mcIdxS: mcIdxE ) + imgs;
        nContributed( mcIdxS: mcIdxE) = nContributed( mcIdxS: mcIdxE) + 1;
    end % end for
    meanCycle = meanCycle  / nContributed;
    z_frames = [z_frames, meanCycle];
    puls_img = seq.imgr;
    puls_img.elem_data = z_frames;
    puls_img.calc_colours.clim = max(meanCycle, [], 'all');
    puls_img.show_slices.img_cols = longestCC;
    figure;
    show_slices(puls_img);
    title(horzcat(seq.name, " Conductivity from individual and ensemble-averaged (bottom) Cardiac Cycles"));

end % end function

function meanCycle = get_mean_cc(seq)

    allignSyst = false;
    allignDiast1 = false;
    allignDiast2 = true;
    if seq.eit.apn < 1
        start = 1;
    else
        start = seq.eit.apn;
    end % end if
    stop = seq.eit.inj;
    use_lndmrk_idx = find( (seq.eit.peaks < stop) + (seq.eit.peaks >= start) == 2);
    use_syst = seq.eit.peaks(use_lndmrk_idx);
    use_diast = seq.eit.vals(:, use_lndmrk_idx);
    nCycles = length(use_syst);
    
    % set up matrix to view each CC as one row in show_sliecs.
    ccLen = diff(use_diast) + 1;
    z_diast1 =  seq.imgr.elem_data(:, use_diast(2,:));
    z_diast2 =  seq.imgr.elem_data(:, use_diast(1,:));
    z_diast = (z_diast1 + z_diast2) ./ 2;
    ls = use_syst - use_diast(1,:);
    rs = use_diast(2,:) - use_syst;
    if allignSyst
        lpad = abs( ls - max(ls) );
    elseif allignDiast1
        lpad = zeros(nCycles, 1);
    elseif allignDiast2
        lpad = abs(ccLen - max(ccLen));
    end % end if
    longestCC = max(lpad + ccLen);
    z_frames = zeros( size(seq.imgr.elem_data, 1), longestCC * nCycles);
    meanCycle = z_frames(:, 1:longestCC);
    nContributed = zeros(1,longestCC);
    
    % create full image set for each CC that is the difference between each
    % frame and the mean diastole for that CC.
    for i = 1:nCycles
        startIdx = (i - 1) * longestCC + 1 + lpad(i);
        endIdx = startIdx + ccLen(i) - 1;
        imgs = seq.imgr.elem_data( :, use_diast(1,i): use_diast(2,i) ) - z_diast(:,i);
        z_frames(:, startIdx:endIdx) = imgs;
        mcIdxS = 1 + lpad(i);
        mcIdxE = lpad(i) + ccLen(i);
        meanCycle(:, mcIdxS: mcIdxE) = meanCycle(:, mcIdxS: mcIdxE ) + imgs;
        nContributed( mcIdxS: mcIdxE) = nContributed( mcIdxS: mcIdxE) + 1;
    end % end for
    
    meanCycle = meanCycle ./ nContributed;
    meanCycle( isnan(meanCycle) ) = 0;

end % end function

% ======================================================================= %

function compare_pre_inj_imgs(D, sel)

    fn = fieldnames(D);
    if nargin == 1
        sel = 1:length(fn);
    end % end if
    MC = struct;
    maxCcLen = 0;
    CcLens = zeros(1, length(sel));
    for j = 1:length(sel)
        i = sel(j);
        MC.(fn{i}) = get_mean_cc( D.(fn{i}) );
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

function plot_av_cc(data, syst, diast)
    
    if size(data, 1) > 1
        data = mean( real(data), 1 );
    end % end if
    
    ls = syst - diast(1,:);
    rs = diast(2,:) - syst;
    lpad = abs( ls - max(ls) );
    rpad = abs( rs - max(rs) );
    tlt_len = max(ls) + max(rs) + 1;
    cc_data = zeros( length(syst), (tlt_len) );
    
    hold on;
    for i = 1:length(syst)
        dd = data( diast(1,i): diast(2,i) );
        dd = dd - min(dd);
        cc_data(i, :) = [zeros(1,lpad(i)), dd, zeros(1, rpad(i))];
    end % end for

    plot(cc_data', 'Color', [0 0 0] + 0.75);
    plot(mean(cc_data, 1), 'k', 'linewidth', 4);

end % end function

% ======================================================================= %

