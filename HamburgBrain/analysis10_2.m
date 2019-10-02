% -------------------------------------------------------------------------
% Description: This script is used for analyzing and visualizing swine stroke data.
% -------------------------------------------------------------------------
% Author:
%   Mark Campbell
%   Carleton University
%   markacampbell@cmail.carleton.ca
% -------------------------------------------------------------------------

run 'myStartup.m';
maxsz= 0.2; maxh= 2; imgsize= [64 64];
[fmdl, imdl]= mk_pighead_fmdl(maxsz, maxh, imgsize);
% Load data
pig= "10.2";
D= load_HamburgBrain_data(pig);
fn= fieldnames(D);
cd 'C:\Users\Mark\Documents\GraduateStudies\LAB\HamburgBrain\Figures\10.2';
suffix= ' rm_1_13_14_18_19_23_28';
%% VIDEO OF RECONSTRUCTED IMAGE AND BRAIN SEGMENTATION

cd 'C:\Users\Mark\Documents\GraduateStudies\LAB\HamburgBrain\Figures\10.2';
for i= 1:numel(fn)
    start= D.(fn{i}).eit.inj;
    stop= D.(fn{i}).eit.inj+ 1000;
    mk_vid(D.(fn{i}), start, stop, suffix);
end % end for
%% EIT, PREFUSION, AND RECONSTRUCTED IMAGE ENSEMBLES IN ONE FIGURE
cd 'C:\Users\Mark\Documents\GraduateStudies\LAB\HamburgBrain\Figures\10.2';
opt.pv= 3;
opt.usefData= 1;
opt.plotLM= 2;
opt.ensemble= 'one';
cd ensemble;
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
ttp2= [];
for i= [1, 3, 5]
%     opt.start= 50;
%     opt.stop= D.(fn{i}).eit.inj;
    opt.start= D.(fn{i}).eit.inj;
    opt.stop= D.(fn{i}).eit.inj+ 1000;
    ttp1= [ttp1, get_ttp(D.(fn{i}), 2, opt)];
%     title("time to peak for " + D.(fn{i}).name);
end % end for
for i= [2, 4, 6]
%     opt.start= 50;
%     opt.stop= D.(fn{i}).eit.inj;
    opt.start= D.(fn{i}).eit.inj;
    opt.stop= D.(fn{i}).eit.inj+ 1000;
    ttp2= [ttp2, get_ttp(D.(fn{i}), 2, opt)];
%     title("time to peak for " + D.(fn{i}).name);
end % end for
figure; imagesc(ttp1); axis 'equal'; colorbar;
figure; imagesc(ttp2); axis 'equal'; colorbar;

%% SHOW MEAN FRAME
cd meanFrame
for i= 1:numel(fn)
    meanFrame = mean(calc_slices(D.(fn{i}).imgr), 3);
    figure('units','normalized','outerposition',[0 0 1 1]);
    imagesc(meanFrame);
    title(D.(fn{i}).name+ suffix);
    colorbar;
    axis equal
    print_convert(char(D.(fn{i}).name+ suffix+ ".png"));
end % end for
close all
cd ../

%% SHOW EIT AND PERFUSION DATA WITH PERFUSION ANNOTATIONS
cd 'C:\Users\Mark\Documents\GraduateStudies\LAB\HamburgBrain\Figures\10.2';
cd seqData
opt.pv= 3;
opt.usefdata= 1;
for i= 1:numel(fn)
    plot_seq_data(D.(fn{i}), opt);
    print_convert(char(D.(fn{i}).name+ suffix+ ".png"));
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
% % Plot ensemble average of pixels in brain segmentation
% opt.plotLM= 2;
% for i= 1:numel(fn)
%     for j= 0:4
%         opt.section= j;
%         plot_ensemble_per_pixel(D.(fn{i}), opt);
%     end % end for
%     close all
% end % end for



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

% cross-reference noisy to stim pattern
% write function that based on measurements removed, return stimulating and
% measuring pairs involved in that measurement

vv= real(D.seq1.eit.fdata);
% remove noisy measurement pairs
meas_std= std(vv,{},2);
lim= mean(meas_std)+ std(meas_std); % remove measurements more than 1 sd of mean
noisy= std(vv,{},2) > lim;
vv(noisy,:)= 0;
noisy= reshape(noisy, 32, 32);
bad_meas= struct;
meas_pattern= imdl.fwd_model.stimulation.meas_pattern;
for i= 1:size(noisy, 2) % for each measurement combination
    was_noisy= noisy(:, i);
    % the stimulating electrodes for this measurement column were
    bad_meas(i).stim_elec= find(imdl.fwd_model.stimulation(i).stim_pattern);
    % the noisy electrodes for these stimulating electrodes were
    meas_select= reshape(imdl.fwd_model.meas_select, 32, 32); % 32 x 32 logical mtx
    this_meas_sel= meas_select(:, i)* 1;
    meas_pattern= imdl.fwd_model.stimulation(i).meas_pattern; % 32 x 32 sparse mtx
    search_rows= find(was_noisy);
    % drop dim from 32 to 29 while preserving idx
    for r= search_rows
        if this_meas_sel(r)~= 0
            this_meas_sel(r)= 2;
        end % end if
    end % end for
    new_search_rows= find(this_meas_sel(this_meas_sel>0) == 2);
    meas_elecs= [];
    if ~isempty(new_search_rows)
        for row= new_search_rows'
            if row > 29
                keyboard;
            else
                look_in= meas_pattern(row, :);    
                meas_elecs= [meas_elecs; find(look_in)];
            end % end if
        end % end for
    end % end if
    bad_meas(i).meas_elecs= meas_elecs;
end % end for

figure; imagesc(meas_select); axis 'equal'
figure; imagesc(meas_pattern); axis 'equal'
%%
% Stuff to show Andy

% show raw data
D= load_HamburgBrain_data(pig);
xax= 1: size(vv, 1);
vv= D.seq1.eit.data;

inspect_eit_elec_and_data(D.seq1, imdl);
