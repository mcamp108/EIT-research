function show_ensemble(seq, opt)
% -------------------------------------------------------------------------
% DESCRIPTION:
%   show_ensemble(seq, opt)
%
%   Plots ensembles.
% -------------------------------------------------------------------------
% PARAMETERS:
%   seq:        
%       A swine sequence struct from load_HamburgBrain_data.
%
%   opt:    
%       An options struct with possible fields:
%
%       plotLM:
%           1: Use perfusion peaks as window start reference
%           2: Use perfusion valleys as window star reference (default)
%
%       section:
%           0: Ensemble for Brain for Full Timeseries (default)
%           1: Ensemble for Brain Before Apnoea
%           2: Ensemble for Brain After Apnoea and Before Injection
%           3: Ensemble for Brain After Injection and Before Apnoea End
%           4: Ensemble for Brain After Apnoea End
%
%       imgCols: 
%           Same as img_cols input in show_slices (default = 10).
%
%       ensemble:
%           'each': show the mean of each ensemble as one frame
%           'one':  show the mean of all ensembles
% -------------------------------------------------------------------------
% RETURNS:
%   Figure showing specified ensembles.
% -------------------------------------------------------------------------  
% AUTHOR:
%   Mark Campbell
%   Carleton University
%   markacampbell@cmail.carleton.ca
% -------------------------------------------------------------------------

% Check inputs and options
imgr= seq.imgr;

% which landmarks to use
if ~exist('opt', 'var') || ~isfield(opt, 'plotLM')
    opt.plotLM= 2;
    opt.imgCols= 10;
end % end if

% which time range to plot
if ~isfield(opt, 'section')
    opt.section= 0;
end % end if

% number of reconstructed images to show in figure
if ~isfield(opt, 'imgCols') || opt.imgCols< 1
    imgCols= 10;
elseif opt.imgCols> ensembleLength
    imgCols= ensembleLength;
else
    imgCols= opt.imgCols;
end % end if

% which dimension of ensemble to average
if ~isfield(opt, 'ensemble') || strcmp(opt.ensemble, 'each')
    take_mean_of_dim= 3; % set to 'each'
elseif strcmp(opt.ensemble, 'one')
    take_mean_of_dim= 2;
else
    disp("Unrecognized ensemble average dimension");
end % end if

% Which segment of time series at which to look
if opt.section== 1
    start= 1;
    stop= seq.eit.apn;
    figureTitle= seq.name+ " 1 - Ensemble for Brain Before Apnoea";
elseif opt.section== 2
    start= seq.eit.apn;
    stop= seq.eit.inj;
    figureTitle= seq.name+ " 2 - Ensemble for Brain After Apnoea and Before Injection";
elseif opt.section== 3
    start= seq.eit.inj;
    stop= seq.eit.vnt;
    figureTitle= seq.name+ " 3 - Ensemble for Brain After Injection and Before Apnoea End";
elseif opt.section== 4
    start= seq.eit.vnt;
    stop= size(imgr.elem_data, 2);
    figureTitle= seq.name+ " 4 - Ensemble for Brain After Apnoea End";
else
    start= 1;
    stop= size(imgr.elem_data, 2);
    figureTitle= seq.name+ " - Ensemble for Brain for Full Timeseries";
end % end if

% get ensembles for EIT, perfusion, and reconstructed images
e_opt.sel= [1 1 1];
e_opt.start= start;
e_opt.stop= stop;
ensemble= get_ensembles(seq, opt.plotLM, e_opt);
eit_e= ensemble.eit_ensemble;
perf_e= ensemble.perf_ensemble;
imgr_e= ensemble.imgr_ensemble;

% average ensembles across selected dimension
eit_e= squeeze(mean(eit_e, take_mean_of_dim));
perf_e= squeeze(mean(perf_e, take_mean_of_dim));
imgr_e= squeeze(mean(imgr_e, take_mean_of_dim));

% partition out brain segmentation from reconstructed image ensemble
brain_seg= get_brain_segmentation_for_pig(seq);
brain_full= imgr_e(brain_seg.idx, :);
brain_topL= imgr_e(brain_seg.topL, :);
brain_topR= imgr_e(brain_seg.topR, :);
brain_botL= imgr_e(brain_seg.botL, :);
brain_botR= imgr_e(brain_seg.botR, :);

figure('units','normalized','outerposition',[0 0 1 1]);
% plot slices
subplot(4, 1, 1);
    sgtitle(figureTitle, 'FontSize', 20);
    imgc= seq.imgr;
    imgc.show_slices.img_cols = imgCols;
    use_frames= round(linspace(1, size(imgr_e, 2), imgCols));
    imgc.elem_data= imgr_e(:, use_frames);
    clim= max(imgc.elem_data(:));
    imgc.calc_colours.ref_level= 0;
    imgc.calc_colours.clim= clim;
    show_slices(imgc);
    colorbar;
    title("Reconstructed Frames",'FontSize', 12);

% Plot EIT TBV data alongside
subplot(4, 1, 2);
    eit_data= sum(eit_e, 1);
    maxEit= max(eit_data);
    minEit= min(eit_data);
    hold on
    plot(eit_data, 'b');
    
    for i= use_frames
        plot([i, i], [minEit, maxEit], 'r');
    end % end for
    hold off
    title("Ensemble EIT Signal",'FontSize', 12);
    xlabel("Frame");
    ylabel("Voltage (mV)");

% Plot brain signal
subplot(4, 1, 3);
    plot(sum(brain_full, 1), 'g');
    hold on;
    plot(sum(brain_topL, 1));
    plot(sum(brain_topR, 1));
    plot(sum(brain_botL, 1));
    plot(sum(brain_botR, 1));
    title("Ensemble EIT Brain Signal",'FontSize', 12);
    legend('Whole Brain', 'Top Left', 'Top Right', 'Bottom Left', 'Bottom Right');
    xlabel("Frame");
    ylabel("Conductivity");
    hold off;

subplot(4, 1, 4);
plot(perf_e);
title("Ensemble Perfusion Signal",'FontSize', 12);
xlabel("Frame");
ylabel("Arterial Pressure (mmHg)");

end % end function