% TODO VQ image: show as log scale. Move axes closer together.

path_name   = 'C:\Users\Mark\Dropbox\Joaquin\data\';
global fig_path;
global titleFS;
global axsFS;
global textFS;
global doPlotting;
global mkPaperFig;
global roi;
titleFS     = 16;
axsFS       = 16;
textFS      = 12;
fig_path    = 'C:\Users\Mark\Dropbox\Joaquin\figures\';
doPlotting  = false;
mkPaperFig  = true;

roi = horse_roi();
% for horse = 1:4
for horse = 2
    close all;
    switch horse
        case 1
            horseID = 296955;
        case 2
            horseID = 296439;
        case 3
            horseID = 296695;
        case 4
            horseID = 296954;
    end
    do_analysis(path_name, horseID);
end


function do_analysis(path_name, horseID)
    global roi;
    global doPlotting;
    global mkPaperFig;
    EIT = load_joaquin_data(path_name, horseID);	% load raw and filtered data
    fn = fieldnames(EIT);
    for i = 1:length(fn)
        if contains(lower(fn{i}), 'baseline')
            [fmdl, imdl] = mk_horse_model(EIT.info.elecShift);    % models fmdl and imdl
            eitFile = EIT.(fn{i});
            msel = find(fmdl.meas_select);
            vhCommon = mean(eitFile.fdata(msel,:), 2);
            break
        end
    end
    % Noise compensation
    specs = struct; specs.n = 90; specs.type = 'meas';
    [worst, scores, ~, ~] = worst_n_elecs(EIT, imdl, specs);
    rmMeas = worst(scores >= 0.25);
    imdl_comp = comp_RM_bad_elec(imdl, rmMeas, 'meas');
    fprintf('\nmeasurements removed: %.0f \n', length(rmMeas));

    % find clim
    CLIM = -inf;
    for i = 1:length(fn)
        eitFile = EIT.(fn{i});
        if ~isfield(eitFile, 'data')
            continue
        end
        eitFile.data = eitFile.data(msel,:);
        eitFile.fdata = eitFile.fdata(msel,:);
        EIT.(fn{i}).triads = horse_breath_finder(eitFile.fdata);
        EIT.(fn{i}).imgr = inv_solve(imdl_comp, vhCommon, eitFile.fdata);
        clim = max(EIT.(fn{i}).imgr.elem_data,[],'all');
        CLIM = max([clim, CLIM]);
    end

    for i = 1:length(fn)
        eitFile = EIT.(fn{i});
        if ~isfield(eitFile, 'data')
            continue
        end
        eitFile.imgr.calc_colours.ref_level = 0;
        eitFile.imgr.calc_colours.clim = CLIM;
        EIT.(fn{i}) = eitFile;
    end
    
    mk_paper_fig(EIT, horseID);
%     plot_ventilation_together(EIT, horseID);
end % end function


function plot_ventilation_together(EIT, horseID)
    global roi;
    fig_path = 'C:\Users\Mark\Dropbox\Joaquin\figures\files_shown_consecutively\';
    x1 = 1;
    bigFig();
    switch horseID
        case 296955 % 1
            fn = {'baseline','IE11_T5_min','IE11_T10_min','IE11_T25_min','IE13_T15_min','IE13_T25_min'};
        case 296439 % 2
            fn = {'baseline','IE13_T15_min','IE13_T25_min','IE11_T15_min','IE11_T25_min'};
        case 296695 % 3
            fn = {'Baseline','IE_11_T15_min','IE11_T25_min','IE13_T15_min','IE13_T25_min'};
        case 296954 % 4
            fn = {'baseline','IE13_T15_min','IE13_T25_min','IE11_T15_min','IE11_T25_min'};
    end
    for i=1:length(fn)
        imgs = calc_slices(EIT.(fn{i}).imgr);
        for j = 1:length(imgs)
            thisImg = imgs(:,:,j);
            thisImg(roi.BothLungs == 0) = 0;
            imgs(:,:,j) = thisImg;
        end
        
        if contains(fn{i}, '_')
            splitFn = strsplit(fn{i}, 'T');
            time = strsplit(splitFn{2}, '_');
            time = str2double(time{1});
            addFrames = time * 60 * EIT.(fn{i}).fs;
            x1 = baselineEnd + addFrames;
            if time == 25
                baselineEnd = x1 + length(imgs) - 1; % e.g. IE_11 after IE_13 time is rel to IE 13 T25
            end
        end
        x2 = x1 + length(imgs) - 1;
        if contains(fn{i}, 'Baseline') || contains(fn{i}, 'baseline')
           baselineEnd = x2;
        end
        xax = (x1:x2) / EIT.(fn{i}).fs / 60;
        plot(xax, -squeeze(sum(imgs, [1,2], 'omitnan')));
        text( x1 / EIT.(fn{i}).fs / 60, 0, remove_underscores(fn{i}));
        hold on;
        xline(x1 / EIT.(fn{i}).fs / 60);
        x1 = x1 + length(imgs);
    end
    sgtitle(sprintf('Total lung Impedance for %.0f ventilation files', horseID));
    xlabel('Time (minutes)');
    ylabel('Total lung impedance');
    printPDF(sprintf('%s%0.f Time Course.pdf', fig_path, horseID));
end


function mk_paper_fig(EIT, horseID)
    switch horseID
        
        case 296955 % 1
            breathSel = 2;
            perfSel = 3;
            %IE11 breathing
            ventFile = EIT.IE_11_Perfusion_breathing;
            perfFile = EIT.IE_11_Perfusion_breathing;
            do_the_figure(ventFile, perfFile, 2, 4);
            %IE13 breathing
            ventFile = EIT.IE13_Perfusion_breathing;
            perfFile = EIT.IE13_Perfusion_breathing;
            do_the_figure(ventFile, perfFile, 1, 4);
            %IE11 apnea
            ventFile = EIT.IE_11_Perfusion_breathing;
            perfFile = EIT.IE_11_perfusion_apnea;
            do_the_figure(ventFile, perfFile, 2, 3);
            %IE13 apnea
            ventFile = EIT.IE13_Perfusion_breathing;
            perfFile = EIT.IE_13_Perfusion_apnea;
            do_the_figure(ventFile, perfFile, 1, 4);
            
        case 296439 % 2
            perfSel = 3;
            %IE11 breathing
            ventFile = EIT.IE11_Perfusion_breathing;
            perfFile = EIT.IE11_Perfusion_breathing;
            do_the_figure(ventFile, perfFile, 1, perfSel);
            %IE13 breathing
            ventFile = EIT.IE13_Perfusion_breathing;
            perfFile = EIT.IE13_Perfusion_breathing;
            do_the_figure(ventFile, perfFile, 3, perfSel);
            %IE11 apnea
            ventFile = EIT.IE11_Perfusion_apnea;
            perfFile = EIT.IE11_Perfusion_apnea;
            do_the_figure(ventFile, perfFile, 1, perfSel);
            %IE13 apnea
            ventFile = EIT.IE13_Perfusion_apnea;
            perfFile = EIT.IE13_Perfusion_apnea;
            do_the_figure(ventFile, perfFile, 1, perfSel);
            
        case 296695 % 3
            breathSel = 2;
            perfSel = 3;
            %IE11 breathing
            ventFile = EIT.IE11_Perfusion_breathing;
            perfFile = EIT.IE11_Perfusion_breathing;
            do_the_figure(ventFile, perfFile, 3, 4);
            %IE13 breathing
            ventFile = EIT.IE13_Perfusion_breathing;
            perfFile = EIT.IE13_Perfusion_breathing;
            do_the_figure(ventFile, perfFile, 1, 3);
            %IE11 apnea
            ventFile = EIT.IE11_Perfusion_apnea;
            perfFile = EIT.IE11_Perfusion_apnea;
            do_the_figure(ventFile, perfFile, 1, 4);
            %IE13 apnea
            ventFile = EIT.IE13_Perfusion_breathing;
            perfFile = EIT.IE13_Perfusion_Apnea;
            do_the_figure(ventFile, perfFile, 1, 3);
        
        case 296954 % 4
            breathSel = 2;
            perfSel = 3;
            %IE11 breathing
            ventFile = EIT.IE11_Perfusion_breathing;
            perfFile = EIT.IE11_Perfusion_breathing;
            do_the_figure(ventFile, perfFile, 2, 2);
            %IE13 breathing
            ventFile = EIT.IE13_Perfusion_breathing;
            perfFile = EIT.IE13_Perfusion_breathing;
            do_the_figure(ventFile, perfFile, 1, 3);
            %IE11 apnea
            ventFile = EIT.IE11_Perfusion_breathing;
            perfFile = EIT.IE11_Perfusion_apnea;
            do_the_figure(ventFile, perfFile, 2, 2);
            %IE13 apnea
            ventFile = EIT.IE_13_Perfusion_apnea;
            perfFile = EIT.IE_13_Perfusion_apnea;
            do_the_figure(ventFile, perfFile, 2, 4);
    end
end


function do_the_figure(ventFile, perfFile, breathSel, perfSel)
    fig_path = 'C:\Users\Mark\Dropbox\Joaquin\figures\VQ_ratio\';
    global ax2;
    global ax4;
    global ax6;
    
    even = 128;
    whiteGrad = linspace(0,1,254)'; % cmap
    colour = ones(254,1);            % cmap
    colour(1:even) = linspace(0, 1, even);
    ventCmap = [whiteGrad, whiteGrad, colour];
    ventCmap(1,:) = 0.9;
    ventCmap(2,:) = 0;
    perfCmap = [colour, whiteGrad, whiteGrad];
    perfCmap(1,:) = 0.9;
    perfCmap(2,:) = 0;
    
    fg = bigFig(3);
    clf();
    imgs11 = calc_slices(ventFile.imgr);
    triads11 = ventFile.triads;
    ventDist11 = ventilation_distribution(imgs11, ventFile.triads);
    close;
    perfDist11 = perfusion_distribution(perfFile);
    
    mk_vent_and_perf_fig(imgs11, triads11(breathSel, :), perfFile, perfSel);
    vent_row_fig(ventDist11, breathSel);
    perf_row_fig(perfDist11, perfSel);
    v_by_q_fig();
    v_by_q_row_fig();
    
    ax2.Colormap = ventCmap;
    ax4.Colormap = perfCmap;
    ax6.Colormap = flipud(ax6.Colormap);
    ax6.Colormap(1,:) = 0.9;
    ax6.Colormap(2,:) = 0;
    sgtitle(remove_underscores( perfFile.name) );
    printPDF( sprintf('%s%s VQ fig', fig_path, perfFile.name) );
end


function mk_vent_and_perf_fig(imgs, triad, perfFile, perfSel)
    global roi;
    global ax2;
    global ax4;
    global ventLungZSum1;
    global perfLungZSum1;
%     global lim;
    MAXLIM = 1;
    bothLungs = roi.BothLungs;
    lungRows = find(sum(bothLungs, 2) > 0);
    lungCols = find(sum(bothLungs, 1) > 0);
    nPixelsInLungs = sum(bothLungs, 'all');
    even = 128;
    
    background = zeros(32,32);
    background(bothLungs == 0) = -126;
    background(roi.Inside == 0) = -127;
    
    breathZ = calc_breath_delta_z(-imgs, triad);  % breath difference image
    ventLungZ = breathZ(bothLungs == 1);
    ventLungZ(ventLungZ <= 0) = nan;
    totalZ = sum(ventLungZ, 'all', 'omitnan');
    unPixVal = totalZ / sum(ventLungZ > 0, 'all', 'omitnan');
    ventLungZSum1 = (ventLungZ ./ unPixVal);
    ventLungZUniform = ventLungZSum1;
    ventLungZUniform(ventLungZUniform > 0) = ventLungZUniform(ventLungZUniform > 0) - 1; % uniform distributed pixel has value 0
    
    ventImg = zeros(32,32);
    ventImg(bothLungs == 1) = ventLungZUniform;
    upperBound = max(ventImg, [], 'all');
    lowerBound = min(ventImg, [], 'all');
    lim = max(abs([upperBound, lowerBound]));

    % PERF INIT
    NPOINTS = 5;
    refXVal = 47;
    start = perfFile.perfStart - refXVal;
    stop = perfFile.perfEnd;
    nFrames = stop  -start + 1;
    timePoints = linspace(refXVal, nFrames, NPOINTS + 2);
    timePoints = round(timePoints(2: NPOINTS+1));
    
    perfSel = timePoints(perfSel);
    
    if strcmp(perfFile.name, '296439 IE11 Perfusion breathing')
        perfSel = 1380 - perfFile.perfStart + 1;
    elseif strcmp(perfFile.name, '296439 IE13 Perfusion breathing')
        perfSel = 2888 - perfFile.perfStart + 1;
    end
    
    perfImgr = perfFile.imgr;
    perfImgr.elem_data = perfFile.imgr.elem_data(:, start: stop);
    clim = max(perfImgr.elem_data,[],'all');
    perfImgr.calc_colours.ref_level = 0;
    perfImgr.calc_colours.clim = clim;
    perfImgs = calc_slices(perfImgr);           % get perfusion images
    
%     rsperfImgs = reshape(perfImgs, [32^2, 5651]);
    perfZ = perfImgs(:,:,perfSel) - perfImgs(:,:,1); % select specific perfusion image
    perfLungZ = perfZ(bothLungs == 1); % only lungs
    
    perfLungZ = perfZ(bothLungs == 1); % only lungs
    perfLungZ(perfLungZ <= 0) = nan; % zero pixels less than 0
    totalZ = sum(perfLungZ,'all', 'omitnan'); % find total
    unPixVal = totalZ / sum(perfLungZ > 0, 'all', 'omitnan'); % pixel as proportion of total
    perfLungZSum1 = (perfLungZ ./ unPixVal); 
    perfLungZUniform = perfLungZSum1;
    perfLungZUniform(perfLungZUniform > 0) = perfLungZUniform(perfLungZUniform > 0) - 1; % uniform distributed pixel has value 0

    perfImg = zeros(32, 32);
    perfImg(bothLungs == 1) = perfLungZUniform;
    
    upperBound = max(perfImg,[],'all');
    lowerBound = min(perfImg,[],'all');
    thisLim = max(abs([upperBound,lowerBound]));
    if thisLim > lim
       lim = thisLim; 
    end
    
    if lim > MAXLIM
        lim = MAXLIM;
        addGTLT = true;
    else
        addGTLT = false;
    end
    
    ventImg(ventImg > lim) = lim;
    ventImg(ventImg < -lim) = -lim;
    perfImg(perfImg > lim) = lim;
    perfImg(perfImg < -lim) = -lim;
    
    % MAKE VENT IMAGE
    ventImg(1) = lim; ventImg(2) = -lim;
    ax2 = subplot( 3, 2, 1 );
    slc = show_slices(ventImg);
    slc(isnan(ventImg)) = 2;
    slc(1) = even; slc(2) = even;
    slc = slc + background;

    [X, Y] = meshgrid(1:size(slc, 2), 1:size(slc, 1));
    [X2, Y2] = meshgrid(1:0.25:size(slc, 2), 1:0.25:size(slc, 1));
    betterSlc = interp2(X, Y, slc, X2, Y2, 'linear');
    betterbkg = interp2(X, Y, background, X2, Y2, 'linear');
    betterbkg(betterbkg ~= -127) = 0;
    image(betterSlc + betterbkg);
    axis image; axis off; axis equal; axis tight;
    
    lowTick = 2;
    highTick = 253;
    highLbl = horzcat(sprintf('%.0f', lim * 100 + 100),'%');
    lowLbl = horzcat(sprintf('%.0f', (0.5 / lim) * 100),'%');
    if addGTLT
        lowLbl = horzcat('<= ', lowLbl);
        highLbl = horzcat('>= ', highLbl);
    end
    
    colorbar;
    fg = gcf();
    fg.Children(1).Ticks = [lowTick, even, highTick];
    fg.Children(1).TickLabels = {lowLbl,'100%',highLbl};
    title('Ventilation image');
    
    % MAKE PERF IMAGE
    perfImg(1) = lim; perfImg(2) = -lim;
    ax4 = subplot(3, 2, 3);
    
    slc = show_slices(perfImg);
    slc(1) = even; slc(2) = even;
    slc(isnan(perfImg)) = 2;
    slc = slc + background;

    [X, Y] = meshgrid(1:size(slc,2), 1:size(slc,1));
    [X2, Y2] = meshgrid(1:0.25:size(slc,2), 1:0.25:size(slc,1));
    betterSlc = interp2(X, Y, slc, X2, Y2, 'linear');
    betterbkg = interp2(X, Y, background, X2, Y2, 'linear');
    betterbkg(betterbkg ~= -127) = 0;
    image(betterSlc + betterbkg);
    axis image; axis off; axis equal; axis tight;
    
    colorbar();
    fg = gcf();
    fg.Children(1).Ticks = [lowTick, even, highTick];
    fg.Children(1).TickLabels = {lowLbl, '100%', highLbl};
    title('Perfusion image');
end


function vent_row_fig(ventDist, breathNumber)
    global ax3;
    global ie11RowsV;
    global roi;
    lungRows = sum(roi.BothLungs,2) ~= 0;
    numRows = length(lungRows);
    lungRowIdx = find(lungRows);
    ax3 = subplot( 3, 2, 2 );
    
    ie11RowsV = nan(numRows, 1);
    for i = 1:length(lungRowIdx)
        idx = lungRowIdx(i);
        dat = ventDist.(sprintf('rowMeanRes%.0f',i));
        ie11RowsV(idx) = dat(breathNumber);
    end
    
    unRowVal = 1 / sum(ie11RowsV > 0);
    ie11RowsV(ie11RowsV == 0) = nan;
    ie11RowsV2 = ie11RowsV;
    ie11RowsV2 = ie11RowsV2 ./ unRowVal; % 1 is uniform
    plot(flipud(ie11RowsV2), 1:numRows, 'LineWidth', 2);
    hold on;
    plot([1, 1],[1, numRows], '--k');
    ax3.YAxis.Limits = [1,numRows];
    ax3.YAxis.TickValues = [3, numRows - 2];
    ax3.YAxis.TickLabels = {'Dorsal','Ventral'};
    title('Relative ventilation per lung row');
    xlabel('Percentage of uniform ventilation');
end


function perf_row_fig(perfDist, number)
    global ax3;
    global ie11RowsP;
    global roi;
    lungRows = sum(roi.BothLungs,2) ~= 0;
    numRows = length(lungRows);
    lungRowIdx = find(lungRows);
    ie11RowsP = nan(numRows, 1);
    for i = 1 : length(lungRowIdx)
        idx = lungRowIdx(i);
        dat = perfDist.(sprintf('rowMeanRes%.0f',i));
        ie11RowsP(idx) = dat(number);
    end
    unRowVal = 1 / sum(ie11RowsP > 0);
    ie11RowsP(ie11RowsP == 0) = nan;
    ie11RowsP2 = ie11RowsP;
    ie11RowsP2 = ie11RowsP2 ./ unRowVal; % 1 is uniform (for ratio)
    ax5 = subplot( 3, 2, 4);
    plot(flipud(ie11RowsP2), 1:numRows, 'r', 'LineWidth', 2);
    hold on;
    plot([1, 1],[1, numRows], '--k');
    
    ax5.YAxis.Limits = [1,numRows];
    ax5.YAxis.TickValues = [3, numRows - 2];
    ax5.YAxis.TickLabels = {'Dorsal','Ventral'};
    
    Mlim = max([ax5.XAxis.Limits(2), ax3.XAxis.Limits(2)]);
    mlim = min([ax5.XAxis.Limits(1), ax3.XAxis.Limits(1)]);
    lim = max(abs([mlim, Mlim] - 1));
    ax3.XAxis.Limits(1) = 1 - lim;
    ax5.XAxis.Limits(1) = 1 - lim;
    ax3.XAxis.Limits(2) = 1 + lim;
    ax5.XAxis.Limits(2) = 1 + lim;
    
    xTickVals = ax3.XAxis.TickValues;
    ax3.XAxis.TickLabels = xTickVals * 100;
    xTickVals = ax5.XAxis.TickValues;
    ax5.XAxis.TickLabels = xTickVals * 100;
    
    title('Relative perfusion per lung row');
    xlabel('Percentage of uniform perfusion');
end


function v_by_q_fig()
    global roi;
    global ventLungZSum1;
    global perfLungZSum1;
    global ax6;
    MAXLIM = 2;
    bothLungs = roi.BothLungs;
    nPixelsInLungs = sum(bothLungs, 'all');
    unPixVal = 1;
    even = 128; 
    bothLungs = roi.BothLungs;
    background = zeros(32,32);
    background(bothLungs == 0) = -126;
    background(roi.Inside == 0) = -127;
        
    vqImg = zeros(32, 32);                       % instantiating to avoid div by 0 of V/Q
    vqImg(bothLungs == 1) = (ventLungZSum1 ./ perfLungZSum1); % subtract 1 to make 0 uniform
    vqImg = log(vqImg);
    vqImg(isinf(vqImg)) = nan;
    upperBound = max(vqImg, [], 'all');
    lowerBound = min(vqImg, [], 'all');
    lim = max(abs([upperBound, lowerBound]));
    if lim > MAXLIM
        lim = MAXLIM;
        addGTLT = true;
    else
        addGTLT = false;
    end
    
    vqImg(vqImg > lim) = lim;
    vqImg(vqImg < -lim) = -lim;
    vqImg(1) = lim;
    vqImg(2) = -lim;
    
    ax6 = subplot( 3, 2, 5);
    slc = show_slices(vqImg);
    slc(1) = even; slc(2) = even;
%     slc = slc(lungRows(1):lungRows(end), lungCols(1):lungCols(end));
    slc(slc == 1) = 2;    
    [X, Y] = meshgrid(1:size(slc,2), 1:size(slc,1));
    [X2, Y2] = meshgrid(1:0.1:size(slc,2), 1:0.1:size(slc,1));
    betterSlc = interp2(X, Y, slc, X2, Y2, 'linear');
    betterbkg = interp2(X, Y, background, X2, Y2, 'linear');
    betterbkg(betterbkg~=-127) = 0;
    
    image(betterSlc + betterbkg);
    axis image;axis off;axis equal;axis tight;
    
    highTick = 253;
    lowTick = 2;
    highLbl = num2str(lim);
    lowLbl = num2str(-lim);
    if addGTLT
        lowLbl = horzcat('<= ', lowLbl);
        highLbl = horzcat('>= ', highLbl);
    end
    colorbar();
    fg = gcf();
    fg.Children(1).Ticks = [lowTick, even, highTick];
    fg.Children(1).TickLabels = {lowLbl, 'uniform', highLbl};
    title('Log V/Q');
end


function v_by_q_row_fig()
    global ie11RowsV;
    global ie11RowsP;
    ax7 = subplot( 3, 2, 6);
    numRows = length(ie11RowsV);
    vbyQDat = ie11RowsV ./ ie11RowsP;
    plot(flipud(vbyQDat), 1:32, 'k', 'LineWidth', 2); hold on;
    plot([1, 1],[1, numRows], '--k');
    lim = max(abs(1-vbyQDat));
    ax7.XAxis.Limits = [1 - lim, 1 + lim];
    ax7.YAxis.TickValues = [3, numRows-2];
    ax7.YAxis.TickLabels = {'Dorsal','Ventral'};
    xlabel('V/Q ratio');
    title('Ventilation / perfusion ratio per lung row');
end


function global_signal()
    % GLOBAL SIGNAL
    ax1 = subplot( 5, 2, 1:4 );
    xax = (1:length(brth11.fdata)) ./ brth11.fs;
    elemSums = squeeze(-sum(imgs11, [1,2], 'omitnan') ); 
    plot(xax, elemSums, 'LineWidth', LW);
    % show perf curve
    hold on;
    endEx = [triads11(:,1)', triads11(end,3), length(xax)];
    perfLine = csapi(xax(endEx), elemSums(endEx), xax(endEx(1):endEx(end)));
    plot( xax(endEx(1):endEx(end)) , perfLine, 'r', 'LineWidth', LW);
    xlabel('Time (s)'); ylabel('Global signal (AU)');
    
    for i=1:size(imgs11, 3)
        temp = imgs11(:,:,i);
        temp(bothLungs==0) = 0;
        imgs11(:,:,i) = temp;
    end
end


function fig = bigFig(num)
    % make a full-screen figure
    if nargin == 0
        fig=figure();
    else
        fig=figure(num);
    end
    fig.Units='normalized';
    fig.OuterPosition=[0 0 1 1];
    clf(); 
    fig.PaperOrientation='landscape';
end