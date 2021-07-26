% -------------------------------------------------------------------------
% DESCRIPTION:
%   Main analysis file for Joaquin 
% -------------------------------------------------------------------------

path_name = 'C:\Users\Mark\Dropbox\Joaquin\data\';
global fig_path;
global titleFS;
global axsFS;
global textFS;
global doPlotting
titleFS = 16;
axsFS = 16;
textFS = 12;
fig_path = 'C:\Users\Mark\Dropbox\Joaquin\figures\_EILI_before_inspiratory_pause\';
doPlotting = false;

calc_colours('cmap_type', 'blue_red');
colormap(calc_colours('colourmap'));
for horse = 1 : 4
% for horse=1
% for horse = 2
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
    global doPlotting
    EIT = load_joaquin_data(path_name, horseID);	% load raw and filtered data
    roi = horse_roi();
%     fn = fieldnames(EIT);
    switch horseID
        case 296955 % 1
            fn = {'baseline','IE11_T5_min','IE11_T10_min','IE11_T25_min','IE13_T15_min','IE13_T25_min'};
%             fn = {'IE11_T5_min','IE11_T10_min','IE11_T25_min','IE13_T15_min','IE13_T25_min'};
            refFile = 'IE11_T5_min';
        case 296439 % 2
            fn = {'baseline','IE13_T15_min','IE13_T25_min','IE11_T15_min','IE11_T25_min'};
%             fn = {'IE13_T15_min','IE13_T25_min','IE11_T15_min','IE11_T25_min'};

%             fn = {'IE13_T15_min', 'IE11_T15_min'};
%             fn = {'IE11_T15_min', 'IE13_T15_min'};
            refFile = 'IE13_T15_min';
        case 296695 % 3
            fn = {'Baseline','IE_11_T15_min','IE11_T25_min','IE13_T15_min','IE13_T25_min'};
%             fn = {'IE_11_T15_min','IE11_T25_min','IE13_T15_min','IE13_T25_min'};
            refFile = 'IE_11_T15_min';
        case 296954 % 4
            fn = {'baseline','IE13_T15_min','IE13_T25_min','IE11_T15_min','IE11_T25_min'};
%             fn = {'IE13_T15_min','IE13_T25_min','IE11_T15_min','IE11_T25_min'};
            refFile = 'IE13_T15_min';
    end % end switch
    EIT2 = struct;
    EIT2.info = EIT.info;
    for i = 1:length(fn)
       EIT2.(fn{i}) =  EIT.(fn{i});
    end
    EIT = EIT2;
    
    % Choose reference EIT FILE
    [fmdl, imdl] = mk_horse_model(EIT.info.elecShift);    % models fmdl and imdl
    refEitFile = EIT.(refFile);
    msel = find(fmdl.meas_select);
    vhCommonV = mean(refEitFile.fdata(msel, :), 2);
    
    % Noise compensation
    specs = struct; 
    specs.n = 90; 
    specs.type = 'meas'; 
    specs.thresh = 0.25;
    [imdl_comp, ~, ~] = eqadr(EIT, imdl, specs);
    
    % find clim
    CLIM = -inf;
    for i = 1: length(fn)
        eitFile = EIT.(fn{i});
        if ~isfield(eitFile, 'data')
            continue
        end
        eitFile.data = eitFile.data(msel, :);
        eitFile.fdata = eitFile.fdata(msel, :);
%         eitFile.triads = horse_breath_finder(eitFile.fdata, false);
        
        if contains(eitFile.name, '11')
            eitFile.triads = horse_breath_finder(eitFile.fdata, true);
        else
            eitFile.triads = horse_breath_finder(eitFile.fdata, false);
        end
        
        eitFile.imgr = inv_solve(imdl_comp, vhCommonV, eitFile.fdata);
        clim = max(eitFile.imgr.elem_data,[], 'all');
        CLIM = max([clim, CLIM]);
        
        EIT.(fn{i}) = eitFile;
    end
    
    % run analysis
    for i = 1:length(fn)
        eitFile = EIT.(fn{i});
        if ~isfield(eitFile, 'data')
            continue
        end
        eitFile.imgr.calc_colours.ref_level = 0;
        eitFile.imgr.calc_colours.clim = CLIM;
        EIT.(fn{i}) = eitFile;
        imgr = eitFile.imgr;
        imgs = calc_slices(imgr);
        
        if doPlotting
            plot_ventilation(eitFile, imgs);
            if contains(fn{i}, 'Perfusion') || contains(fn{i}, 'perfusion')
                plot_perfusion(eitFile, imgr, roi);
            end
        end % end if
        output_params_to_csv(eitFile, imgs);
%         output_params_to_csv(EIT);
    end % end for i

end % end function

function output_params_to_csv(eitFile, imgs)
    % Parameter Calculation
    
    if nargin == 1
        EIT = eitFile;
        try
            id = EIT.baseline.horseID;
        catch
            id = EIT.Baseline.horseID;
        end
        % ======================================================
%         saveDir = 'C:\Users\Mark\Dropbox\Joaquin\features\EELI_slopes\first_IE_ref\';
        
        % ======================================================
%         saveDir = 'C:\Users\Mark\Dropbox\Joaquin\features\mean_EELI_first_IE_as_reference\';
%         save_name = sprintf('%s%0.f_EELI_slopes.csv', saveDir, id);
%         header = {'filename', 'EELI_slope_AU_min'};
        
        % ======================================================
%         saveDir = 'C:\Users\Mark\Dropbox\Joaquin\features\mean_EELI_first_IE_as_reference\etc\';
%         save_name = sprintf('%s%0.f_etc.xlsx', saveDir, id);
%         header = {'filename', 'ect'};
        
        % ======================================================
        saveDir = 'C:\Users\Mark\Dropbox\Joaquin\features\__first_IE_as_reference\rvd\';
        save_name = sprintf('%s%0.f_rvd.xlsx', saveDir, id);
        header = {'filename', 'ect'};
        
        % ======================================================
        
        fileNames = fieldnames(EIT);
        fileNames = fileNames(2: end); % exclude baseline!

        Excel = actxserver('Excel.Application'); 
        set(Excel, 'Visible', 0);
        set(Excel,'DisplayAlerts', 0);
        Workbooks = Excel.Workbooks;
        
        for i = 1: length(fileNames)

            eitFile = EIT.(fileNames{i});
            imgs = calc_slices(eitFile.imgr);
            
            % ======================================================
%             res = [];
            % find means of reference IE ratio recording
%             thisEELI = EELI1_EELI2_mean_Z(imgs, eitFile.triads);
%             if i == 1
%                 paramFn = fieldnames(thisEELI);
%                 refMeans = zeros(length(paramFn), 1);
%             end
%             for j = 1: length(paramFn)
%                 if i == 1
%                     refMeans(j) = mean(thisEELI.(paramFn{j}));
%                 end
%                 res = [res, thisEELI.(paramFn{j}) - refMeans(j)];
%             end % end for j
            
            % ======================================================
%             ETC = expiratory_time_constant(imgs, eitFile.triads, eitFile.fs);
%             res = ETC.image;
%             paramFn = [];
            
            % ======================================================
            RVD = ventilation_delay(imgs, eitFile.triads, false);
            res = RVD.image;
            paramFn = [];
            
            % ======================================================
            
            if isempty(paramFn)
                col_ = res;
                table_out = array2table(col_);
            else
                table_out = array2table(res, 'VariableNames', paramFn);
            end
            writetable(table_out, save_name, 'Sheet', fileNames{i});
        end % end for i
        
        % delete the first 3 empty sheets
        [~, sheet_names] = xlsfinfo(save_name);
        Workbook = Workbooks.Open(save_name);
        Sheets = Excel.ActiveWorkBook.Sheets;
        index_adjust = 0;
        for i = 1:max(size(sheet_names))
            current_sheet = get(Sheets, 'Item', (i - index_adjust));
            if contains(sheet_names{i}, 'Sheet')
                invoke(current_sheet, 'Delete');
                index_adjust = index_adjust + 1;
            end 
        end
        Workbook.Save();
        Workbook.Close();
        Workbooks.Close;
        Excel.Quit();
        delete(Excel);        
        return
    end    
%     perfSave = sprintf('%s%s_features.csv', saveDir, eitFile.name);

%     saveDir = 'C:\Users\Mark\Dropbox\Joaquin\centre_of_ventilation\';

    saveDir = 'C:\Users\Mark\Dropbox\Joaquin\features\expiratory_time_constant\';
    save_name = sprintf('%s%s_etc.csv', saveDir, eitFile.name);

%     saveDir = 'C:\Users\Mark\Dropbox\Joaquin\features\mean_end_ex_deltas\';
%     save_name = sprintf('%s%s_meed.csv', saveDir, eitFile.name);

%     saveDir = 'C:\Users\Mark\Dropbox\Joaquin\features\mean_EELI_first_IE_as_reference\';
%     save_name = sprintf('%s%s_mean_EELI.csv', saveDir, eitFile.name);
    
%     saveDir = 'C:\Users\Mark\Dropbox\Joaquin\features\RVDi\';
%     save_name = sprintf('%s%s_RVDi.csv', saveDir, eitFile.name);
    
%     saveDir = 'C:\Users\Mark\Dropbox\Joaquin\features\_EILI_before_inspiratory_pause\RVD\';
%     save_name = sprintf('%s%s_sd_RVD.csv', saveDir, eitFile.name);
    output = [];
    header = {};
    count = 1;
    params = struct;
    if ~isempty(eitFile.triads)
%         params.RVD = ventilation_delay(imgs, eitFile.triads, false);
%         params.sdRVD = ventilation_delay(imgs, eitFile.triads, true);
        
%         return
%         sgtitle(eitFile.name);
%         printPDF( sprintf('%s%s', 'C:\Users\Mark\Dropbox\Joaquin\figures\ventilation_delay\', eitFile.name) );
        
        params.ETC = expiratory_time_constant(imgs, eitFile.triads, eitFile.fs);
%         sgtitle(eitFile.name);
%         printPDF( sprintf('%s%s', 'C:\Users\Mark\Dropbox\Joaquin\expiratory_time_constant\figures\', eitFile.name) );
        
%         params.cov = centre_of_ventilation(imgs, eitFile.triads);        
%         params.LungZ = lung_impedances(imgs, eitFile.triads);
        
%         params.ventDist = ventilation_distribution(imgs, eitFile.triads);
%         sgtitle(eitFile.name);
%         printPDF( sprintf('%s%s', 'C:\Users\Mark\Dropbox\Joaquin\figures\ventilation_distribution\', eitFile.name) );
        
%         params.meed = mean_end_ex_deltas(imgs, eitFile.triads);
        
%         params.EELIMean = end_ex2_mean_Z(imgs, eitFile.triads);
%         params.meanEELI = EELI1_EELI2_mean_Z(imgs, eitFile.triads);
%         [params.openingPhenom50, params.closingPhenom50] = opening_closing_phenom(imgs, eitFile.triads, eitFile.fs, 50);
%         [params.openingPhenom90, params.closingPhenom90] = opening_closing_phenom(imgs, eitFile.triads, eitFile.fs, 90);
%         [params.onOpening50, params.onClosing50] = pixels_not_exceeding_pct_lung_mean(imgs, eitFile.triads, 50);
%         [params.onOpening90, params.onClosing90] = pixels_not_exceeding_pct_lung_mean(imgs, eitFile.triads, 90);
%         params.timeToBaseline = time_to_baseline_after_exhale(eitFile.triads, eitFile.fs);
        
        % output to .csv
        paramNames = fieldnames(params);
        for i = 1:length(paramNames)
            thisParamName = paramNames{i};
            thisParam = params.(thisParamName);
            if isstruct(thisParam)
                subParams = fieldnames(thisParam);
                for j = 1:length(subParams)
                    header{count} = horzcat(thisParamName, '_', subParams{j});
                    output = [output, thisParam.(subParams{j})];
                    count = count + 1; 
                end % end for j
            else
                header{count} = thisParamName;
                output = [output, thisParam];
                count = count + 1;
            end % end if
        end % end for i
    end % end if
    
    table_out = array2table(output, 'VariableNames', header);
    writetable(table_out, save_name);
    
%     if eitFile.perfStart ~= 0
%         perfDist = perfusion_distribution(eitFile);
%         fn = fieldnames(perfDist);
%         perfOut = [];
%         perfHdr = {};
%         pCount = 1;
%         for i=1:length(fn)
%             perfHdr{pCount} = horzcat('perfDist','_',fn{i});
%             perfOut = [perfOut, perfDist.(fn{i})];
%             pCount = pCount + 1;
%         end % end for
%         pTable_out = array2table(perfOut, 'VariableNames', perfHdr);
%         writetable(pTable_out, perfSave);
%     end % end if

end % end function


function triads = insert_custom_triads(name)
    EELI1 = nan;
    EILI = nan;
    EELI2 = nan;
    switch name
        case '296439 IE11 T15 min'
            EELI1 = [178	614	1050	1482	1907	2338	2769	3179	3617	4046	4473	4871	5320	5760	6198	6618	7049	7483	7905	8323	8764];
            EILI = [327	766	1180	1627	2036	2467	2902	3350	3755	4191	4609	5045	5497	5896	6319	6746	7211	7632	8062	8467	8901];
            EELI2 = [614	1050	1482	1907	2338	2769	3179	3617	4046	4473	4871	5320	5760	6198	6618	7049	7483	7905	8323	8764	9188];
    end
    triads = [EELI1', EILI', EELI2'];
end

function plot_perfusion(eitFile, imgr, roi)
    global fig_path;
    global axsFS;
    global titleFS;
    global textFS;
    bigFig(2);
    refXVal     = 47;
    bothLungs   = (roi.LeftLung + roi.RightLung) > 0;
    lungMask    = find(bothLungs);
    start       = eitFile.perfStart - refXVal;
    stop        = eitFile.perfEnd;
    nFrames     = stop-start+1;
    
    perfImgr    = imgr;
    perfImgr.elem_data = imgr.elem_data(:,start:stop);
    clim        = max(perfImgr.elem_data,[],'all');
    perfImgr.calc_colours.ref_level  = 0;
    perfImgr.calc_colours.clim       = clim;
    perfImgs    = calc_slices(perfImgr);
    
    rsImgs      = reshape(perfImgs, [32^2, nFrames]);
    lungZ       = -mean(rsImgs(lungMask, :), 1);
    xax         = (1:nFrames) / eitFile.fs;
    subplot(2,5,1:5)
        plot(xax, lungZ);xlabel('Time (s)','FontSize',axsFS);ylabel('Lung Impedance (RU)','FontSize',axsFS);
    hold on; axis tight;
    timePoints = linspace(refXVal,nFrames, 7);
    timePoints = round(timePoints(2:6));
    
    % text Y position
    textY = min(lungZ);
    % plot vertical lines
    xVal = refXVal/eitFile.fs;
    xline(xVal);
    text(xVal+.1, textY, 'R', 'FontSize',textFS);
    for i=1:5
       xVal = timePoints(i)/eitFile.fs;
       xline(xVal);
       text(xVal+.1, textY, num2str(i), 'FontSize',textFS);
    end
    % show slices at time points
    for i=1:5
        subplot(2,5,i+5)
        show_slices(perfImgs(:,:,timePoints(i)) - perfImgs(:,:,1));
        title(num2str(i), 'FontSize',titleFS);
    end
    sgtitle(eitFile.name, 'FontSize',titleFS);
    % save
%     printPDF( sprintf('%s%s Perfusion Images', fig_path, eitFile.name) );
end


function plot_ventilation(eitFile, imgs)
    global fig_path;
%     global axsFS;
    global titleFS;
    fig     = bigFig(1);
    triads  = eitFile.triads;
    nBreaths= size(triads,1);
    if nBreaths > 0
        nRows = 4;
        subplot( nRows, nBreaths, 1:nBreaths*nRows-nBreaths);
        show_breath_boundaries(eitFile.fdata, triads);
        calc_colours('cmap_type', 'blue_red');
        colormap(calc_colours('colourmap'));
        hold on;
        for j=1:nBreaths
            subplot( nRows, nBreaths, nBreaths*nRows-nBreaths+j);
            image(calc_breath_delta_z(imgs, triads(j,:)).*1000 + 128);
            title(num2str(j));
            axis equal; axis off;
        end
        hold off;
    else
        show_breath_boundaries(eitFile.fdata, triads);
    end
    sgtitle(remove_underscores(eitFile.name), 'FontSize', titleFS);
    printPDF( sprintf('%s%s Global Z', fig_path, eitFile.name) );
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

