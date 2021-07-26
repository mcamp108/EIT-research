global excel_dir;
global customRoi;
path_name   = 'C:\Users\Mark\Dropbox\Joaquin\data\';
titleFS     = 16;
axsFS       = 16;
textFS      = 12;
fig_path    = 'C:\Users\Mark\Documents\EIT\proj\vq_apnea_breathing\';
excel_dir   = 'C:\Users\Mark\Documents\EIT\proj\vq_apnea_breathing\';
doPlotting  = false;
mkPaperFig  = true;
do_apnea    = true;
roi = horse_roi('vq_full');
customRoi = zeros(32, 32);

% model specs
imgsize = [32 32];
radius = 0.2; % - requested weighting matrix  (recommend 0.2 for 16 electrodes)
weight = []; % - weighting matrix (weighting of noise vs signal). Can be empty options.noise_figure is specified
opt.noise_figure = 0.5;
opt.imgsz = imgsize;
opt.keep_intermediate_results = true;
opt.save_memory = 1;
opt.square_pixels = true;

%%
% dat = [1.15, 1.31, 1.03, 1.18, 1.16, 1.01];
% err = [0.07, 0.13, 0.04, 0.13, 0.09, 0.07];
% 
% % plot (dat, 'k');
% hold on;
% for i = 1:length(err)
%     errorbar(i, dat(i), err(i), 's', 'LineWidth', 2, 'Color', 'k', 'CapSize', 12, 'MarkerSize', 10);
% end
% xlim([0, 7]);
% ylim([0.4, 1.6]);
% ylabel('$\dot{V}/Q$', 'Interpreter', 'latex');
% xlabel('Lung region');
% lbls = {''; 'Right ventral'; 'Right intermediate'; 'Right dorsal'; 'Left ventral'; 'Left intermediate'; 'Left dorsal'; ''};
% xticklabels(lbls);
% xtickangle(-15);

%%
% horse IDs
IDS = [296439, 296695, 296954, 296955];

for i = 1:length(IDS)
    for do_apnea = [true false]
        close all;
        horseID = IDS(i);
        disp(horseID);

        % load data
        EIT = load_joaquin_data(path_name, horseID);

        % make models
        [fmdl, ~] = mk_horse_model(EIT.info.elecShift);
        img = mk_image(fmdl, 1); % background conductivity
        img.elem_data([fmdl.mat_idx{2};fmdl.mat_idx{3}]) = 0.3;
        img.fwd_model.normalize_measurements = 0;
        imdl = mk_GREIT_model(img, radius, weight, opt);
        msel = find(fmdl.meas_select); % used for plotting

        % Noise compensation
        [dd, aux] = eidors_readdata('C:\Users\Mark\Dropbox\Joaquin\data\296439\296439_IE13_Perfusion_breathing.eit');
        specs = struct; specs.n = 90; specs.type = 'meas';
        [imdl_comp, removed, scores] = eqadr(dd, imdl, opt);

        % do images
        switch horseID
            case 296439
                if do_apnea % perfusion apnea
                    t25 = EIT.IE13_Perfusion_apnea; blm_25 = [298, 381]; plm_25 = [1016, 1822, 2124];
                    t50 = EIT.IE11_Perfusion_apnea; blm_50 = [175, 356]; plm_50 = [993, 1804, 2334];
                else % perfusion breathing
                    t25 = EIT.IE13_Perfusion_breathing; blm_25 = [2885, 2973]; plm_25 = [1975, 2885, 3760];
                    t50 = EIT.IE11_Perfusion_breathing; blm_50 = [1196, 1364]; plm_50 = [767, 1196, 2055];
                end % end if
            case 296695
                if do_apnea
                    t25 = EIT.IE11_Perfusion_apnea; blm_25 = [97, 247]; plm_25 = [797, 2186, 3165];
                    t50 = EIT.IE13_Perfusion_Apnea; blm_50 = [142, 1]; plm_50 = [142, 1720, 3335]; % half breath
                else
                    t25 = EIT.IE11_Perfusion_breathing; blm_25 = [2018, 2195]; plm_25 = [824, 2017, 3208];
                    t50 = EIT.IE13_Perfusion_breathing; blm_50 = [1585, 1674]; plm_50 = [785, 1585, 3195];
                end % end if
            case 296954
                if do_apnea
                    t25 = EIT.IE_13_Perfusion_apnea; blm_25 = [219, 330]; plm_25 = [916, 2215, 3118];
                    t50 = EIT.IE11_Perfusion_apnea; blm_50 = [14, 199]; plm_50 = [416, 1583, 2639];
                else
                    t25 = EIT.IE13_Perfusion_breathing; blm_25 = [1598, 1687]; plm_25 = [795, 1598, 2799];
                    t50 = EIT.IE11_Perfusion_breathing; blm_50 = [1321, 1473]; plm_50 = [525, 1321, 2511];
                end % end if
            case 296955
                if do_apnea
                    t25 = EIT.IE_11_perfusion_apnea; blm_25 = [167, 354]; plm_25 = [605, 1737, 1959];
                    t50 = EIT.IE_13_Perfusion_apnea; blm_50 = [109, 192]; plm_50 = [890, 1811, 2227];
                else
                    t25 = EIT.IE_11_Perfusion_breathing; blm_25 = [2384, 2582]; plm_25 = [808, 2384, 3972];
                    t50 = EIT.IE13_Perfusion_breathing; blm_50 = [2069, 2155]; plm_50 = [839, 2069, 3270];
                end % end if

        end % end switch
        if do_apnea
            t25.name = sprintf('%s T25 apnea', num2str(horseID));
            t50.name = sprintf('%s T50 apnea', num2str(horseID));
        else
            t25.name = sprintf('%s T25 breathing', num2str(horseID));
            t50.name = sprintf('%s T50 breathing', num2str(horseID));
        end % end if

%         show_vq_landmarks( sum(t25.fdata(msel, :), 1), blm_25, plm_25);
%         sgtitle(t25.name);
%         saveas(gcf, horzcat(fig_path, t25.name, ' landmarks.png'));
        [vqImg11, v_raw11, q_raw11, vq_raw11] = vq_img(t25, imdl_comp, blm_25, plm_25);

%         show_vq_landmarks( sum(t50.fdata(msel, :), 1), blm_50, plm_50);
%         sgtitle(t50.name);
%         saveas(gcf, horzcat(fig_path, t50.name, ' landmarks.png'));
        [vqImg13, v_raw13, q_raw13, vq_raw13] = vq_img(t50, imdl_comp, blm_50, plm_50);

        % make vq comparison figure
        figure(2); 
        figImg = [vq_raw11, vq_raw13];
        figImg(figImg == -1) = 1; % fill inside
        padding = nan(1, size(figImg, 2));
        figImg = [padding; figImg; padding];

        showImg = figImg;
        showImg(showImg > 1.6) = 1.6;
        showImg(showImg < 0.4) = 0.4;

        showImg = draw_boundary(showImg);
        clim = max(abs(showImg), [], 'all', 'omitnan') - 1;
        imagesc(showImg);
        cmap = colormap;
        cmap = flipud([cmap(2: end, :); [0.5, 0.5, 0.5]]);
        colormap(cmap);
        caxis([0.35, 1.65]);
        axis equal; axis image; axis tight; axis off;
        colorbar();
        saveas(gcf, horzcat(fig_path, num2str(horseID), '_apnea=', num2str(do_apnea), '.png'));

        % print stats (lazy)
    %     [means11, sd11, medians11] = calc_roi_stats(vq_raw11);
    %     [means13, sd13, medians13] = calc_roi_stats(vq_raw13);
    %     disp([means11, means13, sd11, sd13, medians11, medians13]);
    % 
        % save vq image to excel 
        vq_raw11(vq_raw11 == -1) = nan; % (set inside, non-lung to nan)
        vq_raw13(vq_raw13 == -1) = nan; % (set inside, non-lung to nan)
        save_img_to_xl(vq_raw11, t25.name);
        save_img_to_xl(vq_raw13, t50.name);
    end % end for do_apnea
end % end for i


function save_img_to_xl(img, sheet)
    global excel_dir
    writematrix(img, sprintf('%sVQ_images_ap_br.xlsx', excel_dir), 'Sheet', sheet);
end


function [vq_show, v_raw, q_raw, vq_raw] = vq_img(eitFile, imdl, breath, perf)
    roi = horse_roi('vq_full');
    roi = roi.BothLungs;
    b1 = breath(1); b2 = breath(2);
    p1 = perf(1); p2 = perf(2); p3 = perf(3);
    
    % find relative indices
    s1 = min([b1, b2, p1, p2, p3]);
    sN = max([b1, b2, p1, p2, p3]);
    b1 = b1 - s1 + 1;
    b2 = b2 - s1 + 1;
    p1 = p1 - s1 + 1;
    p2 = p2 - s1 + 1;
    
    % make images
    eitFile.data = eitFile.data(:, s1: sN);
    eitFile.fdata = eitFile.fdata(:, s1: sN);
    eitFile.triads = horse_breath_finder(eitFile.fdata);
    vh = eitFile.fdata(:, p1);
    eitFile.imgr = inv_solve(imdl, vh, eitFile.fdata);
    imgr = eitFile.imgr;
    
    % ventilation image
    vImg = imgr;
    vImg.elem_data = imgr.elem_data(:, b2) - imgr.elem_data(:, b1);
    vSlc = calc_slices(vImg) .* roi;
    vSlc(vSlc > 0) = 0;
    
    % perfusion image
    qImg = imgr;
    qImg.elem_data = imgr.elem_data(:, p2);
    qSlc = calc_slices(qImg) .* roi;
    qSlc(qSlc < 0) = 0;
    
%     refImg = qSlc;
%     tmp = show_slices(refImg);
%     tmp(refImg==0) = 128;
%     image(tmp); axis equal; axis image; axis tight; axis off;
    
    inside = ~isnan(qSlc);
    bkg = inside == 0;
    inside = logical(inside - roi);
    
    clim = max([-vSlc(:); qSlc(:)], [], 'omitnan');
    vImg.calc_colours.clim = clim;
    qImg.calc_colours.clim = clim;
    
    v_raw = -vSlc;
    q_raw = qSlc;
    
    v_ = calc_colours(v_raw);
    q_ = calc_colours(q_raw);
    vq_raw = (v_ ./ q_);

    % -- find 0 and inf pixels --
    temp = v_raw./ q_raw;
    temp(inside) = 1;
    vq_raw(isinf(temp)) = 2;
    vq_raw(isnan(temp)) = 1; % call 0/0 even VQ.
    vq_raw(temp == 0) = 0.05;
    
    % -- V, Q, and VQ iamge --
    v_show = calc_colours(-v_raw);
    v_show(inside) = 128;
    v_show(bkg) = nan;

    q_(inside) = 128;
    q_(bkg) = nan;
    
    vq_show = vq_raw - 1;
    vq_show(bkg) = nan;
    vq_show(inside) = 0;
    vq_show = calc_colours(-vq_show);
    vq_show(bkg) = nan;
    
    combined = [v_show, q_, vq_show];
    padding = nan(1, size(combined, 2));
    figImg = [padding; combined; padding];
    figImg = draw_boundary(figImg, 128);
    
    figure();
    imagesc(figImg);
    axis equal; axis image; axis tight;
    ax = gca;
    imgCols = size(vq_show, 2);
    imgC = imgCols / 2;
    ax.XTick = [imgC, imgC + imgCols, imgC + imgCols * 2];
    ax.TickLabelInterpreter = 'latex';
    ax.XTickLabel{1} = '$\dot{V}$';
    ax.XTickLabel{2} = '$Q$';
    ax.XTickLabel{3} = '$\dot{V}/Q$';
    title(eitFile.name);
    calc_colours('cmap_type', 'blue_red');
    colormap(calc_colours('colourmap'));
    colorbar();
    % set inside of raw images to nan
    vq_raw(inside) = -1;
    vq_raw(bkg) = nan;
    
end

function show_vq_landmarks(tbv, bl, pl)
    figure();
    plot(tbv, 'LineWidth', 2);
    xline(bl(1), '--', 'Color', '#D95319','LineWidth', 2);
    xline(bl(2), '--', 'Color', '#D95319','LineWidth', 2);
    xline(pl(1), '--', 'Color', '#D95319','LineWidth', 2);
    xline(pl(2), '--', 'Color', '#EDB120','LineWidth', 2);
    xline(pl(3), '--', 'Color', '#7E2F8E','LineWidth', 2);
end

function lung = isolate_lung(img, roi)
    lung = img .* roi;
    lung(lung == 0) = nan;
end

function [reg1, reg2, reg3] = divide_in_three_rows(lung, roi)
    rowIdx = find(sum(roi, 2));
    r1 = rowIdx(1:6);
    r2 = rowIdx(7:12);
    r3 = rowIdx(13:18);
    reg1 = lung(r1, :);
    reg2 = lung(r2, :);
    reg3 = lung(r3, :);
end

function [r1, r2, r3] = divide_in_three_rows_idx(roi)
    rowIdx = find(sum(roi, 2));
    r1 = rowIdx(1:6);
    r2 = rowIdx(7:12);
    r3 = rowIdx(13:18);
end


function [means, sds, medians] = calc_roi_stats(img)
    roi = horse_roi();
    RL = isolate_lung(img, roi.RightLung);
    means = zeros(6, 1);
    sds = zeros(6, 1);
    medians = zeros(6, 1);
    [r1, r2, r3] = divide_in_three_rows(RL, roi.BothLungs);
    means(1) = mean(r1(:), 'omitnan');
    means(2) = mean(r2(:), 'omitnan');
    means(3) = mean(r3(:), 'omitnan');
    sds(1) = std(r1(:), 'omitnan');
    sds(2) = std(r2(:), 'omitnan');
    sds(3) = std(r3(:), 'omitnan');
    medians(1) = median(r1(:), 'omitnan');
    medians(2) = median(r2(:), 'omitnan');
    medians(3) = median(r3(:), 'omitnan');
    
    LL = isolate_lung(img, roi.LeftLung);
    [l1, l2, l3] = divide_in_three_rows(LL, roi.BothLungs);
    means(4) = mean(l1(:), 'omitnan');
    means(5) = mean(l2(:), 'omitnan');
    means(6) = mean(l3(:), 'omitnan');
    sds(4) = std(l1(:), 'omitnan');
    sds(5) = std(l2(:), 'omitnan');
    sds(6) = std(l3(:), 'omitnan');
    medians(4) = median(l1(:), 'omitnan');
    medians(5) = median(l2(:), 'omitnan');
    medians(6) = median(l3(:), 'omitnan');
    

%     [r1, r2, r3] = divide_in_three_rows_idx(RL);
%     r1 = r1 + 4; r2=r2+4;r3=r3+4;
%     r1 = [4;r1];
%     r3 = [r3;23];
%     tempR = roi.RightLung * 1;
%     tempL = roi.LeftLung * 1;
%     tempR(r1, :) = tempR(r1, :) * 1;
%     tempR(r2, :) = tempR(r2, :) * 2;
%     tempR(r3, :) = tempR(r3, :) * 3;
%     tempL(r1, :) = tempL(r1, :) * 4;
%     tempL(r2, :) = tempL(r2, :) * 5;
%     tempL(r3, :) = tempL(r3, :) * 6;
%     temp = tempR+tempL;
%     temp(isnan(img)) = nan;
%     temp(temp==0) = 7;
%     padding = nan(1, size(temp, 2));
%     temp = [padding; temp; padding];
%     temp = draw_boundary(temp, 7);
%     imagesc(temp);
%     colormap('bone');
%     axis equal; axis image; axis off;
    
end
   

function roi = find_roi(imgr)
    roi = horse_roi();
    roi = roi.BothLungs;
    roi(11:17, 16:19) = 0;
    return
    tmp = calc_slices(imgr);
    sdImg = std(tmp, [], 3);
    sdImg(isnan(sdImg)) = min(sdImg(:));
    sdImg = ( sdImg - min(sdImg(:)) ) ./ ( max(sdImg(:)) - min(sdImg(:)) );
    tmp = sdImg;
%     roi = tmp >= 0.4;
%     roi(:, 16:19) = 0;

    roi = tmp >= 0.4;
    roi(11:17, 16:19) = 0;
end

function out = mk_fig_slc(slc, roi, inside, bkg)

    slc_ = show_slices(slc);
    slc_ = slc_ .* roi;
%     mu = mean(slc_(roi));
%     out = slc_ ./ mu;
%     out(inside) = 0;
    out(bkg) = nan;
end


function newImg = draw_boundary(img, bb)
    if nargin == 1
        bb = 1;
    end
    height = size(img, 1);
    width = size(img, 2);
    bkg = find(isnan(img));
    newImg = img;
    for i = 1: length(bkg)
        elemIdx = bkg(i);
        elem = img(elemIdx);
        neigh = find_neighbours(elemIdx, width, height);
        if all(  isnan(img(neigh))  )
            newImg(elemIdx) = bb;
        else
            continue
        end
    end
end


function neigh = find_neighbours(idx, width, height)
    neigh = [];
    if mod(idx - 1, height) ~= 0 % U
        neigh = [neigh, idx - 1];
    end
    if mod(idx, height) ~= 0 % D
        neigh = [neigh, idx + 1];
    end
    if idx > height % L
        neigh = [neigh, idx - height];
    end
    if idx <= height * (width - 1) % R
        neigh = [neigh, idx + height];
    end
end