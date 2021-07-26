files = {'C:\Users\Mark\Documents\EIT\Will_Ditcham\3deit\PB_Breath_Hold_3D_Trial.eit', 'C:\Users\Mark\Documents\EIT\Will_Ditcham\3deit\PB_Upright_3D_Trial.eit'};

% make models
thorax = shape_library('get','adult_male','boundary');
shape = { 2, {thorax}, [4, 40], 0.5};

% Set up the electrodes with 1:16 on bottom and 17:32 on top
elec_pos = [16, 0, 0.75, 1.25];
elec_spec = [0.01];
% elec_spec = [0.5, 0, 0.1]; % Radius of circular electrodes
fmdl = ng_mk_extruded_model(shape, elec_pos, elec_spec);

% Set up square pattern
idx = reshape([1: 32], 2, [])'; %matrix of 32 elem numbered from 1:32, shaped with row length= L then transposed so numbering goes across rows
idx(2: 2: end, :) = fliplr(idx(2: 2: end, :)); % flips every second row to get proper square electrode arrangement
idx = fliplr(idx);
idx = [idx(end - 3: end, :); idx(1: end - 4, :)]; 
fmdl.electrode(idx) = fmdl.electrode(:);

% Stim pattern
[fmdl.stimulation, fmdl.meas_select] = mk_stim_patterns(32, 1, [0, 5], [0, 5], {'no_meas_current_next2'}, 1); % Skip 4
fmdl.normalize = 0;

% img = mk_image(fmdl, 1);
% img.elem_data([fmdl.mat_idx{2}; fmdl.mat_idx{3}]) = 0.3;
% img.fwd_model.normalize_measurements = 0;

% imdl
vopt.imgsz = [32 32];
vopt.square_pixels = true;
vopt.zvec = linspace(0.5, 1.5, 4);
vopt.save_memory = 1;
opt.noise_figure = 0.5;
opt.keep_intermediate_results = true;
[imdl_t, opt.distr] = GREIT3D_distribution(fmdl, vopt);
imdl = mk_GREIT_model(imdl_t, 0.2, [], opt);

% pre-process
[D, imdl_comp] = wr_pp(files, imdl, 300);

% make images
dataFn = fieldnames(D);
for j = 1:2
    use_data = D.(dataFn{j}).useData;
    imgr = inv_solve(imdl_comp, mean(use_data,2), use_data);
    imgs = calc_slices(imgr, [inf inf 1]);
    triads = horse_breath_finder(use_data);
    vImg = img_slices(:, :, triads(1, 2)) - img_slices(:, :, triads(1, 1));
    show_slices(vImg)
    eitFile.triads = triads;
    eitFile.fdata = use_data;
    nameslpit = strsplit(files{j}, '\');
    eitFile.name = nameslpit{end};
    plot_ventilation(eitFile, imgs);
end


function plot_ventilation(eitFile, imgs)
    global fig_path;
%     global axsFS;
    global titleFS;
    fig     = figure();
    triads  = eitFile.triads;
    nBreaths = size(triads,1);
    if nBreaths > 0
        subplot( 2, 1, 1);
        gsig = sum(eitFile.fdata, 1);
        plot(gsig, 'LineWidth', 2); hold on;
        title('Global signal');
        for i = 1:size(triads,1)
            xline(triads(i, 1), 'r');
            xline(triads(i, 2), 'b');
        end
        calc_colours('cmap_type', 'blue_red');
        colormap(calc_colours('colourmap'));
        hold on;
        showimgs = zeros(32, 32, nBreaths);
        for j=1:nBreaths
            showimgs(:, :, j) = calc_breath_delta_z(imgs, triads(j,:));
        end
        subplot(2, 1, 2);
        imagesc(mk_mosaic(showimgs, 3, [], round(nBreaths/4)));
        title('Tidal images');
        axis equal; axis off;
        hold off;
    else
        show_breath_boundaries(eitFile.fdata, triads);
    end
    sgtitle(remove_underscores(eitFile.name), 'FontSize', 16);
end