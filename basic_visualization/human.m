% Reconstruct human EIT data
[fmdl, imdl] = mk_model('human');
% files = {'C:\Users\Mark\Documents\EIT\Will_Ditcham\Peter_Test_2_-_Left_recumbency.eit',...
%         'C:\Users\Mark\Documents\EIT\Will_Ditcham\Peter_Test_3_-_Right_Recumbency.eit'};

% files = {'C:\Users\Mark\Dropbox\EIT 2021\Study files\16-02-21_Richter_DIAMOND_KG2337342_nachLonge.eit',...
%         'C:\Users\Mark\Dropbox\EIT 2021\Study files\16-02-21_Richter_DIAMOND_KG2337342_vorLonge.eit',...
%         'C:\Users\Mark\Dropbox\EIT 2021\Study files\18-02-21_Pianta_ORIAS_KG2361200_vorLonge.eit',...
%         'C:\Users\Mark\Dropbox\EIT 2021\Study files\18-02-21_Pianta_ORIAS_KG2361200_VT_nachLonge.eit',...
%         'C:\Users\Mark\Dropbox\EIT 2021\Study files\2021_23_02_Brun_SONYIVCH_KG2240697_nachLonge.eit',...
%         'C:\Users\Mark\Dropbox\EIT 2021\Study files\2021_23_02_Brun_SONYIVCH_KG2240697_vorLonge.eit'};

% files = {'C:\Users\Mark\Documents\EIT\Will_Ditcham\#2 Will Mitcham (22-04-2021)-20210506T001539Z-001\Will_-_Left_Laterel_Recumbency.eit',...
%     'C:\Users\Mark\Documents\EIT\Will_Ditcham\#2 Will Mitcham (22-04-2021)-20210506T001539Z-001\Will_-_Prone_(on_front).eit',...
%     'C:\Users\Mark\Documents\EIT\Will_Ditcham\#2 Will Mitcham (22-04-2021)-20210506T001539Z-001\Will_-_Right_Lateral_Recumbency.eit',...
%     'C:\Users\Mark\Documents\EIT\Will_Ditcham\#2 Will Mitcham (22-04-2021)-20210506T001539Z-001\Will_-_Seated_upright_Post-Protocol.eit',...
%     'C:\Users\Mark\Documents\EIT\Will_Ditcham\#2 Will Mitcham (22-04-2021)-20210506T001539Z-001\WIll_-_Seated_Upright_Pre-Protocol.eit',...
%     'C:\Users\Mark\Documents\EIT\Will_Ditcham\#2 Will Mitcham (22-04-2021)-20210506T001539Z-001\Will_-_Supine_(on_back).eit'};
files = {'C:\Users\Mark\Documents\EIT\Will_Ditcham\#2 Will Mitcham (22-04-2021)-20210506T001539Z-001\Will_-_Left_Laterel_Recumbency.eit'};

saveDir = 'C:\Users\Mark\Documents\EIT\Will_Ditcham\#2 Will Mitcham (22-04-2021)-20210506T001539Z-001\';
[D, imdl_comp] = wr_pp(files, imdl, 60);
fn = fieldnames(D);

for i = 1:length(files)
    imgr = inv_solve(imdl_comp, mean(D.(fn{i}).useData, 2), D.(fn{i}).useData);
    img_slices = calc_slices(imgr, [inf inf 1]);
    img_slices(isnan(img_slices)) = 0;
    img_slices = img_slices .* lung_segmentation();
    lungZ = squeeze( sum( img_slices, [1 2]) )';
    [end_in, end_ex, Ti_by_Tt, BF] = select_breaths(-lungZ, D.seq_1.fs);

    tvImg = imgr;
    tvImg.elem_data = imgr.elem_data(:, end_in) - imgr.elem_data(:, end_ex(1, :));
%     tvImg.elem_data(isnan(tvImg.elem_data)) = 0;

    tvImg.calc_colours.clim = max(tvImg.elem_data, [], 'all');
    imgs = calc_colours(calc_slices(tvImg, [inf inf 1]));
    
    fig = figure();
    fig.Units = 'normalized';
    fig.OuterPosition = [0 0 1 1];
    
    show_slices(imgs);
    ttl = strsplit(files{i}, '\');
    ttl = ttl{end};
    ttl = strrep(ttl, '.eit', '');
    title(remove_underscores(ttl));
    printPDF(horzcat(saveDir, ttl));
end
for i = 1:length(files)
    figure();
    plot_file(files{i});
end


function plot_file(file)
    [vv, aux]= eidors_readdata(file);
    startIdx = 1;
    endIdx = size(vv, 2);
    t_rel = aux.t_rel(:, startIdx: endIdx);
    fs = 1e6 ./ median( diff(t_rel) ); %framerate is median dif of time points/ 1000000 (convert to s)
    vv = vv(:, startIdx: endIdx);
    rvv = real(vv);
    xax = (1:length(vv)) ./ fs;

    impedanceFactor =  2.048 / (2^12 * 0.173 * 0.003) / 2^15; % = 0.9633 / 2^15;
    elec_impedance = aux.elec_impedance * impedanceFactor;
    elec_impedance = abs(elec_impedance);
    ei = mean(elec_impedance, 2);
    
    fig = figure();
    fig.Units = 'normalized';
    fig.OuterPosition = [0 0 1 1];
    
    subplot(2, 1, 1)
    plot(xax, sum(rvv,1));
    ttl = strsplit(file, '\');
    ttl = ttl{end};
    ttl = strrep(ttl, '.eit', '');
    title(remove_underscores(ttl));
    xlabel('Time (s)'); ylabel('Global Signal (V)');

    subplot(2,1,2)
    b = bar(ei,'FaceColor','flat'); hold on;
    b.CData(:,:) = repmat([0 0 1], 32, 1);
    eb = errorbar(1:32, ei, max(elec_impedance, [], 2)-ei, min(elec_impedance, [], 2)- ei);
    set(eb,'Color',[0, 0, 0],'LineStyle','none');
    xlim([0, 33]);
    xlabel('Electrode Number'); ylabel('Impedance (Ohms)');
    title 'Median Electrode Impedance';
    hold off
end


function [fmdl, imdl] = mk_model(animal)
    switch animal
        case 'human'
            thorax = shape_library('get','adult_male','boundary');
            shape = { 2, {thorax}, [4,40], 0.5};

            % 32 electrodes evenly spaced
            elec_pos = [32, 0, 1];
            elec_spec = [0.01];
            fmdl = ng_mk_extruded_model(shape, elec_pos, elec_spec);
            % shift belt
%             fmdl = shift_electrodes(fmdl, 16);
%             fmdl = shift_electrodes(fmdl, 'lr');
            % Stim pattern
            [fmdl.stimulation, fmdl.meas_select] = mk_stim_patterns(32, 1, [0, 5], [0, 5], {'no_meas_current_next2'}, 1); % Skip 4
            fmdl.normalize = 0;

            % imdl
            vopt.imgsz = [32 32];
            vopt.square_pixels = true;
            vopt.zvec = linspace(0.5, 1.5, 4);
            vopt.save_memory = 1;
            opt.noise_figure = 0.5;
            opt.keep_intermediate_results = true;
            [imdl_t, opt.distr] = GREIT3D_distribution(fmdl, vopt);
            imdl = mk_GREIT_model(imdl_t, 0.2, [], opt);
        case 'horse'
            [fmdl, imdl] = mk_horse_model(10);
    end     
end


function [end_in, end_ex, Ti_by_Tt, BF] = select_breaths(tbv, FR)

    % Breath boundaries
    breaths= find_breaths(tbv);
    end_in= breaths.ins_idx;
    end_ex= breaths.exp_idx;
    exp_to_exp= (diff(end_ex, 1))';
    Te= ((end_ex(2,:)- end_in))';
    Ti= (abs(end_ex(1,:)- end_in))';
    
    reject1 = abs( (end_in - end_ex(1,:)) ./ (end_in - end_ex(2,:)) );
    reject2 = abs( (end_in - end_ex(2,:)) ./ (end_in - end_ex(1,:)) );
    reject3 = (tbv(end_in) - tbv(min(end_ex,[],1))) ./ abs(diff(tbv(end_ex),1));
    
    keep = ((reject1 < 2.5) + (reject2 < 2) + (reject3 > 3)) == 3;
    
    end_in = end_in(keep);
    end_ex = end_ex(:, keep);
    Ti = Ti(keep);
    Te = Te(keep);
    exp_to_exp = exp_to_exp(keep);
    
    Ti_by_Tt= Ti ./ (Ti + Te);
    BF = 60./ (exp_to_exp./ FR); % instantaneous breaths per minute

end % end function

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
