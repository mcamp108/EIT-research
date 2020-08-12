% Chapter 5 Figures
%--------------------------------------------------------------------------
figure('units','normalized','outerposition',[0 0 1 1]); clf;
cd 'C:\Users\Mark\Documents\GraduateStudies\LAB\EIT-restraint\zzMC\data\Mali Weighted Restraint';
SAVEDIR = 'C:\Users\Mark\Documents\GraduateStudies\LAB\Thesis\files\figures\Chap_5\';
dirs= ls;
% get model
[fmdl, imdl]= mk_weighted_restraint_model(); % model
lung_roi= lung_segmentation(); % lung roi
mid_col= size(lung_roi, 2)/2;
lung_divide_idx= mid_col* size(lung_roi, 1) + 1;
is_lung= find(lung_roi);
left_lung_idx= is_lung(is_lung>= lung_divide_idx);
right_lung_idx= is_lung(is_lung< lung_divide_idx);
recordings= {'sref', 'pref', 'wref', 'wepos', 'epos'}; % change these with your own abbreviations for trials you have.
bsln = 2; % prone unweighted default is recording #2. This is used as the reference for all recordings.
i=5;
folder= dirs(i, :);
%--------------------------------------------------------------------------
% Weighted restraint pre-processing
clip= 50;
D = struct;
RMTHRESHOLD = 0.25;
msel= imdl.fwd_model.meas_select;
mm = find(msel);
% sequences = {};
file = 2;
f = files{file};
field = horzcat('seq', num2str(file));
[ dd, D.(field).aux ] = eidors_readdata(f);
D.(field).aux.elec_impedance = D.(field).aux.elec_impedance(:, 1: size(dd, 2));
t_rel = D.(field).aux.t_rel(:, 1: size(dd, 2));
FR = 1e6 ./ median( diff(t_rel) ); %framerate is median dif of time points/ 1000000 (convert to s)
dd = dd( :, 1: min(size(dd, 2), round(FR * timeSeriesLength)) ); % ensure all time series are timeSeriesLength or less
D.(field).data = dd;
useData = real(dd(mm, :));
useData = lowpass(useData', 1, FR)'; % 1 Hz lowpass
useData = useData(:, clip:(size(useData,2)-clip) ); % trim filter edge artifacts
D.(field).useData = useData;
%--------------------------------------------------------------------------
% find worst N electrodes
[rmElecs, scores] = worst_n_elecs(dd, imdl, 6);
rmElecs = rmElecs(scores >= RMTHRESHOLD);
imdl_comp_resolve = comp_RM_bad_elec(imdl, 26);
imgr_resolve = inv_solve(imdl_comp_resolve, mean(useData,2), useData);
imdl_comp_noResolve = comp_RM_bad_elec(imdl, [21 26 31]);
imgr_noResolve = inv_solve(imdl_comp_noResolve, mean(useData,2), useData);
%--------------------------------------------------------------------------  
% Data pre-processing and breath selection
img_slices_resolve(isnan(img_slices_resolve))= 0;
img_slices = img_slices_resolve .* lung_segmentation();
lungZ = squeeze( sum( img_slices, [1 2]) )';
[end_in, end_ex, Ti_by_Tt, BF]= select_breaths(-lungZ, FR);
%--------------------------------------------------------------------------        
% make figures
frm = 1;
imgr_compare = imgr_resolve;
imgr_compare.elem_data = [  imgr_noResolve.elem_data(:,end_ex(1,frm)),  imgr_resolve.elem_data(:,end_ex(1,frm)),zeros(1878,1),...
                            imgr_noResolve.elem_data(:,end_in(frm)),    imgr_resolve.elem_data(:,end_in(frm)),  zeros(1878,1),...
                            imgr_noResolve.elem_data(:,end_ex(2,frm)),  imgr_resolve.elem_data(:,end_ex(2,frm)),zeros(1878,1)
                            ];
slices_compare = calc_colours( calc_slices(imgr_compare, [inf inf 1]) );
bkg = slices_compare(:,:,1)==1;
slices_compare(:,:,3) = slices_compare(:,:,2) - slices_compare(:,:,1) + 125;
slices_compare(:,:,6) = slices_compare(:,:,5) - slices_compare(:,:,4) + 125;
slices_compare(:,:,9) = slices_compare(:,:,8) - slices_compare(:,:,7) + 125;
for i=[3 6 9]
    temp = slices_compare(:,:,i); temp(bkg) = 1;
    slices_compare(:,:,i) = temp;
end
final_img = mk_mosaic(slices_compare(6:27,:,:),[1, 1],[], 3); 
final_img(isnan(final_img)) = 0;
image(final_img);axis equal; axis image; axis tight;
colorbar;
fg = gcf();
ax = fg.Children(2);
incrX = linspace(1,98,7);
ax.XTick = incrX([2,4,6]);
ax.XTickLabel = {'No Score Resolution';'Score Resolution';'Difference'};
incrY = linspace(1,65,7);
ax.YTick = incrY([2,4,6]);
ax.YTickLabel = {'End Expiration';'End Inspiration';'End Expiration'};
ax.XAxis.FontSize=20; ax.XAxis.FontWeight='bold';
ax.YAxis.FontSize=20; ax.YAxis.FontWeight='bold'; ax.YAxis.TickLabelRotation=-45;
printPDF(sprintf('%sResolveImgs',SAVEDIR));

% =========================================================================

function [end_in, end_ex, Ti_by_Tt, BF]= select_breaths(tbv, FR)

    % Breath boundaries
    breaths= find_breaths(tbv);
    end_in= breaths.ins_idx;
    end_ex= breaths.exp_idx;
    exp_to_exp= (diff(end_ex, 1))';
    Te= ((end_ex(2,:)- end_in))';
    Ti= (abs(end_ex(1,:)- end_in))';
    
    reject1 = abs( (end_in - end_ex(1,:)) ./ (end_in - end_ex(2,:)) );
    reject2 = abs( (end_in - end_ex(2,:)) ./ (end_in - end_ex(1,:)) );
    reject3 = abs(diff(tbv(end_ex), 1));
    
    keep = ((reject1 < 2.5) + (reject2 < 2)) == 2;
    
    end_in= end_in(keep);
    end_ex= end_ex(:,keep);
    Ti= Ti(keep);
    Te= Te(keep);
    exp_to_exp= exp_to_exp(keep);
    
    Ti_by_Tt= Ti ./ (Ti + Te);
    BF= 60./ (exp_to_exp./ FR); % instantaneous breaths per minute

end % end function
