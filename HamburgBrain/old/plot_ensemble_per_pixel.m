function plot_ensemble_per_pixel(seq, opt)
% This function was written by :
%                               Mark Campbell
%                               Carleton University
% opt can have fields:
%     plotLM
%     section
imgr= seq.imgr;
if ~exist('opt', 'var') || ~isfield(opt, 'plotLM')
    opt.plotLM= 1;
end % end if

if opt.plotLM~= 2
    intervals= diff(seq.eit.peaks);
    idx= seq.eit.peaks;
else
    intervals= diff(seq.eit.vals);
    idx= seq.eit.vals;
end % end if
% http://www.bg.ic.ac.uk/research/k.parker/guide_to_wia/03_ensemble_average.html
ensembleLength= round( 1.2* mean(intervals(intervals< 50)) );

if ~isfield(opt, 'section')
    opt.section= 0;
end % end if
% Which segment of time series at which to look
if opt.section== 1
    start= 1;
    lim= seq.eit.apn;
    figureTitle= seq.name+ " 1 - Ensemble for Brain Before Apnoea";
elseif opt.section== 2
    start= seq.eit.apn;
    lim= seq.eit.inj;
    figureTitle= seq.name+ " 2 - Ensemble for Brain After Apnoea and Before Injection";
elseif opt.section== 3
    start= seq.eit.inj;
    lim= seq.eit.vnt;
    figureTitle= seq.name+ " 3 - Ensemble for Brain After Injection and Before Apnoea End";
elseif opt.section== 4
    start= seq.eit.vnt;
    lim= size(imgr.elem_data, 2);
    figureTitle= seq.name+ " 4 - Ensemble for Brain After Apnoea End";
else
    start= 1;
    lim= size(imgr.elem_data, 2);
    figureTitle= seq.name+ " - Ensemble for Brain for Full Timeseries";
end % end if
idx= idx((idx> start) + (idx< lim)== 2);
n= 0;
slices= calc_slices(imgr);

for ispart= idx
    if ispart+ ensembleLength< lim
        if n== 0
            ensemble= slices(:, :, ispart: ispart+ ensembleLength);
        else
            ensemble= ensemble+ slices(:, :, ispart: ispart+ ensembleLength);
        end % end if
        n= n+1;
    end % end if
end % end for
ensemble= ensemble./ ensembleLength;
brainSeg= get_brain_segmentation_for_pig(seq);

sz_ensemble= size(ensemble);
ensemble= reshape(ensemble, sz_ensemble(1)* sz_ensemble(2), sz_ensemble(3)); 
roi_idx= [];
idx_ref= (brainSeg.upperL(2)-1)* sz_ensemble(1) + brainSeg.upperL(1)-1;
for col= 1: size(brainSeg.mtx, 2)
    ispart= find(brainSeg.mtx(:, col))+ sz_ensemble(1)*(col-1)+ idx_ref;
    roi_idx= horzcat(roi_idx, ispart');
end % end for
xax= 1:size(ensemble, 2);
figure('units','normalized','outerposition',[0 0 1 1]);
sgtitle(figureTitle);
axis equal;
axis square;
r= size(brainSeg.mtx, 1);
c= size(brainSeg.mtx, 2);
dat= ensemble(roi_idx, :);

j= 1;
mtx= brainSeg.mtx';
for i= 1:length(mtx(:))
    if mtx(i)== 1
        subplot(r, c, i);
        plot(xax, dat(j, :));
        j= j+ 1;
    end % end if
end % end for
cd 'C:\Users\Mark\Documents\GraduateStudies\LAB\HamburgBrain\Figures\10.2\valRef';
print_convert(char(figureTitle+ ".png"));
end % end function