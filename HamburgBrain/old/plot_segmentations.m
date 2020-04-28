function roi= plot_segmentations(seq, opt)
% This function was written by :
%                               Mark Campbell
%                               Carleton University
% Returns: Plots the ensemble-averaged conductivity values for the
% specified pixels over the specified timeframe
%
% Parameters:
%    seq: Pig data sequence
%
%   opt:
%     Options:
%       mode: "show each", "show_mean", "vs mean", "show ROC", "each mean
%       vs all mean"
%       title: Title of figure
%       roiTitle: Name of ROI to display in figure legend
%       pixels: The indices of pixels in the segmentation to be plotted
%       start
%       stop
%       plotLM: 
%             1: (plot perfusion peaks) 
%             2: (plot perfusion valleys)
%             3: (plot perfusion peaks and valleys)

if ~exist('opt', 'var')
    opt.mode= "show each";
end % end if
if ~isfield(opt, 'mode')
    opt.mode= "show each";
end % end if
if ~isfield(opt, 'start')
    opt.start= 1;
end % end if
if ~isfield(opt, 'stop')
    opt.stop= size(seq.imgr.elem_data, 2);
end % end if
opt.sel= [0 0 0 1];

seq2= seq;
brainSeg= get_brain_segmentation_for_pig(seq);
roi= brainSeg.idx;
seq2.imgr.elem_data= seq.imgr.elem_data(roi,:);
ensemble= get_ensembles(seq2, 1, opt);
roi= ensemble.imgr_ensemble_t;
plot(roi');


% 
% if opt.mode== "each mean vs all mean"
%     brainSeg.mat( brainSeg.mat==1 )= segsAcv;
%     figure; imagesc(brainSeg.mat);
%     title("Mean value of pixel over timeframe vs mean value of all pixels over timeframe");
% else
%     acvMin= min(segsAcv(:));
%     acvMax= max(segsAcv(:));
%     % If only certain pixels of the segmentation are to be plotted
%     if ~isfield(opt, 'pixels')
%         pixels= 1:size(segsAcv, 1);
%     else
%         pixels= opt.pixels;
%     end % end if
%     xax = (1: size(segsAcv, 2))/seq.eit.fs;
% 
%     figure; hold on;
%     plot([apn, apn], [acvMin, acvMax]); 
%     plot([inj, inj], [acvMin, acvMax]); 
%     plot([vnt, vnt], [acvMin, acvMax]); 
%     plot(xax, segsAcv(pixels, :));
%     title("Brain Segmentation: "+ seq.name);     
%     xlabel("Time (s)");
%     ylabel("Conductivity (Au)");
%     hold off
%     % if rate of change of each pixel is to be plotted
%     if opt.mode== "show ROC"
%         ylabel("Difference (Xn+1 - Xn)");
%     else
%         ylabel("Voltage (mV)");
%     end % end if
% end % end if
% if isfield(opt, 'plotLM')
%     if opt.plotLM < 4 && opt.plotLM > 0
%         plot_perf_landmarks(seq, acvMin, acvMax, opt.plotLM);
%     else
%         disp("Unrecognized landmark plotting code.");
%     end % end if
% end % end if
% 
% legend('Start of Apnoea' , 'Saline Injected', 'End of Apnoea');
% end % end function
% 
% function acv= get_acv(slices, brainSeg, Mode)
% % Get time series of average conductivity value for a region of interest
% % upperL is the row and column of the pixel in the upper left corner of the
% % rectangle that encompasses the ROI.
% if ~exist('Mode', 'var')
%     Mode= "show mean";
% end % end if
% do_function= true;
% upperL= brainSeg.upperL;
% mat= brainSeg.mat;
% 
% if isempty(upperL) || isempty(mat)
%     do_function= false;
%     acv= [];
% end % end if
% 
% if do_function
%     sz_slices= size(slices);
% %     roi= slices(
%     idx_ref= (upperL(2)-1)* sz_slices(1) + upperL(1)-1;
%     slices= reshape(real(slices), sz_slices(1)* sz_slices(2), sz_slices(3)); 
%     roi_idx= [];
%     for col= 1: size(mat, 2)
%         ispart= find(mat(:, col))+ sz_slices(1)*(col-1)+ idx_ref;
%         roi_idx= horzcat(roi_idx, ispart');
%     end % end for
%     
%     acv= slices(roi_idx, :);
%     
%     if Mode== "show each"
%     elseif Mode== "vs mean"
%         meanAcv= mean(slices(roi_idx, :), 1);
%         acv= acv- meanAcv;
%     elseif Mode== "show ROC"
%         acv= diff(acv);
%     elseif Mode== "each mean vs all mean"
%         acv= mean(acv, 2);
%         acv= acv- mean(acv);
%     else
%         acv= mean(acv, 1);
%     end % end if
% end % end if
% end % end function
% 
% function plot_perf_landmarks(seq, acvMin, acvMax, pOpt)
% mid= (acvMin+ acvMax)/ 2;
% hold on
% if pOpt== 1 || 3
%     for i= 1:length(seq.eit.peaks)
%         loc= seq.eit.peaks(i);
%         y= seq.eit.data(loc);
%         plot([loc/ seq.eit.fs, loc/ seq.eit.fs], [acvMax, mid], 'r');
%     end % end for
% end % end if
% if pOpt== 2 || 3
%     for i= 1:length(seq.eit.vals)
%         loc= seq.eit.vals(i);
%         y= seq.eit.data(loc);
%         plot([loc/ seq.eit.fs, loc/ seq.eit.fs], [acvMin, mid], 'g');
%     end % end for
% end % end if
% hold off
% end % end function