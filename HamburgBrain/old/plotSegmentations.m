function plot_segmentations(imgr, brainSeg, timePoints, acvOpt)
% This function was written by :
%                               Mark Campbell
%                               Carleton University
% Parameters:
%    imgr: Inverse solution object
% 
%    brainSeg: String of model name whose segmentations will be plotted.
%    Currently accepted are: "mri-ct" and "pighead"
%
% timePoints:
%     includes stopVent, inject, stopVent relative to first frame included
%     in the time series. Includes the framerate of data acquisition.
% acvOpt:
%     Options:
%       mode: "show each", "show_mean", "vs mean", "show ROC", "each mean
%       vs all mean"
%       title: Title of figure
%       roiTitle: Name of ROI to display in figure legend
%       pixels: The indices of pixels in the segmentation to be plotted
% Returns:
%     A plot containing the voltage data for each pixel in the segmentation
%     for imgr.

segs= getSegmentationsForModel(brainSeg);
slices= calc_slices(imgr);    

fsTrue= timePoints.fsTrue;
stopVent= timePoints.stopVent /fsTrue;
inject= timePoints.inject /fsTrue;
startVent= timePoints.startVent /fsTrue;
if isfield(timePoints, 'start') && isfield(timePoints, 'stop')
    slices= slices(:, :, timePoints.start: timePoints.stop);
end % end if

segsAcv= get_acv(slices, segs.brain.upperL, segs.brain.mat, acvOpt.mode);
if acvOpt.mode== "each mean vs all mean"
    segs.brain.mat( find(segs.brain.mat) )= segsAcv;
    figure; imagesc(segs.brain.mat);
    title("Mean value of pixel over timeframe vs mean value of all pixels over timeframe");
else
    acvMin= min(segsAcv(:));
    avcMax= max(segsAcv(:));
    if ~isfield(acvOpt, 'pixels')
        pixels= 1:size(segsAcv, 1);
    else
        pixels= acvOpt.pixels;
    end % end if
    xax = (1: size(segsAcv, 2))/fsTrue;

    figure; hold on;
    plot([stopVent, stopVent], [acvMin, avcMax]); 
    plot([inject, inject], [acvMin, avcMax]); 
    plot([startVent, startVent], [acvMin, avcMax]); 
    plot(xax, segsAcv(pixels, :));
    title(acvOpt.title); legend('Stop Ventilation' , 'Inject Saline', 'Start Ventilation'); hold off;
    xlabel("Time (s)");
    if acvOpt.mode== "show ROC"
        ylabel("Difference (Xn+1 - Xn)");
    else
        ylabel("Voltage");
    end % end if
end % end if
end % end function

function segs= getSegmentationsForModel(model)
% returns: 
%           brainUpperL: The upper left reference corner for segmentatin of the
%           brain from reconstructed 64 x 64 image of airchill pig.
%           
%           brainSeg: The logical matrix representing the pixels in the 64
%           x 64 image that are included in the brain, referenced from
%           upperL.

segs= struct;
if model== "mri-ct"
    segs.brain.upperL= [43, 20]; % row, column
    segs.brain.mat= [
                      0 0 0 0 1 1 1 0 0 0 0 0;
                      0 0 1 1 1 1 1 1 1 0 0 0;
                      0 1 1 1 1 1 1 1 1 1 0 0;
                      1 1 1 1 1 1 1 1 1 1 1 0;
                      1 1 1 1 1 1 1 1 1 1 1 1;
                      1 1 1 1 1 1 1 1 1 1 1 1;
                      0 1 1 1 1 1 1 1 1 1 1 0;
                      0 0 1 1 1 1 1 1 1 1 0 0;
                      0 0 0 1 1 1 1 1 1 0 0 0;
                      0 0 0 0 1 1 1 0 0 0 0 0;
                                               ];
elseif model== "pighead"
    r_hem_upprL=       {[], [16,25], [], [22,21], [8,14],  [7,15], [], [6,20], [], [22,23], [25,19]};
    l_hem_upprL=       {[], [4,22],  [], [24,18], [9,10],  [9,9],  [], [6,13], [], [25,19], [25,14]};
    brainstem_upprL= {[], [12,19], [], [20,14], [16,16], [3,12], [], [19,14],[], [],      [19,21]};

    r_hem_seg=       {[],... 
                      [0 1 1 1 1 0;1 1 1 1 1 1;1 1 1 1 1 1;0 1 1 1 1 1;0 0 1 1 1 1;0 0 0 1 1 1],... 
                      [],... 
                      [0 1 1 1 1 0 0 0;1 1 1 1 1 1 0 0;0 1 1 1 1 1 1 0;0 0 1 1 1 1 1 1;0 0 0 1 1 1 0 0;0 0 0 0 1 0 0 0],... 
                      [1 1 1 1 1 1 0;0 1 1 1 1 1 1;0 1 1 1 1 1 1;0 0 1 1 1 1 0;0 0 1 0 0 0 0],... 
                      [0 0 0 1 1 1;0 1 1 1 1 1;1 1 1 1 1 1;0 1 1 1 1 1;0 1 1 1 1 1],... 
                      [],... 
                      [1 1 1 1 1;1 1 1 1 1;1 1 1 1 1;1 1 1 1 1],... 
                      [],... 
                      [0 0 1 1 0 0; 0 1 1 1 1 0; 1 1 1 1 1 1; 0 1 1 1 1 0; 0 0 1 1 0 0],...  
                      [1 1 1 1 1;1 1 1 1 1;1 1 1 1 1;1 1 1 1 1;1 1 1 1 1]};
    l_hem_seg=       {[],... 
                      [1 1 0 0 0 0;1 1 1 1 0 0; 1 1 1 1 1 0;1 1 1 1 1 1;1 1 1 1 1 1;0 1 1 1 1 1],... 
                      [], ...
                      [0 0 1 0 0 0 0;0 1 1 1 0 0 0;1 1 1 1 1 0 0;1 1 1 1 1 1 0;1 1 1 1 1 1 1;0 1 1 1 1 0 0],...
                      [0 0 1 1 0;0 1 1 1 0;1 1 1 1 1;1 1 1 1 1;1 1 1 1 0;1 1 1 0 0],...
                      [0 0 1 1 1 0;0 1 1 1 1 1;1 1 1 1 1 1;1 1 1 1 1 1;0 1 1 1 1 0],...
                      [],...
                      [1 1 1 1 1;1 1 1 1 1;1 1 1 1 1;1 1 1 1 1],...
                      [],...
                      [0 0 1 1 0 0;0 1 1 1 1 0;1 1 1 1 1 1;0 1 1 1 1 0;0 0 1 1 0 0],...
                      [1 1 1 1 1;1 1 1 1 1;1 1 1 1 1;1 1 1 1 1;1 1 1 1 1]};
    brainstem_seg= {[],...
                      [1 1 1 1 0;1 1 1 1 1;1 1 1 1 1;1 1 1 1 0],...
                      [],...
                      [1 1 1 1 1;1 1 1 1 1;1 1 1 1 1;1 1 1 1 1],...
                      [1 1 1;1 1 1;1 1 1],...
                      [1 1 1;1 1 1;1 1 1],...
                      [],...
                      [1 1 1 1 1;1 1 1 1 1;1 1 1 1 1],...
                      [],...
                      [],...
                      [1 1 1 1 1;1 1 1 1 1;1 1 1 1 1;1 1 1 1 1]};
    for i= 1: length(r_hem_upprL)
        segs(i).r_hem.upperL= r_hem_upprL{i};
        segs(i).r_hem.mat= r_hem_seg{i};
        segs(i).l_hem.upperL= l_hem_upprL{i};
        segs(i).l_hem.mat= l_hem_seg{i};
        segs(i).brainstem.upperL= brainstem_upprL{i};
        segs(i).brainstem.mat= brainstem_seg{i};
    end % end for
end % end if

end % end function

function acv= get_acv(slices, upperL, mat, Mode)
% Get time series of average conductivity value for a region of interest
% upperL is the row and column of the pixel in the upper left corner of the
% rectangle that encompasses the ROI.
if ~exist('Mode', 'var')
    Mode= "show mean";
end % end if

do_function= true;
if isempty(upperL) || isempty(mat)
    do_function= false;
    acv= [];
end % end if
if do_function
    sz_slices= size(slices);
    idx_ref= (upperL(2)-1)* sz_slices(1) + upperL(1)-1;
    slices= reshape(real(slices), sz_slices(1)* sz_slices(2), sz_slices(3)); 
    roi_idx= [];
    for col= 1: size(mat, 2)
        ispart= find(mat(:, col))+ sz_slices(1)*(col-1)+ idx_ref;
        roi_idx= horzcat(roi_idx, ispart');
    end % end for
    
    acv= slices(roi_idx, :);
    
    if Mode== "show each"
    elseif Mode== "vs mean"
        meanAcv= mean(slices(roi_idx, :), 1);
        acv= acv- meanAcv;
    elseif Mode== "show ROC"
        acv= diff(acv);
    elseif Mode== "each mean vs all mean"
        acv= mean(acv, 2);
        acv= acv- mean(acv);
    else
        acv= mean(acv, 1);
    end % end if
end % end if
end % end function