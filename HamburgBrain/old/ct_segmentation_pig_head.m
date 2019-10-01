function [bounds] = ct_segmentation_pig_head() 
% Summary of this function goes here
% Detailed explanation goes here
	
% Get the image in...
dcm_dir = 'C:\Users\Mark Campbell\Documents\GraduateStudies\LAB\HamburgBrain\Models\Pig Skull CT\ScalarVolume_9';
[vol, spatial, dim] = dicomreadVolume(dcm_dir);

% Image Structure 
img = squeeze(vol);
layers = size(img,3); % To 512x512xlength or 1024x1024xlength
ext_bound    = {};
%% Segment the outer boundary and lungs
reverseStr = '';
found_empty= 0;
for i=1:layers
    loc = i- found_empty;
    
    % Progress message
    percentDone = 100 * i / layers; 
    msg = sprintf('Percent done: %3.1f', percentDone);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
    
    % Use the binarised image as the mask - threshold set based on background intensity
    I = img(:, :, i)*255; % Multiply pixel values of int16 to be on scale of 65535
    mask = imbinarize(I, 0.7); %imshowpair(I, mask, 'montage');
    
    % segment region
   
    bw = activecontour(I,mask,'edge',300); 
    
    % remove small external regions
    bw2 = bwareaopen(bw, 1000);
    bw3 = imerode(bw2,strel('disk',4));
    ext_skull_bounds = bwboundaries(bw3); % ext_skull_bounds = bwboundaries(bw3, 'no holes');
        
    ext_bound{loc}= [];
    for bnds= 1:length(ext_skull_bounds)
        if size( ext_skull_bounds{bnds}, 1) > 100
            ext_bound{loc}= cat(1, ext_bound{loc}, ext_skull_bounds{bnds});
        end % end if
    end % end for

    if isempty(ext_bound{loc})
        found_empty= found_empty+ 1;
    continue
    end % end if
end % end for

bounds.exterior = ext_bound;
bounds.length = loc;

end % end function


