function [bounds] = ct_segmentation_pig_thorax() 
% Summary of this function goes here
% Detailed explanation goes here
	
	% Get the image in...
	dcm_dir = 'C:\Users\Mark Campbell\Documents\GraduateStudies\LAB\HamburgBrain\Models\Pig Skull CT\ScalarVolume_9';
	[vol, spatial, dim] = dicomreadVolume(dcm_dir);

	% Image Structure 
	img = squeeze(vol);
	layers = size(img,3); % To 512x512xlength
	ext_bound    = {};
	l_lung_bound = {};
	r_lung_bound = {};
	c_lung_bound = {}; % Centre lung picks up the trachea or other internal air structures
	%% Segment the outer boundary and lungs
	for i=1:layers
	    loc = i;
	    I = img(:,:,loc);
	    % Use the binarised image as the mask - threshold set based on
	    % background intensity
	    mask = imbinarize(I);%,255*0.40);
	    % segment region
	    bw = activecontour(I,mask,'edge',300);
	    % remove small external regions
	    bw2 = bwareaopen(bw, 5000);
   		bw3 = imerode(bw2,strel('disk',4));
	    % ASSUMES LUNGS AND EXTERIOR ARE THREE LARGEST THINGS
		% if the trachea is largest it can get added as another object
	    ext_lung_bounds = bwboundaries(bw3);
	    n_ext_lung_bounds = size(ext_lung_bounds, 1);
	    ext_bound{i} = ext_lung_bounds{1};
	    for j=2:n_ext_lung_bounds
	        if mean(ext_lung_bounds{j}(:,2)) < 225
	            l_lung_bound{i} = ext_lung_bounds{j};
			else 
				l_lung_bound{i} = 0;
	        end
	        if mean(ext_lung_bounds{j}(:,2)) > 300
	            r_lung_bound{i} = ext_lung_bounds{j};
			else 
				r_lung_bound{i} = 0;
	        end
	        if (mean(ext_lung_bounds{j}(:,2)) > 225) && (mean(ext_lung_bounds{j}(:,2)) < 300)
	            c_lung_bound{i} = ext_lung_bounds{j};
			else 
				c_lung_bound{i} = 0;
	        end
	    end
	end

	%bounds.heart = heart_bound;
    bounds.left_lung = l_lung_bound;
    bounds.right_lung = r_lung_bound;
    bounds.exterior = ext_bound;
    bounds.esophagus = c_lung_bound;
    bounds.length = layers;

	keyboard
end


