function [fmdl, mat_idx] = mk_3D_human_model(num_elecs, rows_elecs) 
% Using many equations for lung, esophagus and heart boundaries 
% Generated from CT segmentation using a watershed algorithm 
% Create a 3D EIT inverse and forward model.

	%% STEP 1 - get the segmented boundaries for each row
	% adapted from image processing course:
	% This could take a while - save these as a separate file one done
	if exist('ct_boundaries.mat', 'file') == 2
		bounds = load('ct_boundaries.mat');
		bounds = bounds.bounds;
	else
		bounds = ct_segmentation(); % Only run if we have to... Takes ages!!
		save('ct_boundaries.mat','bounds'); % Save it to save time
	end
	% Fix the heart region... Just remove the high mistakes from segmentation
	bounds.heart(1,1:185) = {0};
	bounds.heart(1,300:end) = {0};
	bounds.right_lung(1,340:end) = {0};
	% For each layer downsample the plots and shift to be centered on 0
	heart  = [];
	l_lung = [];
	r_lung = [];
	ext    = [];
	esoph  = [];
	n = 10; % Keep every nth element only
	offset = 512/2; % Centre on 0
	scale = 512/2/10; % Scale between -10 and 10 % TODO 0 to 1
	vert_scale = bounds.length/15;
	bounds.left_lung(356:384) = {[]}; % Segmentation output workaround
	bounds.esophagus(333:384) = {[]}; % Segmentation output workaround
	for i=1:bounds.length % For each layer
		j = 384-i;
		temp = bounds.heart{i};
		if temp ~= 0
			temp = (temp-offset)./scale;
			%bounds.heart{i}      = temp(1:n:end,:);
			temp = temp(1:n:end,:);
			temp(:,3) = j/vert_scale; % some vertical scaling!!
			heart = [heart; temp];
		end
		temp = bounds.left_lung{i};
		if temp ~= 0
			temp = (temp-offset)./scale;
			%bounds.left_lung{i}  = temp(1:n:end,:);
			temp = temp(1:n:end,:);
			temp(:,3) = j/vert_scale; % some vertical scaling!!
			l_lung = [l_lung; temp];
		end
		temp = bounds.right_lung{i};
		if temp ~= 0
			temp = (temp-offset)./scale;
			%bounds.right_lung{i} = temp(1:n:end,:);
			temp = temp(1:n:end,:);
			temp(:,3) = j/vert_scale; % some vertical scaling!!
			r_lung = [r_lung; temp];
		end
		temp = bounds.exterior{i};
		if temp ~= 0
			temp = (temp-offset)./scale;
			%bounds.exterior{i}   = temp(1:n:end,:);
			temp = temp(1:n:end,:);
			temp(:,3) = j/vert_scale; % some vertical scaling!!
			ext = [ext; temp];
		end
		temp = bounds.esophagus{i};
		if temp ~= 0
			temp = (temp-offset)./scale;
			%bounds.esophagus{i}  = temp(1:n:end,:);
			temp = temp(1:n:end,:);
			temp(:,3) = j/vert_scale; % some vertical scaling!!
			esoph = [esoph; temp];
		end
	end

	if 0 % Do a nice nice plot
		figure(1)
		clf;
		scatter3(ext(:,1),ext(:,2),ext(:,3),'.')
		hold on
		scatter3(heart(:,1),heart(:,2),heart(:,3),'.')
		scatter3(l_lung(:,1),l_lung(:,2),l_lung(:,3),'.')
		scatter3(r_lung(:,1),r_lung(:,2),r_lung(:,3),'.')
		scatter3(esoph(:,1),esoph(:,2),esoph(:,3),'.')
	end
	%% STEP 2 - Save pointcloud to .XYZ files
	dlmwrite('heart.XYZ',heart,'delimiter',' ','precision',5)
	dlmwrite('left_lung.XYZ',l_lung,'delimiter',' ','precision',5)
	dlmwrite('right_lung.XYZ',r_lung,'delimiter',' ','precision',5)
	dlmwrite('external.XYZ',ext,'delimiter',' ','precision',5)
	dlmwrite('esophagus.XYZ',esoph,'delimiter',' ','precision',5)
end

