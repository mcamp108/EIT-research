clear 
close all
path    = 'C:\Users\Mark\Dropbox\Joaquin\data\296439\';
contents = ls(path);

alpha = 0;
for i=3:size(contents,1)
    fName = strtrim(contents(i,:));
    if contains(fName, '.mat')
        fullPath= horzcat(path,fName);
        show_rotated_image(fullPath, alpha);
    end
end

% % update_rot_file(path);
% fullPath= horzcat(path,fileName);
% happy = false;
% alpha = 40;
% if ~happy
%     show_rotated_image(fullPath, alpha);
% else
%     rotate_ibex_images(path, alpha);
% end % end if
% % 
fullPath = 'C:\Users\Mark\Documents\EIT\Will_Ditcham\22-04-2021_#2\Will_-_Left_Laterel_Recumbency_STEM.mat';
data = load(fullPath);
data = data.data;
img = data.measurement.ZeroRef;





% data.measurement.CompositValue(1)
% img2 = movmedian(img, 50, 3);
% data.measurement.ZeroRef = img2;
% data.measurement.CompositValue = squeeze(sum(img2, [1,2]));
% save([fullPath]);
% 
% clf();
% hold on;

function update_rot_file(path)
    files = ls(path);
    ibexFiles = {};
    count=1;
    for i = 1:size(files,1)
        f = strtrim(files(i,:));
        if contains(f, '_ibex_rot.mat')
            ibexFiles{count} = f;
            count = count + 1;
        end % end if
    end % end for
    for i = 1:length(ibexFiles)
        f = ibexFiles{i};
        fullPath = horzcat(path,'\',f);
        data = load(fullPath);
        data = data.data;
        img = data.measurement.ZeroRef;
        data.measurement.CompositValue = squeeze(sum(img, [1,2]));
        save([fullPath(1:end-4),'_global.mat']);
    end
end

function ibexFiles = ibex_iterator(path)
    files = ls(path);
    ibexFiles = {};
    count=1;
    for i = 1:size(files,1)
        f = strtrim(files(i,:));
        if contains(f, '_ibex.mat')
            ibexFiles{count} = f;
            count = count + 1;
        end % end if
    end % end for
    for i = 1:length(files)
        
    end
end % end function


function imgRotated = rotate_img(img, angle)
    % can also set angle to 'flip-vert' or 'flip-horz' to flip about the
    % vertical and horizontal axes.
    if isstr(angle)
        if strcmp(angle, 'flip-vert')
            for i = 1: length(img) 
                imgRotated(:,:,i) = fliplr(img(:,:,i));
            end % end for
        elseif strcmp(angle, 'flip-horz')
            for i = 1: length(img) 
                imgRotated(:,:,i) = flipud(img(:,:,i));
            end % end for
        end % end if
    else
        imgScaled = imresize(img,4);
        for i = 1: length(imgScaled) 
            scaledRotated(:,:,i)= imrotate(imgScaled(:,:,i),angle);
        end % end for
        % measure how dimensions have changed
        oldNonZeroCol   = sum(img(:,:,1), 1) > 0;
        oldNonZeroRow   = sum(img(:,:,1), 2) > 0;
        newNonZeroCol   = sum(scaledRotated(:,:,1), 1) > 0;
        newNonZeroRow   = sum(scaledRotated(:,:,1), 2) > 0;
        % stretch image to match original body dimensions
        imgStretch  = scaledRotated(newNonZeroRow,newNonZeroCol,:);
        scaledRotated = imresize(imgStretch, [sum(oldNonZeroRow),sum(oldNonZeroCol)]);
        imgRotated  = img;
        imgRotated(oldNonZeroRow,oldNonZeroCol,:) = scaledRotated;        
    end % end if
end % end function

function show_rotated_image(fullPath, angle)
    fileName = strsplit(fullPath, '\');
    fileName = fileName{end};
    data = load(fullPath);
    data = data.data;
    if nargin == 1
        angle = 0;
    end % end if
    img = data.measurement.ZeroRef;
    imgRot = rotate_img(img, angle);
    figure('units','normalized','outerposition',[0.5 0 .5 1]);
    subplot(2,1,1); imagesc(std(img(:,:,:),1,3));   axis equal; title('Original');
    stdImgRot = std(imgRot(:,:,:),1,3);
    subplot(2,1,2); imagesc(stdImgRot); axis equal; title('Rotated');
end % end function

function rotate_ibex_images(path, angle)
    files = ibex_iterator(path);
    for i = 1:length(files)
        f = files{i};
        fullPath = horzcat(path,'\',f);
        load(fullPath);
        img = data.measurement.ZeroRef;
        imgRot = rotate_img(img, angle);    
        data.measurement.ZeroRef = imgRot;
        data.measurement.CompositValue = squeeze(sum(img, [1,2]));
        save([fullPath(1:end-4),'_rot.mat']);
    end % end for
end % end function


