function mk_vid(seq, start, stop, suffix)

if ~exist('start', 'var')
    start= 1;
end % end if
if ~exist('stop', 'var')
    stop= size(seq.eit.data, 2);
end % end if

% Make video of reconstructed image over time
cd vid
vidfile = VideoWriter(seq.name+ suffix+ ".mp4",'MPEG-4');
vidfile.FrameRate= round(seq.eit.fs);
brainSeg= get_brain_segmentation_for_pig(seq);

imgr2= seq.imgr; % used as imgr for desired time period
imgr2.elem_data= seq.imgr.elem_data(:, start:stop);
clim= max(imgr2.elem_data(:));
cmin= min(imgr2.elem_data(:));

imgr3= imgr2; % used as copy for video reel
imgr3.calc_colours.ref_level= 0;
imgr3.calc_colours.clim= clim;

showSegmentation= imgr3; % Show elem data from only segmentation region
showSegmentation.elem_data= zeros(size(imgr3.elem_data))- 50; % set low to have contrast
showSegmentation.elem_data(brainSeg.idx, :)= imgr3.elem_data(brainSeg.idx, :); % transfer over only brain region element data
onlyBrain= showSegmentation; % make copy for video reel

topL= sum(imgr2.elem_data(brainSeg.topL, :), 1);
topR= sum(imgr2.elem_data(brainSeg.topR, :), 1);
botL= sum(imgr2.elem_data(brainSeg.botL, :), 1);
botR= sum(imgr2.elem_data(brainSeg.botR, :), 1);
holdingCell= [topL, topR, botL, botR];
highest= max(holdingCell); % values for plotting progress line
lowest= min(holdingCell);

figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2, 2, [3,4]);
xax= linspace(1, length(topL)/seq.eit.fs, length(topL));
hold on
plot(xax, topL);
plot(xax, botL);
plot(xax, topR);
plot(xax, botR);
% legend("top left", "bottom left", "top right", "bottom right", "Frame");
hold off
xlabel("Time (s)");
ylabel("Voltage (uV)");
open(vidfile);
for frame= 1:(stop- start)
    imgr3.elem_data= imgr2.elem_data(:, frame);
    onlyBrain.elem_data= showSegmentation.elem_data(:, frame);
    subplot(2, 2, 1)
    show_fem(imgr3);
    eidors_colourbar(imgr3); caxis([cmin clim]);
    title(seq.name);
    
    subplot(2, 2, 2)
    show_fem(onlyBrain);
    eidors_colourbar(imgr3); caxis([cmin clim]);
    title(seq.name + " brain segmentation");
    
    subplot(2, 2, [3, 4])
    hold on
    progress= plot([xax(frame), xax(frame)], [highest, lowest], 'g');
    hold off
    legend("top left", "bottom left", "top right", "bottom right", "Frame");
    title(seq.name + " regional conductivity");
    writeVideo(vidfile,getframe(gcf));
    delete(progress);
end % end for
close(vidfile);
cd ../

end % end function

% function
% 
% end