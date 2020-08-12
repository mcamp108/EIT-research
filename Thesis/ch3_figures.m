pig = '11-2';
[fmdl, imdl] = get_pig_mdl(pig);
% Load data
ref = 'self';
D = load_HamburgBrain_data(pig, ref);

%%
seq = D.seq4;
% close all;
% bigFig();
clf();
subplot(2, 1, 1);
    plot_perf_data(seq);
    xlim([25,29.75]);
subplot(2, 1, 2);
    plot_sum_Z(seq);
    xlim([25,29.75]);
saveas( gcf, 'C:\Users\Mark\Documents\GraduateStudies\LAB\Thesis\files\figures\Chap_3\methods11-2seq4.svg');

%%
pig='8-2';
elec_z_plane=106;
cd(horzcat('C:\Users\Mark\Documents\GraduateStudies\LAB\HamburgBrain\Models\', pig,'\mesh'));

seg3dOutMatFile = horzcat(pig, '_seg3D_out');
load(seg3dOutMatFile);
V = scirunnrrd.data;
V = rot90(V, 1); % looking through front of pig
colormap(flipud(parula()))
elecPlaneImg = squeeze( V(:, :, elec_z_plane) );
img1 = elecPlaneImg;

rowsum = sum(elecPlaneImg, 2);
colsum = sum(elecPlaneImg, 1);
elecPlaneImg = elecPlaneImg(rowsum>0, colsum>0);
img2 = elecPlaneImg;

elecPlaneImg(elecPlaneImg~=3) = 0;
img3 = elecPlaneImg;

seg = imresize(elecPlaneImg, [64 64]);
seg(seg>0) = 1;
img4 = seg;

    % find brain centroid
rows = find(sum(seg, 2)>0);
cols = find(sum(seg, 1)>0);
cX = mean(cols);
cY = mean(rows);
r=40;
cnt = 0;

for i = -90:60:210
    cnt = cnt + 1;
    t = linspace(-i,-i-60,128);
    x = [cX, cX+r*cosd(t), cX];
    y = [cY, cY+r*sind(t), cY];
    bw = poly2mask( x, y, size(seg,1),size(seg,2));
    seg(bw) = seg(bw)*cnt;
end
img5 = seg;

subplot(1,5,1);
    plot_img(img1);
subplot(1,5,2);
    plot_img(img2);
subplot(1,5,3);
    plot_img(img3);
subplot(1,5,4);
    plot_img(img4);
subplot(1,5,5);
    plot_img(img5);
fg = gcf();

fg.Children(5).Colormap = [ .9 .9 .9;...
    1 0 0;...
    1 1 165/255;...
    108/255 0 212/255;...
    0 0 1];
fg.Children(4).Colormap = fg.Children(5).Colormap;
fg.Children(3).Colormap = [ .9 .9 .9;...
    108/255 0 212/255];
fg.Children(2).Colormap = fg.Children(3).Colormap;
fg.Children(1).Colormap = [ .9 .9 .9;...
    0.8500 0.3250 0.0980;...
    0.9290 0.6940 0.1250;...
    0.4940 0.1840 0.5560;...
    0.4660 0.6740 0.1880;...
    0.3010 0.7450 0.9330;...
    0.6350 0.0780 0.1840];
% saveas( gcf, 'C:\Users\Mark\Documents\GraduateStudies\LAB\Thesis\files\figures\Chap_3\methodsBrainSeg.svg');
%%

function plot_img(img)
img=flipud(img);
imagesc('XData',linspace(0,1,size(img,2)),'YData',linspace(0,1,size(img,1)),'CData',img);
axis equal;axis image; axis off; xlim([0,1]);ylim([0,1]);
end % end function
% ======================================================================= %
function plot_perf_data(seq)
    bot = 5;
    fs = seq.perf.tickrate;
    x = seq.perf.data(seq.perf.apn:seq.perf.inj);
    xax = (1: length(x))/fs;
    plot(xax, x,'LineWidth',2);
    xlabel("Time (s)",'FontSize',16);
    ylabel("Arterial Pressure (mmHg)",'FontSize',16);
%     title("Arterial Blood Pressure: "+ seq.name);
    use = ((seq.perf.peaks >= seq.perf.apn) + (seq.perf.peaks <= seq.perf.inj) == 2);
    seq.perf.peaks = seq.perf.peaks(use) - seq.perf.apn+ 1;
    seq.perf.vals = seq.perf.vals(:, use) - seq.perf.apn + 1;
    hold on;
    for i = 1:length(seq.perf.peaks)
        loc = seq.perf.peaks(i);
        y = x(loc);
        plot([loc/ fs, loc/ fs], [y-bot, 120], '--r','LineWidth',1);
    end % end for
    vals = unique( seq.perf.vals(:) );
    for i = 1:length(vals)
        loc = vals(i);
        y = x(loc);
        plot([loc/ fs, loc/ fs], [y-bot/2, 120], 'g','LineWidth',1);
    end % end for
    hold off;
end % end function
% ======================================================================= %
function plot_sum_Z(seq)
    bot = 0.001;
    upper = -0.005;
    lower = -0.02;
    fs = seq.eit.fs;
    x = sum(seq.imgr.elem_data, 1);
    x = x(seq.eit.apn:seq.eit.inj);
    xax = (1: size(x, 2))/fs;
    plot(xax, x,'LineWidth',2);
    xlabel("Time (s)",'FontSize',16);
    ylabel("EIT Global Signal (Delta Z)",'FontSize',16);
%     title("Global Signal: "+ seq.name);
    use = ((seq.eit.peaks >= seq.eit.apn) + (seq.eit.peaks <= seq.eit.inj) == 2);
    seq.eit.peaks = seq.eit.peaks(use) - seq.eit.apn+ 1;
    seq.eit.vals = seq.eit.vals(:, use) - seq.eit.apn + 1;
    hold on;
    
    vals = unique( seq.eit.vals(:) );
    for i = 1:length(vals)
        loc = vals(i);
        y = x(loc);
        plot([loc/fs, loc/fs], [y+bot, upper], 'g','LineWidth',1);
        plot([loc/fs, loc/fs], [y-bot, lower+bot], '--k','LineWidth',1);
        plot([loc/fs, loc/fs], [lower+0.0005, lower-0.0005], 'k','LineWidth',1);
        if i < length(vals)
            plot([vals(i)/fs, vals(i+1)/fs],[lower,lower], 'k','LineWidth',2);
        end
    end % end for
    for i = 1:length(seq.eit.peaks)
        loc = seq.eit.peaks(i);
        y = x(loc);
        plot([loc/fs, loc/fs], [y+bot, upper], '--r','LineWidth',1);
        plot([loc/fs, loc/fs], [y-bot, lower+bot], '--k','LineWidth',1); % line to systole in frame line
        plot(loc/fs, lower, 'xr','MarkerSize',16);
    end % end for

    
    hold off;
end % end function
% ======================================================================= %
function bigFig()
    figure('units','normalized','outerposition',[0 0 1 1]);
end % end function