fullPath = 'C:\Users\Mark\Documents\EIT\Will_Ditcham\22-04-2021_#2\Will_-_Left_Laterel_Recumbency_STEM.mat';
data = load(fullPath);
data = data.data;
img = data.measurement.ZeroRef;

rl = data.patient.ROI.RightLung > 50;
ll = data.patient.ROI.LeftLung > 50;
bl = rl + ll;

plot(squeeze(sum(img, [1,2])));
p1 = 52;
p2 = 125;
p3 = 243;
img(img < 0) = 0;
slcs = img .* bl;
slcs = img;
tv = (img(:,:,p2) - img(:,:,p1)) .* bl;
tv = img(:,:,p2);
tv = tv .* bl;
% tv = img(:,:,p2) - (img(:,:,p1) + img(:,:,p3)) / 2;
sd = std(img(:,:,p1:p2), [], 3);

ee_time = [1.061738092 4.876130496 9.59496646 13.92056609 17.32206035 20.78254006 24.87219789,...
28.84388483 32.83523358 36.3153751 38.89107307 43.35430542 46.99174147 51.25835566,...
54.48289357 57.78607874 60.9516312 64.8839945 68.97365234 72.41447023 75.77664085,...
79.86629869 85.05701825 89.59889786 98.64333346 107.8843872 116.8501756 125.3440803] .* 50.86;
ee_time = round(ee_time);
slcs = slcs(:, :, ee_time);
slcs(slcs == 0) = nan;

subplot(2, 1, 1);
plot(squeeze(sum(img, [1,2]))); hold on;
for i = 1:length(ee_time)
    xline(ee_time(i));
end
subplot(2, 1, 2);
imagesc(mk_mosaic(-slcs, 0, [], i/2));
axis equal;
axis image;
axis off;
axis tight;
sgtitle('Will - Left Laterel Recumbency');

calc_colours('cmap_type', 'blue_red');
colormap(calc_colours('colourmap'));