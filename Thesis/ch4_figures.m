% This function was written by :
%                               Mark Campbell
%                               Carleton University
[fmdl, imdl] = get_hum_head_mdl(); % load models
img = mk_image(fmdl, 0.41); % Background conductivity is scalp
img.elem_data([fmdl.mat_idx{1}]) = 0.41;    %   1: scalp
img.elem_data([fmdl.mat_idx{2}]) = 0.016;   %   2: skull
img.elem_data([fmdl.mat_idx{3}]) = 1.71;    %   3: CSF
img.elem_data([fmdl.mat_idx{4}]) = 0.47;    %   4: grey matter
img.elem_data([fmdl.mat_idx{5}]) = 0.22;    %   5: white matter
img.elem_data([fmdl.mat_idx{6}]) = 0.7;     %   6: diploe
img.elem_data([fmdl.mat_idx{7}]) = 0.0001;  %   7: air
img.fwd_solve.get_all_meas = 1;
vh = fwd_solve(img); % reference data (assumes no noise)
J = calc_jacobian( calc_jacobian_bkgnd( imdl) );
iRtR = inv(prior_noser( imdl ));
hp = 0.17;
iRN = hp^2 * speye(size(J,1));
RM = iRtR*J'/(J*iRtR*J' + iRN);
imdl.solve = @solve_use_matrix; 
imdl.solve_use_matrix.RM  = RM;
figure('units','normalized','outerposition',[0 0 1 1]);

% parameters for setting up brain foci size and position.
% constants
TIMES   = 10;
for DELTA = [0.9, 0.99, 1.01, 1.1]
    run_sims(DELTA, TIMES, fmdl, imdl, img, vh)
end

%%

% Simulation Plots
simDir = 'C:\Users\Mark\Documents\GraduateStudies\LAB\EIT-neuroimaging\Figures\simulations\';
% foci distance (3)    radial position (8)    foci size (3)
df1 = readmatrix(sprintf('%slocErr0.90Delta', simDir));
df2 = readmatrix(sprintf('%slocErr0.99Delta', simDir));
df3 = readmatrix(sprintf('%slocErr1.01Delta', simDir));
df4 = readmatrix(sprintf('%slocErr1.10Delta', simDir));
%%
sim_plots(df1,1);
cd 'C:\Users\Mark\Documents\GraduateStudies\LAB\Thesis\files\figures\Chap_4';
sgtitle('Localization Accuracy for 10% Conductivity Decrease', 'FontSize', 20);
printPDF("10decrease");
sim_plots(df2,2);
sgtitle('Localization Accuracy for 1% Conductivity Decrease', 'FontSize', 20);
printPDF("1decrease");
sim_plots(df3,3);
sgtitle('Localization Accuracy for 1% Conductivity Increase', 'FontSize', 20);
printPDF("1increase");
sim_plots(df4,4);
sgtitle('Localization Accuracy for 10% Conductivity Increase', 'FontSize', 20);
printPDF("10increase");

%%
function sim_plots(df1, n)
figure(n);clf();
radPos = linspace(-pi, pi, 9);
radPos = round(radPos(1:8), 2); % 8 positions

d1 = df1(1:24,:);
d2 = df1(25:48,:);
d3 = df1(49:72,:);

d1_r1to8_s1 = d1(1:3:24, 1);
d1_r1to8_s2 = d1(2:3:24, 1);
d1_r1to8_s3 = d1(3:3:24, 1);

d2_r1to8_s1 = d2(1:3:24, 1);
d2_r1to8_s2 = d2(2:3:24, 1);
d2_r1to8_s3 = d2(3:3:24, 1);

d3_r1to8_s1 = d3(1:3:24, 1);
d3_r1to8_s2 = d3(2:3:24, 1);
d3_r1to8_s3 = d3(3:3:24, 1);

subplot(3,1,1);
    plot(d1_r1to8_s1, '-o', 'LineWidth',2,'MarkerSize',8); ylim([0 1]); hold on
    plot(d1_r1to8_s2, '-x', 'LineWidth',2,'MarkerSize',8); ylim([0 1]);
    plot(d1_r1to8_s3, '-v', 'LineWidth',2,'MarkerSize',8); ylim([0 1]);
    ax=gca;ax.XTickLabels=radPos;
    title('Foci Distance From Brain Centre = 1/8 Brain Radius', 'FontSize', 18);
subplot(3,1,2);
    plot(d2_r1to8_s1, '-o', 'LineWidth',2,'MarkerSize',8); ylim([0 1]); hold on
    plot(d2_r1to8_s2, '-x', 'LineWidth',2,'MarkerSize',8); ylim([0 1]);
    plot(d2_r1to8_s3, '-v', 'LineWidth',2,'MarkerSize',8); ylim([0 1]);
    ax=gca;ax.XTickLabels=radPos;
    title('Foci Distance From Brain Centre = 1/4 Brain Radius', 'FontSize', 18);
    ylabel('Localization Error, X-Y Plane (% Brain Radius)', 'FontSize', 20);
subplot(3,1,3);
    plot(d3_r1to8_s1, '-o', 'LineWidth',2,'MarkerSize',8); ylim([0 1]); hold on
    plot(d3_r1to8_s2, '-x', 'LineWidth',2,'MarkerSize',8); ylim([0 1]);
    plot(d3_r1to8_s3, '-v', 'LineWidth',2,'MarkerSize',8); ylim([0 1]);
    ax=gca;ax.XTickLabels=radPos;
    title('Foci Distance From Brain Centre = 1/2 Brain Radius', 'FontSize', 18);
    xlabel('Radial Position from Baseline', 'FontSize', 20);
end

function run_sims(DELTA, TIMES, fmdl,imdl, img, vh)

SNR     = 80; % measurement signal to noise ratio. ref for SNR: https://www.ncbi.nlm.nih.gov/books/NBK549564/
IMGCOLS = 8;
SAVEDIR = 'C:\Users\Mark\Documents\GraduateStudies\LAB\EIT-neuroimaging\Figures\simulations\';
targetMat   = 4;
matRef      = 0.47;
mdlLims     = [min(fmdl.nodes); max(fmdl.nodes)];
brainNodes  = fmdl.nodes( fmdl.elems( fmdl.mat_idx{targetMat},: ), : );
brainC      = mean(brainNodes); brainC(3) = 0.03; % adjust brain z center to be within detection of all electrodes.
brainD      = mean( max(brainNodes) - min(brainNodes) );
brainR      = brainD/2;
fociR       = [1/16, 1/8, 1/4] * brainR;
fociDist    = [1/8,  1/4, 1/2] * brainR;
levels      = inf(32, 3);
levels(:,3) = imdl.rec_model.mdl_slice_mapper.z_pts;
radPos      = linspace(-pi, pi, 9);
radPos      = radPos(1:8); % 8 positions
res         = zeros(length(fociDist), length(radPos), length(fociR), 2);    % store loc and amp acc results
mRes        = zeros(length(fociDist), length(radPos), length(fociR), 544);  % store measurements
mResN       = zeros(length(fociDist), length(radPos), length(fociR), 544);  % store noisy measurements

for i=1:length(fociDist)
    dist = [1 1 0] * fociDist(i);
    for j=1:length(radPos) % correctly plots in circle.
        vx      = cos(radPos(j)); 
        vy      = sin(radPos(j));
        fociCtr = brainC + dist .* [vx vy 1];
        for k=1:length(fociR)
            img2            = img;
            
            % set delta sigma for foci elements
            eq              = sprintf('(x-(%f)).^2 + (y-(%f)).^2 + (z-(%f)).^2 < (%f)^2', fociCtr(1),fociCtr(2),fociCtr(3),fociR(k));
            select_fcn      = inline(eq, 'x','y','z');
            memb_frac       = find(elem_select( img2.fwd_model, select_fcn));
            img2.elem_data(memb_frac) = img2.elem_data(memb_frac) * DELTA;
            
            % simulate noisy measurement data
            vi              = fwd_solve(img2);
            mRes(i,j,k,:)   = vi.meas;
            viN = rmfield(vi, 'volt'); viN.meas = zeros(size(vi.meas));
            for m=1:TIMES
                viT      = add_noise( SNR, vi);
                viN.meas = viN.meas + viT.meas;
            end % end for m
            viN.meas        = viN.meas / TIMES;
            imgr = inv_solve(imdl, vh, viN);
            mResN(i,j,k,:)  = viN.meas;
            figure(2);clf();
            % show slices
            slcs3d      = squeeze(calc_slices(imgr, levels));
            r_img       = mk_mosaic( slcs3d, 0, [], IMGCOLS);
            c_img       = calc_colours( r_img, imgr);
            out_img     = reshape(c_img, size(r_img,1), size(r_img,2) ,[]);
            image(out_img); cb=colorbar; cb.Label.String='Relative Conductivity'; cb.Label.FontSize=18; axis equal; axis off; axis image;
            % adjust colorbar position to be same heigth as figure
            stp(out_img, mdlLims, fociCtr, fociR(k)); % show true position of foci on reconstruction
            ttl = sprintf('%s%.2f delta %.1f pctFociR %.1f pctDist %.2f radPos n=%.0f.pdf', SAVEDIR, DELTA, fociR(k)/brainR*100, (dist(i)/brainR*100), radPos(j),TIMES);
            % calc image stats
            if DELTA < 1
                posErr = img_quality(-slcs3d, brainR, fociCtr, mdlLims);
            else
                posErr = img_quality(slcs3d, brainR, fociCtr, mdlLims);
            end
            res(i,j,k,:)    = posErr;
            printPDF(ttl);
        end % end for k
    end % end for j
end % end for i

mResN2 = reshape(mResN, 544, 72);
res2 = reshape(res, 72, 2);

table_out = array2table(res2,'VariableNames',{'xy','z'});
writetable(table_out, sprintf('%s/%s%.2fDelta.csv',SAVEDIR, 'locErr',DELTA));

header = {};
for i=1:72
    header{i} = sprintf('meas%.0f',i);
end
table_out = array2table(mResN2,'VariableNames', header);
writetable(table_out, sprintf('%s/%s%.2fDelta.csv',SAVEDIR, 'Measurements',DELTA));

end % end function
% ======================================================================= %
function stp(img, mdlLims, fociCtr, fociR)
drawCtr = [0 0 0];
extent  = [0 0 0];
COLS = size(img,2)/32;
mdlDims = [26,32,32];
for i=1:3
    drawCtr(i)  = (fociCtr(i) - mdlLims(1,i)) / (mdlLims(2,i)- mdlLims(1,i)) * mdlDims(i); % which slice to draw foci with full radius
    extent(i)   = fociR / (mdlLims(2,i) - mdlLims(1,i)) * mdlDims(i); % number of pixels equal to sphere radius
end % end for
extent = ceil(mean(extent));
drawCtr = drawCtr + [3 0 0];
frameLL = LLIdx(round(drawCtr(3)), COLS);
viscircles([frameLL(2)+drawCtr(1), frameLL(1)-drawCtr(2)], extent, 'Color', 'k', 'LineWidth',1);
for ext=1:extent
    % draw lower extent
    frameLL = LLIdx(round(drawCtr(3))-ext, COLS);
    viscircles([frameLL(2)+drawCtr(1), frameLL(1)-drawCtr(2)], extent, 'Color', 'k', 'LineWidth',1);
    % draw upper extent
    frameLL = LLIdx(round(drawCtr(3))+ext, COLS);
    viscircles([frameLL(2)+drawCtr(1), frameLL(1)-drawCtr(2)], extent, 'Color', 'k', 'LineWidth',1);
end % end for
end % end function
% ======================================================================= %
function locErr = img_quality(slcs3d, brainR, fociCtr, mdlLims)

szSlc = size(slcs3d);
cog = zeros(3, 2); % xy, xz, yz
% x-y COG
count = [0 0 0];
for i = 1: szSlc(3) % each z slice
    a = COG(slcs3d(:,:,i)); % xy
    b = COG(slcs3d(:,i,:)); % xz
    c = COG(slcs3d(i,:,:)); % yz
    if sum(a) ~= -2
        cog(1,:) = cog(1,:) + a;
        count(1) = count(1) + 1;
    end
    if sum(b) ~= -2
        cog(2,:) = cog(2,:) + b;
        count(2) = count(2) + 1;
    end
    if sum(c) ~= -2
        cog(3,:) = cog(3,:) + c;
        count(3) = count(3) + 1;
    end
end
cog = cog ./ count';
x = (cog(1,1) + cog(2,1))/2;
y = (cog(1,2) + cog(3,1))/2;
z = (cog(2,2) + cog(3,2))/2;
cog = [x y z];

% draw cog
frameLL = LLIdx(round(cog(3)), 8);
drawEx(frameLL(2)+cog(1), frameLL(1)-cog(2), 3);

% loc error
trueCtr = [0 0 0];
mdlDims = [26,32,32];
brainR = brainR / (mdlLims(2,2)- mdlLims(1,2)) * 32;
for i=1:3
    trueCtr(i)  = (fociCtr(i) - mdlLims(1,i)) / (mdlLims(2,i)- mdlLims(1,i)) * mdlDims(i); % which slice to draw foci with full radius
end % end for
trueCtr(1) = trueCtr(1) + 3;
locErr          = sqrt( sum( (cog(1:2)-trueCtr(1:2)).^2, 2) ) / brainR; % xy component
locErr(2)       = abs(cog(3) - trueCtr(3)) / brainR;
end % end function
% ======================================================================= %
function cog = COG(slc)
cog = [-1 -1];
slc = squeeze(slc);
nPerRow = sum(~isnan(slc), 2); nPerRow(nPerRow==0)=1;
nPerCol = sum(~isnan(slc), 1); nPerCol(nPerCol==0)=1;
try
    slc(isnan(slc)) = min(slc(~isnan(slc)));
catch
    return % these -1 -1 cogs are excluded in the image quality function
end

slc = (slc - min(slc(:))) ./ (max(slc(:)) - min(slc(:))); % normalize
% slc = (slc./(sum(slc,'all'))) *100;

% COG X
col_sums = sum(slc, 1) ./ nPerCol; % represent as average per-pixel value
col_sums = (col_sums - min(col_sums(:))) ./ (max(col_sums(:)) - min(col_sums(:))); % normalize
col_sums = (col_sums./(sum(col_sums,'all'))) *100; % scale from 0 - 100
c = cumsum(col_sums);
cn = find(c <=50,1,'last');
if isempty(cn)
    cog(:) = 1;
end % end if
k = (50 - sum(col_sums(1: cn))) / cn;
cog(1) = round( ((cn+k+0.5) / (32 + 1)) * 32); % X coor

% COG Y
row_sums = sum(slc, 2) ./ nPerRow; % represent as average per-pixel value
row_sums = (row_sums - min(row_sums(:))) ./ (max(row_sums(:)) - min(row_sums(:))); % normalize
row_sums = (row_sums./(sum(row_sums,'all'))) *100; % scale from 0 - 100
r = cumsum(row_sums);
rn = find(r <= 50, 1,'last');
k = (50 - sum(row_sums(1: rn))) / rn;
cog(2) = round( ((rn + k + 0.5) / (32 + 1)) *32); % Y coor
if isnan(cog(1)) || isnan(cog(2))
    keyboard;
end % end if

end % end function
% ======================================================================= %
% num = [0 0 0];
%     den = 0;
%     for idx = 1:size(recon.Faces, 1)
%         x = mean(recon.XData(:,idx));
%         y = mean(recon.YData(:,idx));
%         z = mean(recon.ZData(:,idx));
%         amp = mean(recon.CData(:, idx)); % delta sigma_i
%         num = num + ([x y z] * amp);
%         den = den + amp;
%     end % end for
%     cog     = num/den;
%     meanAmp = 0;
%     N       = 0;
%     for idx = 1:size(recon.Faces, 1)
%         fc  = recon.Faces(idx, :); % gives vertice index
%         vtx = recon.Vertices(fc, :); % vertexes of face
%         pos = mean(vtx); % mean position of vertexes in face
%         che = sqrt( sum((cog-pos).^2, 2) ); % check if distance of pos is within radius
%         if che <= fociR
%             meanAmp = meanAmp + mean(recon.CData(:, idx));
%             N = N + 1;
%         end % end if
%     end % end for
%     meanAmp = meanAmp / N;
% end % end function

function frameLL = LLIdx(frame, COLS)
frameLL = [0 0];
r = ceil(frame/COLS);
c = frame - (r-1)*COLS;
frameLL(1) = r*32;
frameLL(2) = (c-1)*32 + 1;
end % end function
% ======================================================================= %
function h = Circle(x,y,r)
hold on
th = 0:pi/50:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
h = plot(xunit, yunit);
hold off
end % end function
% ======================================================================= %
function [h1,h2] = drawEx(x,y,r)
hold on
incr = -r:0.1:r;
xunit = incr + x;
yunit = incr + y;
h1 = plot(xunit, repmat(y,length(xunit),1), 'm-','LineWidth',2);
h2 = plot(repmat(x,length(xunit),1), yunit, 'm-','LineWidth',2);
hold off
end