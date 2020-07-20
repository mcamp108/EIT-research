% This function was written by :
%                               Mark Campbell
%                               Carleton University
% model x y z ranges from 0 - 0.1.

[fmdl, imdl] = get_hum_head_mdl(); % load models
% constants
DELTA   = 1.1; % ratio of brain conductivity between vh and vi.
SNR     = 80; % measurement signal to noise ratio. ref for SNR: https://www.ncbi.nlm.nih.gov/books/NBK549564/
IMGCOLS = 8;
SAVEDIR = 'C:\Users\Mark\Documents\GraduateStudies\LAB\EIT-neuroimaging\Figures\simulations\';
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

% parameters for setting up brain foci size and position.
targetMat   = 4;
matRef      = 0.47;
brainNodes  = fmdl.nodes( fmdl.elems( fmdl.mat_idx{targetMat},: ), : );
brainC      = mean(brainNodes); 
brainD      = mean( max(brainNodes) - min(brainNodes) );
brainR      = brainD/2;
fociR       = [1/16, 1/8, 1/4] * brainR;
fociDist    = [1/8, 1/4, 1/2] * brainR;
levels      = inf(32, 3);
levels(:,3) = imdl.rec_model.mdl_slice_mapper.z_pts;
% levels = flipud(levels);
radPos      = linspace(-pi, pi, 9);
radPos      = radPos(1:8); % 8 positions
res         = zeros(length(fociDist), length(radPos), length(fociR), 4);    % store loc and amp acc results
mRes        = zeros(length(fociDist), length(radPos), length(fociR), 544);  % store measurements
mResN       = zeros(length(fociDist), length(radPos), length(fociR), 544);  % store noisy measurements
%%
for i=1:length(fociDist)
    dist = [1 1 0] * fociDist(i);
    for j=1:length(radPos) % correctly plots in circle.
        vx      = cos(radPos(j)); 
        vy      = sin(radPos(j));
        fociCtr = brainC + dist .* [vx vy 1];
        for k=1:length(fociR)
            img2        = img;
            % set delta sigma for foci elements
            eq          = sprintf('(x-(%f)).^2 + (y-(%f)).^2 + (z-(%f)).^2 < (%f)^2', fociCtr(1),fociCtr(2),fociCtr(3),fociR(k));
            select_fcn  = inline(eq, 'x','y','z');
            memb_frac   = find(elem_select( img2.fwd_model, select_fcn));
            img2.elem_data(memb_frac) = img2.elem_data(memb_frac) * DELTA;
            % simulate noisy measurement data
            vi = fwd_solve(img2); 
            mRes(i,j,k,:)   = vi.meas;
            vi              = add_noise( SNR, vi);
            mResN(i,j,k,:)  = vi.meas;
            imgr            = inv_solve(imdl, vh, vi);
            figure(1);clf();
            recon           = show_fem(imgr);
            figure(2);clf();
            % show slices
            r_img           = mk_mosaic( squeeze(calc_slices(imgr, levels)) , 0, [], IMGCOLS);
            c_img           = calc_colours( r_img, imgr);
            out_img         = reshape(c_img, size(r_img,1), size(r_img,2) ,[]);
            image(out_img);axis equal;
            ttl = sprintf('%s%.1f pctFociR %.1f pctDist %.2f radPos', SAVEDIR, fociR(i)/brainR*100, (dist(i)/brainR*100), radPos(j));
            printPDF(ttl);
            % calc image stats
            [posErr,meanAmp]= img_quality(recon, brainR, fociCtr, fociR(k), DELTA-1);
            res(i,j,k,:)    = [posErr,meanAmp];
        end % end for k
    end % end for j
end % end for i
%%

function [locErr,meanAmp] = img_quality(recon, brainR, fociCtr, fociR, amp)
    [cog, meanAmp]  = COG(recon, fociR);
    locErr          = sqrt( sum( (cog(1:2)-fociCtr(1:2)).^2, 2) ) / brainR; % xy component
    locErr(3)       = abs(cog(3) - fociCtr(3)) / brainR;
    meanAmp         = meanAmp / amp;
%     szErr = 
%     ampErr = 
end % end function

function [cog, meanAmp] = COG(recon, fociR)
    num = [0 0 0];
    den = 0;
    for idx = 1:size(recon.Faces, 1)
        x = mean(recon.XData(:,idx));
        y = mean(recon.YData(:,idx));
        z = mean(recon.ZData(:,idx));
        amp = mean(recon.CData(:, idx)); % delta sigma_i
        num = num + ([x y z] * amp);
        den = den + amp;
    end % end for
    cog     = num/den;
    meanAmp = 0;
    N       = 0;
    for idx = 1:size(recon.Faces, 1)
        fc  = recon.Faces(idx, :); % gives vertice index
        vtx = recon.Vertices(fc, :); % vertexes of face
        pos = mean(vtx); % mean position of vertexes in face
        che = sqrt( sum((cog-pos).^2, 2) ); % check if distance of pos is within radius
        if che <= fociR
            meanAmp = meanAmp + mean(recon.CData(:, idx));
            N = N + 1;
        end % end if
    end % end for
    meanAmp = meanAmp / N;
end % end function

function printPDF(filename)
    h = gcf;
    set(h, 'PaperUnits','centimeters');
    set(h, 'Units','centimeters');
    pos=get(h,'Position');
    set(h, 'PaperSize', [pos(3) pos(4)]);
    set(h, 'PaperPositionMode', 'manual');
    set(h, 'PaperPosition',[0 0 pos(3) pos(4)]);
    print('-dpdf',filename);
end % end function

