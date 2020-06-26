saveDir = 'C:\Users\Mark\Documents\GraduateStudies\LAB\Thesis\ch1\figures_1\';
% make model
% clf();
nElec = 8;
imdl = mk_common_model('d2d3c',nElec);
fmdl = imdl.fwd_model;

fmdl.stimulation= mk_stim_patterns(nElec,1,[0,2],[0,2],{},1);
e = size(fmdl.elems, 1);
bkgnd = 1;
% inclusion = 1;
inclusion = 5;
% Solve Homogeneous model
img = mk_image(fmdl, bkgnd);
C = [0.5,0];
r = 0.3;

% add inclusion
for i = 1:size(img.fwd_model.elems, 1)
    nodes = img.fwd_model.nodes( img.fwd_model.elems(i,:) , :);
    dist = sqrt(sum( (nodes-C).^2, 2));
    if sum(dist <= r) == 3
        img.elem_data(i) = bkgnd * inclusion;
    end % end if
end
% end for
img.fwd_solve.get_all_meas = 1;
img.fwd_solve.get_all_nodes = true;
% fmdl.solve       = @fwd_solve_1st_order;
% fmdl.system_mat  = @system_mat_1st_order;
vh = fwd_solve(img);

% voltage image setup
res = 128;
msm = struct('x_pts', linspace(-1,1,res), ...
            'y_pts', linspace(-1,1,res) );
imgv = rmfield(img,'elem_data');
imgv.fwd_model.mdl_slice_mapper = msm;
imgv.show_slices.contour_levels = linspace(-1,1,res);
imgv.show_slices.contour_properties = {'Color',[0,0,0]};
imgv.show_slices.axes_msm = true;
imgv.node_data = vh.volt;
imgv.fwd_model.mdl_slice_mapper.npx = res;
imgv.fwd_model.mdl_slice_mapper.npy = res;

% current image setup
imgc = img;
imgc.fwd_model.mdl_slice_mapper.npx = res;
imgc.fwd_model.mdl_slice_mapper.npy = res;
% show_fem(fmdl, [0 1 0]);
% axis equal;axis off
% fg=gcf;
% ax=fg.Children;
% for i=1:length(ax.Children)
%     ch = ax.Children(i);
%     if strcmp(class(ch), 'matlab.graphics.primitive.Text')
%         ch.FontSize = 20;
%     end
% end
% printPDF( sprintf('%sFEM', saveDir) );

% show current density
cbSpace = 1.1;

figure(1);
% clf();
m=inf;
M=-inf;
for i=1:nElec
    if inclusion == 1
        subplot(4,4,i);
    else
        subplot(4,4,i+nElec);
    end
    hold on;
    % show electrode numbering
    hh = show_fem(fmdl, [0 1 0]);
    set(hh,'EdgeColor',[1,1,1]);
    c = show_current(imgc, vh.volt(:,i));
    
    % show current density
    mags = log10( sqrt(c.xc.^2 + c.yc.^2) );
    mags(mags==-Inf)= nan;
    hContour = contourf(c.xp,c.yp, mags);
    
    % show current lines
    ss = streamslice(c.xp,c.yp,c.xc,c.yc);
    set(ss, 'Color', 'k','LineWidth',1.2 );

    % change color mapping
    currentColormap = colormap('jet'); %// Get the current colormap
    [~, ~, ind] = histcounts(mags, size(currentColormap, 1)); %// Now determine the color to make each arrow using a colormap
    cmap = uint8(ind2rgb(ind(:), currentColormap) * 255); %// Now map this to a colormap to get RGB
    cmap(:,:,4) = 255;
    cmap = permute(repmat(cmap, [1 3 1]), [2 1 3]);
    
    f = gcf;
    ylim([-1.1,1.1]);
    xlim([-1.1,1.1]);
    f.Children(i).Color=[1,1,1]*0.75;
    ax=gca;
    if ax.CLim(1) < m
        m = ax.CLim(1);
    end
    if ax.CLim(2) > M
        M = ax.CLim(2);
    end    
    axis equal; axis off; axis tight; box on;
end % end for

% set color limits
fg = gcf;
for i=1:nElec
    ax = fg.Children(i);
    ax.CLim = [m M];
end

% set colorbar
cbW = 0.015;
cbTop = fg.Children(5).Position(2) + fg.Children(5).Position(4);
cbBot = fg.Children(1).Position(2);
cbRig = (fg.Children(1).Position(1) + fg.Children(1).Position(3)) + cbW * cbSpace; 
cb = colorbar;
for i=1:length(cb.TickLabels)
    cb.TickLabels{i} = num2str( round( str2double(cb.TickLabels{i}) * M/64, 3)  );
end
cb.Position = [cbRig cbBot cbW cbTop-cbBot]; % from R from B W H
ylabel(cb,'Log10 Current Magnitude');
% printPDF( sprintf('%scurrentSim', saveDir) );

% equipotential
figure(2);

colormap('jet')
vfmdl = imgv.fwd_model;
% make meshgrid
npx  = vfmdl.mdl_slice_mapper.npx;
npy  = vfmdl.mdl_slice_mapper.npy;
xmin = min(vfmdl.nodes(:,1));    xmax = max(vfmdl.nodes(:,1));
xmean= mean([xmin,xmax]); xrange= xmax-xmin;
ymin = min(vfmdl.nodes(:,2));    ymax = max(vfmdl.nodes(:,2));
ymean= mean([ymin,ymax]); yrange= ymax-ymin;
range= max([xrange, yrange]);
xspace = linspace( xmean - range*0.5, xmean + range*0.5, npx );
yspace = linspace( ymean + range*0.5, ymean - range*0.5, npy );
[X,Y]=meshgrid( xspace, yspace );
imgv2 = imgv;
m=inf;
M=-inf;
for i = 1:nElec
    if inclusion == 1
        subplot(4,4,i);
    else
        subplot(4,4,i+nElec);
    end
    hold on;
    % show electrode numbering
    hh = show_fem(fmdl, [0 1 0]);
    set(hh,'EdgeColor',[1,1,1]);
    
    % show equipotential lines
    imgv2.node_data = vh.volt(:,i);
    V = log(abs(flipud(calc_slices( imgv2 ))));
    hContour = contourf(X,Y,V);
    axis equal; axis off; axis tight;
    ax=gca;
    if ax.CLim(1) < m
        m = ax.CLim(1);
    end
    if ax.CLim(2) > M
        M = ax.CLim(2);
    end
end % end for

% set color limits
fg = gcf;
for i=1:nElec
    ax = fg.Children(i);
    ax.CLim = [m M];
end

% set colorbar
cbW = 0.015;
cbTop = fg.Children(5).Position(2) + fg.Children(5).Position(4);
cbBot = fg.Children(1).Position(2);
cbRig = (fg.Children(1).Position(1) + fg.Children(1).Position(3)) + cbW * cbSpace; 
ax = gca;
cb = colorbar;
for i=1:length(cb.TickLabels)
    cb.TickLabels{i} = num2str( round( str2double(cb.TickLabels{i}) * M/64, 3)  );
end
cb.Position = [cbRig cbBot cbW cbTop-cbBot]; % from R from B W H
ylabel(cb,'Log10 Voltage Magnitude');
%%
    printPDF( sprintf('%sQlines', saveDir) );
%%
clf();
bkgnd = 1;
% inclusion = 1;
inclusion = 2;
% Solve Homogeneous model
img = mk_image(fmdl, bkgnd);
C = [0.5,0];
r = 0.3;
img2 = img;
% add inclusion
for i = 1:size(img2.fwd_model.elems, 1)
    nodes = img2.fwd_model.nodes( img2.fwd_model.elems(i,:) , :);
    dist = sqrt(sum( (nodes-C).^2, 2));
    if sum(dist <= r) == 3
        img2.elem_data(i) = bkgnd * inclusion;
    end % end if
end
% end for
img.fwd_solve.get_all_meas = 1;
img.fwd_solve.get_all_nodes = true;
% img2.fwd_solve.get_all_meas = 1;
% img2.fwd_solve.get_all_nodes = true;
vh = fwd_solve(img);
vv = fwd_solve(img2);
    
imgr = inv_solve(imdl, vh, vv);
figure(3);
subplot(1,2,2); ss = show_slices(imgr); colorbar; axis equal; axis off; axis tight; box on;
ss(ss==1)=0; f = gcf; f.Children(1).Color=[1,1,1]*0.75;
cb = colorbar; cb.Position(1) = 0.92;
subplot(1,2,1); show_fem(img2); axis equal; axis off; axis tight; box on;
figure(3); printPDF( sprintf('%sres', saveDir) );

%% Sensitivity map $Id: sensitivity_map01.m 4070 2013-05-26 21:22:03Z bgrychtol $
% 
% J = calc_jacobian(calc_jacobian_bkgnd(imdl));
% Sens = J(6,:)'./get_elem_volume(imdl.fwd_model);
% % Sens = J(5,:)'./get_elem_volume(imdl.fwd_model);
% img = mk_image(imdl, Sens);
% 
% img.calc_colours.npoints= 256;
% img.calc_slices.filter = conv2(ones(3),ones(3));
% img.calc_colours.clim = 0.5;
% show_slices(img);
% print_convert('sensitivity_map01a.png','-density 60');
% 
% img.calc_colours.cb_shrink_move = [0.5,0.7,0];
% clf;axis square,show_fem(img,[1,1]);axis off
% 
% 
% %%
% figure;hold on;
% for i=1:16
%     subplot(4,4,i);
%     q = show_current(imgc, vh.volt(:,i));
%     sy = linspace(-1,1 ,25); sx= sy;
%     hh = streamline(q.xp,q.yp, q.xc, q.yc, sx,sy); set(hh,'Linewidth',2);
%     hh = streamline(q.xp,q.yp,-q.xc,-q.yc,-sx,-sy); set(hh,'Linewidth',2);
% end
%%
% cd 'C:\Users\Mark\Pictures\PEI Photos';
% [I] = imread('PEI_2019 (19).jpg');  % Read the indexed image
% I = imresize(I, 0.25);
% J1 = imnoise(I,'gaussian',0,0.5);
% H = J1-I;
% subplot(1,3,1);image(I); axis equal; axis image; axis off;
% title('Input');
% subplot(1,3,2);image(H); axis equal; axis image; axis off;
% title('System');
% subplot(1,3,3);image(J1); axis equal; axis image; axis off;
% title('Output');
% printPDF( sprintf('%sinvprob1', saveDir) );
% 
% J2 = imnoise(I,'gaussian',0,0.5);
% subplot(1,3,1);image(J2); axis equal; axis image; axis off;
% title('Output');
% subplot(1,3,2);image(H); axis equal; axis image; axis off;
% title('Knowledge of System');
% subplot(1,3,3);image(J2-H); axis equal; axis image; axis off;
% title('Reconstructed Output');
% printPDF( sprintf('%sinvprob2', saveDir) );



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