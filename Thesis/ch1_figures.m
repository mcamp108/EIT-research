saveDir = 'C:\Users\Mark\Documents\GraduateStudies\LAB\Thesis\files\figures\Chap_2\';

% constants
nElec = 8;
INCR = 20;
SNR = 80;
% make model
imdl = mk_common_model('d2d3c',nElec);
fmdl = imdl.fwd_model;
fmdl.stimulation= mk_stim_patterns(nElec,1,[0,2],[0,2],{},1);
e = size(fmdl.elems, 1);
bkgnd = 1;
inclusions = [1,5];

for iter=1:length(inclusions)
    inclusion = inclusions(iter);
    
    % Solve Homogeneous model
    img = mk_image(fmdl, bkgnd);
    C = [0.5,0];
    r = 0.3;
    
    % add inclusion - set delta sigma for foci elements
    EQ              = sprintf('(x-(%f)).^2+(y-(%f)).^2<(%f)^2', C(1), C(2), r);
    select_fcn      = inline(EQ,'x','y','z');
    memb_frac       = find(elem_select( img.fwd_model, select_fcn));
    img.elem_data(memb_frac) = img.elem_data(memb_frac) * inclusion;
    img.fwd_solve.get_all_meas = 1;
    img.fwd_solve.get_all_nodes = true;
    fmdl.solve       = @fwd_solve_1st_order;
    fmdl.system_mat  = @system_mat_1st_order;
    vh = fwd_solve(img);
% ======================================================================= %
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
% ======================================================================= %
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
% ======================================================================= %
    % show current density
    % reset variables
    if iter==1
        close all
        figure('units','normalized','outerposition',[0 0 0.5 1]);
        figure('units','normalized','outerposition',[0.5 0 0.5 1]);
        m   = inf; M    = -inf;
        m2  = inf; M2   = -inf;
        vhR = vh;
    else
        % =============================================================== %
        vvR = add_noise(SNR, vh);
        imgr = inv_solve(imdl, vhR, vvR);
        figure('units','normalized','outerposition',[0 0 1 1]);
        subplot(1,2,1); show_fem(img); axis equal; axis off; axis tight; box on;
        ax=gca; pos = ax.Position;
        title(sprintf('Ground Truth - FEM with Background of %.1f and Inclusion of %.1f', inclusions(1),inclusions(2)),'FontSize',20);
        
        subplot(1,2,2); show_slices(imgr); axis equal; axis off; axis tight; box on;
        title(sprintf('Reconstruction Image with an SNR of %.1f dB', SNR),'FontSize',20);
        fg = gcf; fg.Colormap(1,:) = [1 1 1];
        colorbar;
        ax=gca; ax.Position(3:4) = pos(3:4);
        printPDF( sprintf('%sRecon', saveDir) );
        % =============================================================== %
    end % end if
    
    figure(1);
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

        % change colour mapping
        currentColormap = colormap(jet(INCR)); %// Get the current colormap
        [~, ~, ind] = histcounts(mags, size(currentColormap, 1)); %// Now determine the color to make each arrow using a colormap
        cmap = uint8(ind2rgb(ind(:), currentColormap) * 255); %// Now map this to a colormap to get RGB
        cmap(:,:,4) = 255;
        cmap = permute(repmat(cmap, [1 3 1]), [2 1 3]);
        
        % set background colour
        f = gcf;
        ylim([-1.1,1.1]);
        xlim([-1.1,1.1]);
        f.Children(i).Color=[1,1,1]*0.75;
        ax=gca;
        
        % find colour limits
        if ax.CLim(1) < m
            m = ax.CLim(1);
        end
        if ax.CLim(2) > M
            M = ax.CLim(2);
        end    
        axis equal; axis off; axis tight; box on;
    end % end for
% ======================================================================= %
    % show isopotentia lines
    figure(2);
    colormap(jet(INCR));
    vfmdl = imgv.fwd_model;
    
    % make meshgrid
    npx     = vfmdl.mdl_slice_mapper.npx;
    npy     = vfmdl.mdl_slice_mapper.npy;
    xmin    = min(vfmdl.nodes(:,1));   xmax = max(vfmdl.nodes(:,1));
    xmean   = mean([xmin,xmax]);       xrange= xmax-xmin;
    ymin    = min(vfmdl.nodes(:,2));   ymax = max(vfmdl.nodes(:,2));
    ymean   = mean([ymin,ymax]);       yrange= ymax-ymin;
    range   = max([xrange, yrange]);
    xspace  = linspace( xmean - range*0.5, xmean + range*0.5, npx );
    yspace  = linspace( ymean + range*0.5, ymean - range*0.5, npy );
    [X,Y]   = meshgrid( xspace, yspace );
    imgv2   = imgv;
    
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
        ax = gca;
        if ax.CLim(1) < m2
            m2 = ax.CLim(1);
        end
        if ax.CLim(2) > M2
            M2 = ax.CLim(2);
        end
    end % end for i

end % end for iter

% set color limits
fg = figure(1);
for i=1:nElec*2
    ax = fg.Children(i);
    ax.CLim = [m M];
end
set_colourbars(fg, nElec, INCR, 1);
% printPDF( sprintf('%sQlines.pdf', saveDir) );

fg = figure(2);
for i=1:nElec*2
    ax = fg.Children(i);
    ax.CLim = [m2 M2];
end
set_colourbars(fg, nElec, INCR, 2);
% printPDF( sprintf('%sVlines.pdf', saveDir) );
    
% ======================================================================= %
% ======================================================================= %

function set_colourbars(fg, nElec, incr, key)
    cbSpace = 1.1;
    ofst = 1;
    for ee = [8,16]
        pnl1 = ee-nElec+ofst;
        pnl2 = pnl1 + nElec/2;
        cbW = 0.015;
        cbTop = fg.Children(pnl2).Position(2) + fg.Children(pnl2).Position(4);
        cbBot = fg.Children(pnl1).Position(2);
        cbRig = (fg.Children(pnl1).Position(1) + fg.Children(pnl1).Position(3)) + cbW * cbSpace; 
        cb = colorbar;
        cb.Position = [cbRig cbBot cbW cbTop-cbBot]; % from R from B W H
        for i=1:length(cb.TickLabels)
            cb.TickLabels{i} = num2str( round( str2double(cb.TickLabels{i}) / incr, 3)  ); % scale from 0-1
        end
        
        if key==1
            ylabel(cb,'Log10 Current Magnitude (NU)');
        elseif key==2
            ylabel(cb,'Log10 Voltage Magnitude (NU)');
        end
        ofst = ofst + 1;
    end % end for
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