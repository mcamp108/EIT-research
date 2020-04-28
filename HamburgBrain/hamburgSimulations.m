savedir = 'C:\Users\Mark\Documents\GraduateStudies\LAB\HamburgBrain\Figures\simulations\';

extra={'out','in', [ ...
       'solid cut = plane(0,0,1.1;0,0,1);' ...
       'solid in  = sphere(0,0,1;0.5) and cut;' ...
       'solid out = sphere(0,0,1;0.65) and (not in) and cut;' ...
       ]};
fmdl= ng_mk_cyl_models(2,[32,1],[0.05],extra);
% fmdl= ng_mk_cyl_models([2,1,0.1],[16,1],[0.05],extra); % fine mesh
stim = mk_stim_patterns(32,1,[0,5],[0,5],{'no_meas_current_next2'},1);
fmdl.stimulation = stim;
img = mk_image(fmdl,0.41); % scalp
img.elem_data(fmdl.mat_idx{2}) = 0.016; % skull
img.elem_data(fmdl.mat_idx{3}) = 0.47; % brain
img.calc_colours.clim = 1;
subplot(211); show_fem(img);
subplot(212); plot([v1.meas, v2.meas, 1e3*(v1.meas- v2.meas)]);
% ======================================================================= %
%% 0. Simulations on actual pig model
% pig= '8-2';
% [fmdl, ~] = get_pig_mdl(pig);
img = mk_image(fmdl, 1); % Background conductivity is scalp
imdl = get_imdl(img);
v1 = fwd_solve(img);

% img.elem_data([fmdl.mat_idx{1}]) = 1.1;    %   1: scalp        0.41
% img.elem_data([fmdl.mat_idx{2}]) = 0.9;   %   2: skull        0.016
img.elem_data([fmdl.mat_idx{3}]) = 1.1;    %   3: grey matter  0.47
% img.elem_data([fmdl.mat_idx{4}]) = 1.4;  %   4: air          0.0001
v2 = fwd_solve(img);

imgr = inv_solve(imdl,v1,v2);
show_slices(imgr);
    
    
%% 1. What is the minimum detectable whole-brain change in conducivity?
savedir = 'C:\Users\Mark\Documents\GraduateStudies\LAB\HamburgBrain\Figures\simulations\';
v1 = fwd_solve(img);
testCase = -10:10;
simMeas = zeros( size(v1.meas,1), length(testCase) );
legendEntries = cell(length(testCase),1);

for i=1:length(testCase)
    delta = testCase(i);
    img.elem_data(fmdl.mat_idx{2}) = 0.47 + 0.0047 * delta; % brain
    v2 = fwd_solve(img);
    simMeas(:, i) = 1e3*(v1.meas- v2.meas);
    legendEntries{i} = sprintf('%i%s', delta,'%');
end % end for

img.elem_data(fmdl.mat_idx{3}) = 0.47; % brain
plot(simMeas);
legend(legendEntries);
title('1000x measurement change with percentage change in whole brain conductivity');
% ======================================================================= %

%% 2. How does sensitivity change with change in position of skull and brain in body? 
% (DONE)
savedir = 'C:\Users\Mark\Documents\GraduateStudies\LAB\HamburgBrain\Figures\simulations\';
DELTAZ = 0.9;
stD = 1.8; % soft tissue diameter
brD = 0.55; % brain diameter
skD = 0.65; % skull diameter
THETA = stD * 0.01; % focus edge from brain edge safety factor
maxChng = stD/2 - (skD/2 + THETA);
increment = maxChng/10;
legendEntries = cell(length(0:10),1);
big_fig();
colormap jet;
cmap = colormap;
cIdx = round(linspace(1,size(colormap,1),11));
colours = cmap(cIdx,:);
nIter = length(0:10);
for i=0:10
%     if i==6;keyboard;end
    newC = -(increment * i);
    img = mk_sim_img(stD, [0, newC, skD], [0, newC, brD]); % st, sk, br
    imdl = get_imdl(img);
    v1 = fwd_solve(img);
    legendEntries{i+1} = sprintf('%0.f%s', abs(newC/maxChng)*100,'%');
    
    % bottom row shows reconstructions
    subplot(5, nIter, 4*nIter+i+1);
    fem_top_view(img);
    xlabel(legendEntries{i+1});
    
    % show reconstructions
    img.elem_data(img.fwd_model.mat_idx{2}) = 0.47 * DELTAZ; % brain
    v2 = fwd_solve(img);
    imgr = inv_solve(imdl,v1,v2);
    if i==0
        elemData = zeros(size(imgr.elem_data,1),nIter);
    end
    elemData(:,i+1) = imgr.elem_data;
    
    % show measurements   
    subplot(5, nIter, 1:nIter*3);
    hold on; plot(1e3*(v1.meas- v2.meas), 'color',colours(i+1,:));
    
end % end for

imgr.clim = max(elemData(:));

for i=1:nIter
    imgr.elem_data = elemData(:,i);
    subplot(5, nIter, 3*nIter+i);
    show_slices(imgr);
    view(180,90);
end % end for

subplot(5, nIter, 1:nIter*3);
legend(legendEntries);
title('Measurement with whole brain conductivity decrease of 10% and change of brain/skull center (% delta radii st-sk)');
xlabel('Measurement number');
ylabel('Voltage change from reference (1000x)');
xlim([0 size(v1.meas,1)]);
saveas( gcf, horzcat(savedir,'skullBrainLocChangeZDownsameclimFine.svg') );
% ======================================================================= %

%% 3. How does sensitivity change wih thickness of skull as percentage of brain radius?
savedir = 'C:\Users\Mark\Documents\GraduateStudies\LAB\HamburgBrain\Figures\simulations\';
%prep
brainR = 0.5;
maxChng = brainR;
increment = maxChng / 20;
legendEntries = cell(length(0:10),1);
stim = mk_stim_patterns(32,1,[0,5],[0,5],{'no_meas_current_next2'},1);
% colours
figure();
colormap jet;
cmap = colormap;
cIdx = round(linspace(1,size(colormap,1),11));
colours = cmap(cIdx,:);

for i=0:10
    newD = brainR + (increment * i);
    extra ={'out','in', [ ...
       'solid cut = plane(0,0,1.1;0,0,1);' ...
       'solid in  = sphere(0,0,1;0.5) and cut;' ...
       sprintf('solid out = sphere(0,0,1;%f) and (not in) and cut;', newD) ...
       ]};
    fmdl = ng_mk_cyl_models(2, [32,1], [0.05], extra);
    fmdl.stimulation = stim;
    img = mk_image(fmdl,0.41); % scalp
    if i==0
        img.elem_data(fmdl.mat_idx{2}) = 0.47; % brain
    else
        img.elem_data(fmdl.mat_idx{2}) = 0.016; % skull
        img.elem_data(fmdl.mat_idx{3}) = 0.47 * 1.5; % brain
    end % end if
    img.calc_colours.clim = 1;
    v1 = fwd_solve(img);
    legendEntries{i+1} = sprintf('%0.f%s', abs((increment * i)/maxChng)*100,'%');
        
    subplot(4, 11, 34+i);
    fem_top_view(img);
    xlabel(legendEntries{i+1});
    
    if i==0
        img.elem_data(fmdl.mat_idx{2}) = 0.47 * 0.9; % brain
    else
        img.elem_data(fmdl.mat_idx{3}) = 0.47 * 0.9; % brain
    end % end if
    v2 = fwd_solve(img);    
    subplot(4, 11, 1:33);
    hold on; plot(1e3*(v1.meas- v2.meas), 'color',colours(i+1,:), 'LineWidth',1);
    
end % end for
legend(legendEntries);
title('Measurements with whole brain conductivity decrease of 10% as a function of skull thickness (% of brain radius)');
xlabel('Measurement number');
ylabel('Voltage change from reference (1000x)');
xlim([0 size(v1.meas,1)]);
savedir = 'C:\Users\Mark\Documents\GraduateStudies\LAB\HamburgBrain\Figures\simulations\';
saveas( gcf, horzcat(savedir,'measChangewSkullThicknessChange.svg') );
% ======================================================================= %

%% 4. How does sensitivity change with focus size?
savedir = 'C:\Users\Mark\Documents\GraduateStudies\LAB\HamburgBrain\Figures\simulations\';
stim = mk_stim_patterns(32,1,[0,5],[0,5],{'no_meas_current_next2'},1);
FCSZ = 0.5; % focus size as a proportion of brain size (0.4 smalles i could get so far with coarse mesh)
stD = 1.8; % soft tissue diameter
brD = 0.55; % brain diameter
skD = 0.65; % skull diameter
fcD = brD * FCSZ; % focus diameter
fcDelta = 0.1; % focus conductivity as a proportion of normal brain conductivity

THETA = brD * 0.01; % focus edge from brain edge safety factor
lim = brD/2 - (fcD/2 + THETA);
brC = [0, -stD/4]; % brain center
skC = brC; % skull center
fcC = [ brC(1), brC(2);
        brC(1), brC(2); 
        brC(1), brC(2)+lim; 
        brC(1)+lim, brC(2);
        brC(1), brC(2)-lim;
        brC(1)-lim, brC(2) ];
legendEntries = {'control','center','bottom','left','top','right'};
nIter = size(fcC, 1);
figure();
colormap jet;
cmap = colormap;
cIdx = round( linspace( 1,size(colormap,1), nIter ) );
colours = cmap(cIdx,:);

for i = 1: nIter
    extra ={'focus', 'brain', 'skull' [ ...
       'solid cut = plane(0,0,1.1;0,0,1);' ...
       sprintf('solid focus = sphere(%f,%f,1;%f) and cut;', fcC(i,1), fcC(i,2), fcD/2) ...
       sprintf('solid brain = sphere(%f,%f,1;%f) and (not focus) and cut;', brC(1),brC(2), brD/2) ...
       sprintf('solid skull = sphere(%f,%f,1;%f)and (not focus) and (not brain) and cut;', skC(1), skC(2), skD/2) ...
       ]};
    fmdl = ng_mk_cyl_models([2, stD/2], [32,1], [0.05], extra);
    fmdl.stimulation = stim;
    img = mk_image(fmdl,0.41); % scalp
    img.elem_data(fmdl.mat_idx{2}) = 0.016; % skull
    img.elem_data(fmdl.mat_idx{3}) = 0.47; % brain
    img.elem_data(fmdl.mat_idx{4}) = 0.47; % ischemic region
    img.calc_colours.clim = 1;
    v1 = fwd_solve(img);
    
    subplot(4, nIter, 3*nIter+i);
    if i==1
        fem_top_view(img,'control');
        img.elem_data(fmdl.mat_idx{4}) = 0.47;
    else
        fem_top_view(img);
        img.elem_data(fmdl.mat_idx{4}) = 0.47 * fcDelta;
    end
    
    xlabel(legendEntries{i});
    v2 = fwd_solve(img);    
    subplot(4, nIter, 1: 3*nIter);
    hold on; plot(1e3*(v1.meas- v2.meas), 'color',colours(i,:));
    
end % end for
legend(legendEntries);
title('Measurements with brain focus at 10% normal brain conductivity');
xlabel('Measurement number');
ylabel('Voltage change from reference (1000x)');
xlim([0 size(v1.meas,1)]);
% ======================================================================= %

%% 5. How does sensitivity change with focus location?
savedir = 'C:\Users\Mark\Documents\GraduateStudies\LAB\HamburgBrain\Figures\simulations\5';
fcDelta = 0.1; % focus conductivity as a proportion of normal brain conductivity
brDelta = 1.2;

FCSZ = 0.5; % focus size as a proportion of brain size (0.4 smalles i could get so far with coarse mesh)
stD = 1.8; % soft tissue diameter
brD = 0.55; % brain diameter
skD = 0.65; % skull diameter
fcD = brD * FCSZ; % focu diameter
THETA = brD * 0.01; % focus edge from brain edge safety factor
lim = brD/2 - (fcD/2 + THETA);
brC = [0, -stD/4]; % brain center
skC = brC; % skull center
fcC = [ brC(1), brC(2); 
        brC(1), brC(2); 
        brC(1), brC(2)+lim; 
        brC(1)+lim, brC(2);
        brC(1), brC(2)-lim;
        brC(1)-lim, brC(2) ];
legendEntries = {'control','center','bottom','left','top','right'};
nIter = size(fcC, 1);
big_fig();
colormap jet;
cmap = colormap;
cIdx = round( linspace( 1,size(colormap,1), nIter ) );
colours = cmap(cIdx,:);

for i = 1: nIter
    img = mk_sim_img(stD, [skC,skD], [brC, brD], [fcC(i,:),fcD], 0);
    imdl = get_imdl(img);
    v1 = fwd_solve(img);
    
    % show forward solution
    subplot(5, nIter, 4*nIter+i);
    fem_top_view(mk_sim_img(stD, [skC,skD], [brC, brD], [fcC(i,:),fcD], 1));
    xlabel(legendEntries{i});
    
    % show reconstructions
    if i == 1
        img.elem_data(img.fwd_model.mat_idx{2}) = 0.47 * brDelta; % changes same as brain
    else
        img.elem_data(img.fwd_model.mat_idx{2}) = 0.47 * fcDelta;
    end % end if
    
    img.elem_data(img.fwd_model.mat_idx{3}) = 0.47 * brDelta;
%     img.elem_data(img.fwd_model.mat_idx{1}) = 0.41 * brDelta;
    v2 = fwd_solve(img);
    imgr = inv_solve(imdl,v1,v2);
    
    if i==1
        elemData = zeros(size(imgr.elem_data,1),nIter);
    end
    elemData(:,i) = imgr.elem_data;
    
    % show measurements
    subplot(5, nIter, 1: 3*nIter);
    hold on; plot(1e3*(v1.meas- v2.meas), 'color',colours(i,:));
    
end % end for

imgr.clim = max(elemData(:));

for i=1:nIter
    imgr.elem_data = elemData(:,i);
    subplot(5, nIter, 3*nIter+i);
    show_slices(imgr);
    view(180,90);
end % end for

subplot(5, nIter, 1: 3*nIter);
legend(legendEntries);
title(sprintf('Measurements with brain focus at %0.f%s and brain at %0.f%s normal brain conductivity', fcDelta*100,'%',brDelta*100,'%'));
xlabel('Measurement number');
ylabel('Voltage change from reference (1000x)');
xlim([0 size(v1.meas,1)]);
saveas( gcf, horzcat(savedir, sprintf('fociLocZ%0.f%sGlobZandSoftTissueZ%0.f%sFine.svg', fcDelta*100,'%',brDelta*100,'%' )));
% ======================================================================= %

%% 6. How does sensitivity change with focus intensity?
savedir = 'C:\Users\Mark\Documents\GraduateStudies\LAB\HamburgBrain\Figures\simulations\6';
big_fig();
loc = {'center','bottom','left','top','right'};
do = 3; % which location to show.
doLoc = loc{do};

FCSZ = 0.5;         % focus size as a proportion of brain size (0.4 smalles i could get so far with coarse mesh)
stD = 1.8;          % soft tissue diameter
brD = 0.55;         % brain diameter
skD = 0.65;         % skull diameter
fcD = brD * FCSZ;   % focus diameter
brC = [0, -stD/4];  % brain center
skC = brC;          % skull center

% brD = 1.4;          % brain diameter human
% skD = 1.6;          % skull diameter human
% fcD = brD * FCSZ;   % focus diameter human
% brC = [0, 0];       % brain center human
% skC = brC;          % skull center human

THETA = brD * 0.05; % focus edge from brain edge safety factor
lim = brD/2 - (fcD/2 + THETA);

fcC = [ brC(1), brC(2); 
        brC(1), brC(2)+lim; 
        brC(1)+lim, brC(2);
        brC(1), brC(2)-lim;
        brC(1)-lim, brC(2) ];
legendEntries = {'center','bottom','left','top','right'};
doRange = fliplr(-0.5 : 0.1 : 0.5);
fcRange = fliplr(-0.9 : 0.1 : 0);
nIter = length(doRange);
fcIter = length(fcRange);
opt.cut=1;
img = mk_sim_img(stD, [skC,skD], [brC, brD], [fcC(do,:),fcD], opt);
imdl = get_imdl(img);
v1 = fwd_solve(img);
count=1;

for i =1:fcIter
    img.elem_data(img.fwd_model.mat_idx{2}) =       0.47 + (0.47 * fcRange(i)); % focus   
    for j =1:nIter
        img.elem_data(img.fwd_model.mat_idx{3}) =   0.47 + (0.47 * doRange(j)); % brain
        img.elem_data(img.fwd_model.mat_idx{1}) =   0.41 + (0.41 * doRange(j)); % soft tissue
        v2 = fwd_solve(img);
        imgr = inv_solve(imdl, v1, v2);
        
        if count == 1
            elemData = zeros(size(imgr.elem_data,1), nIter * fcIter);
        end % end if
        if j==1
            bsln = imgr.elem_data; % baseline simulating second reference at average diastole.
        end % end if
        elemData(:,count) = imgr.elem_data - bsln;
        elemData(:,count) = imgr.elem_data;
        count = count + 1;
    end % end for j
end % end for i

imgr.elem_data = elemData;
imgr.show_slices.img_cols = nIter;
imgr.calc_slices.clim = max(imgr.elem_data(:));
show_slices(imgr);
title(sprintf('Simulated pulsatile image set with %s side ischemia', doLoc));

% format axes
ax = gcf;
ax.Visible = 'on';
axC = ax.Children;
axC.Visible = 'on';
% Y labels
seqYTicksRef = linspace(axC.YLim(1), axC.YLim(2), fcIter*2+1);
axC.YTick = seqYTicksRef(2:2:length(seqYTicksRef));

for j = 1:fcIter
    zVal = (0.47 + (0.47 * fcRange(j)) ) / 0.47 * 100;
    axC.YTickLabel{j} = sprintf('%0.f%s', zVal, '%');
end % end for

% X labels
seqXTicksRef = linspace(axC.XLim(1), axC.XLim(2), nIter*2+1);
axC.XTick = seqXTicksRef(2:2:length(seqXTicksRef));

for j = 1:nIter
    zVal = (0.47 + (0.47 * doRange(j)) ) / 0.47 * 100;
    axC.XTickLabel{j} = sprintf('%0.f%s', zVal, '%');
end % end for

axC.TickLength = [0,0];
axC.FontSize = 15;
ylabel('Focus Conductivity (% Normal Brain)');
xlabel('Brain Conductivity (% Normal Brain)');
view(180,90)
% saveas( gcf, horzcat(savedir, sprintf('SubDiastolePulsatilitySimulation %0.1fFocSize %swithSofttisuedeltaZFine.svg', FCSZ, doLoc)) );
saveas( gcf, horzcat(savedir, sprintf('SubDiastolePulsatilitySimulation %0.1fFocSize %sFine.svg', FCSZ, doLoc)) );

% ======================================================================= %
% ======================================================================= %

function fem_top_view(img, sel)

control = false;

if nargin == 2 
    if strcmp(sel, 'control')
        control = true;
    end % end if
end % end if 

switch length(img.fwd_model.mat_idx)
    case 2 
        img.elem_data(img.fwd_model.mat_idx{2}) = 0.47 * 1.5; % brain
    case 3 
        img.elem_data(img.fwd_model.mat_idx{2}) = 0.47 * 1.5; % brain
    case 4
        img.elem_data(img.fwd_model.mat_idx{2}) = 0; % focus
        img.elem_data(img.fwd_model.mat_idx{3}) = 0.47 * 1.5; % brain
        if control
            img.elem_data(img.fwd_model.mat_idx{2}) = 0.47 * 1.5; % focus
            img.elem_data(img.fwd_model.mat_idx{3}) = 0.47 * 1.5; % brain
        else
            
        end % end if
end % end switch
show_fem(img);
view(180,90);
set(gca,'xtick',[]);
set(gca,'xticklabel',[]);
set(gca,'ytick',[]);
set(gca,'yticklabel',[]);
end % end function

% ======================================================================= %

function imdl = get_imdl(img)
opt.imgsz= [64, 64];
imdl = mk_GREIT_model(img, 0.2, 1, opt);
end

% ======================================================================= %

function img = mk_sim_img(stD, sk, br, fc, opt)
% soft tissue, skull, brain focus arrays. first is center x and y, second
% is diameter. Enter only diameter for soft tissue.
if isfield(opt, 'cut')
    cut=opt.cut;
else
    cut=0;
end
if isfield(opt, 'skip')
    skip=opt.skip;
else
    skip=4;
end

if cut==0
    if nargin == 4
       extra ={'brain', 'skull', [ ...
       sprintf('solid brain = sphere(%f,%f,1;%f);', br(1),br(2), br(3)/2) ...
       sprintf('solid skull = sphere(%f,%f,1;%f) and (not brain);', sk(1), sk(2), sk(3)/2) ...
       ]};           
    elseif nargin == 5
       extra ={'focus', 'brain', 'skull', [ ...
       sprintf('solid focus = sphere(%f,%f,1;%f);', fc(1), fc(2), fc(3)/2) ...
       sprintf('solid brain = sphere(%f,%f,1;%f) and (not focus);', br(1),br(2), br(3)/2) ...
       sprintf('solid skull = sphere(%f,%f,1;%f) and (not focus) and (not brain) ;', sk(1), sk(2), sk(3)/2) ...
       ]}; 
    end % end if
else % have cut
    if nargin == 4
       extra ={'brain', 'skull', [ ...
       'solid cut = plane(0,0,1.1;0,0,1);' ...
       sprintf('solid brain = sphere(%f,%f,1;%f) and cut;', br(1),br(2), br(3)/2) ...
       sprintf('solid skull = sphere(%f,%f,1;%f) and (not brain) and cut;', sk(1), sk(2), sk(3)/2) ...
       ]};           
    elseif nargin == 5
       extra ={'focus', 'brain', 'skull', [ ...
       'solid cut = plane(0,0,1.1;0,0,1);' ...
       sprintf('solid focus = sphere(%f,%f,1;%f) and cut;', fc(1), fc(2), fc(3)/2) ...
       sprintf('solid brain = sphere(%f,%f,1;%f) and (not focus) and cut;', br(1),br(2), br(3)/2) ...
       sprintf('solid skull = sphere(%f,%f,1;%f) and (not focus) and (not brain) and cut;', sk(1), sk(2), sk(3)/2) ...
       ]}; 
    end % end if
end

stim = mk_stim_patterns(32,1,[0,skip+1],[0,skip+1],{'no_meas_current_next2'},1);
fmdl= ng_mk_cyl_models([2, stD/2, 0.1], [32,1], [0.05], extra); % fine
% fmdl = ng_mk_cyl_models([2, stD/2], [32,1], [0.05], extra); % coarse
fmdl.stimulation = stim;
img = mk_image(fmdl,0.41); % scalp

if nargin == 4
    img.elem_data(fmdl.mat_idx{2}) = 0.47; % brain
    img.elem_data(fmdl.mat_idx{3}) = 0.016; % skull    
elseif nargin == 5
    img.elem_data(fmdl.mat_idx{2}) = 0.47; % focus
    img.elem_data(fmdl.mat_idx{3}) = 0.47; % brain
    img.elem_data(fmdl.mat_idx{4}) = 0.016; % skull
end

img.calc_colours.clim = 1;

end % end function

% ======================================================================= %

function big_fig()
    figure('units','normalized','outerposition',[0 0 1 1]);
end % end function

% ======================================================================= %

function wireframe(fmdl, material)
nodeIdx = fmdl.elems( fmdl.mat_idx{material},: );
hold on;
h = trimesh(nodeIdx, fmdl.nodes(:,1),fmdl.nodes(:,2),fmdl.nodes(:,3));
h.FaceAlpha = 0;
hold off;
end % end function

function scat3(fmdl, mat)
nodes = fmdl.nodes( fmdl.elems( fmdl.mat_idx{mat},: ), : );
hold on;
scatter3(nodes(:,1), nodes(:,2), nodes(:,3));
hold off;
end % end function