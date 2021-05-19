% load data
POS1    = 'P07,P08,P10,P11,P12'; % 5
POS2    = 'P04,P05,P06,P09,P15,P18'; % 6
POS3    = 'P01,P02,P03,P13,P17,P19'; % 6
bsln = 2;
file = 'E:\University\Masters\EIT-restraint\zzMC\features\per_observation_features_selfRef.csv';
data = readtable(file);

phases = {'U','R','W1','W2','W3','W4','W5','X1','X2','X3','X4','X5','P'};

%%
close all
dat = data.FRCMRELind;
nSubjs = 19;
colours = linspace(0, 1, nSubjs)';
colours = [colours, flipud(colours), colours];

xshift = linspace(-0.25, 0.25, 19);
figure();
for p = 1:19
    
%     figure();
    if p < 10
        subj = sprintf('P0%.0f', p);
    else
        subj = sprintf('P%.0f', p);
    end
    idx1 = strcmp(data.subj, subj);
    for i = 1:length(phases)
        phase = phases{i};
        idx2 = strcmp(data.phase, phase);

        selection = (idx1 + idx2) == 2;
        idx = find(selection);
        xax = repmat(i + xshift(p), length(idx), 1);
        plot(xax, dat(idx), 'o', 'MarkerFaceColor', colours(p, :), 'MarkerEdgeColor', 'k');
        hold on;
    end
    title(subj);
end
%%
% rename P to P1 and P2
for row =1:size(data,1)
    thisRow = data(row,:);
    if strcmp(thisRow.phase, 'P')
        if thisRow.time <= 60
            data.phase{row} = 'P1';
        else
            data.phase{row} = 'P2';
        end
    end
    
end
writetable(data, file);

%%
for i = 1:size(data, 1)
    row = data(i, :);
    if strcmp(row(1), 'P01')
        
    end
end
    

statdir = 'C:\Users\Mark\Documents\GraduateStudies\LAB\EIT-restraint\zzMC\data\phys_meas\features';
[D, header] = load_files(statdir);
params = {'BF','BFREL','TV','TVREL','FRCB','FRCBREL','FRCM','FRCMREL','FRCT','FRCTREL','MINVNT','MINVNTREL'};
E = get_position_data(D, header, params);

saveDir = 'C:\Users\Mark\Documents\GraduateStudies\LAB\EIT-restraint\zzMC\data\phys_meas\features\Roslyn';
fn = fieldnames(E);
for i = 1:length(fn)
    tbl = participant_data_to_table(E.(fn{i}));
    save_name = sprintf('%s%s%s_features.csv',saveDir,'\',fn{i});
    writetable(tbl, save_name);
end


% results_figs(E);
%%
% tests    = {'U:R';'W1:R';'W5:R';'W1:W5';'X1:X5';'W1:X1';'W5:X5'};
% tests    = {'U:R'};
% F = struct;
% F.pos1=E.pos1;
% F.pos2=E.pos2;
% F.pos3=E.pos3;
% F.all=E.all;
% outtable(F, tests, statdir);
% %%
% age = [22,41,22,20,38,28,27,21,21,21,19,21,21,21,44,21,42];
% bmi = [22.4,34.0,22.9,29.3,28.0,29.4,31.5,23.9,21.0,25.9,28.5,24.0,25.8,22.6,36.9,31.3,21.6];
% E = get_position_data(D, header, params);
% U = 1;
% R = 2;
% W = 3:7;
% X = 8:12;
% P = 13;
% params  = {'MINVNTREL','FRCMREL'};
% param   = params{2};
% posCond = X(5);
% allCond = W(5);
% 
% [pooledw5, pos] = pos_vs_all(E, param, posCond, X(1));
% y1 = pooledw5(:,4);
% [pooledx5, ~] = pos_vs_all(E, param, posCond, W(5));
% y2 = pooledx5(:,4);
% y = y2 - y1;
% % groups = {pos,age',bmi'};
% % test if posture at X5 was different than group mean at W5.
% % [pooled, pos] = pos_vs_all(E, param, posCond, allCond);
% % p = anovan(y,{pos}); fprintf('\npos1 vs all for %s p = %.5f',param,p);
% % p = anovan(y,{age},'continuous',1); fprintf('\npos1 vs all for %s p = %.5f',param,p);
% % p = anovan(y,{bmi},'continuous',1); fprintf('\npos1 vs all for %s p = %.5f',param,p);
% % p = anovan(y1,groups); fprintf('\npos1 vs all for %s p = %.5f',param,p);
% % p = anovan(y2,groups); fprintf('\npos1 vs all for %s p = %.5f',param,p);
% % p = anova1(pooled); fprintf('\npos1 vs all for %s p = %.5f',param,p);
% y1 = y1(~isnan(y2));
% y2 = y2(~isnan(y2));
% p = anova1([y1,y2]); fprintf('\npos1 vs all for %s p = %.5f',param,p);
% 
% 
% % p = anova1(pooled(:,[1,4])); fprintf('\npos1 vs all for %s p = %.5f',param,p);
% % p = anova1(pooled(:,[2,4])); fprintf('\npos2 vs all for %s p = %.5f',param,p);
% % p = anova1(pooled(:,[3,4])); fprintf('\npos3 vs all for %s p = %.5f',param,p);
% % % test if posture at X5 was different than control.
% % [p,anovatab,stats] = anova1(pooled(:,[2,1])); fprintf('\npos2 vs pos1 for %s p = %.5f',param,p);
% % [p,anovatab,stats] = anova1(pooled(:,[3,1])); fprintf('\npos3 vs pos1 for %s p = %.5f',param,p);
% 
% %%
% % 6.	X1 for all postures
% % 7.	X5 for all postures – is there significance where previously there was none?
% % One-way anova shows that cannot reject the null that means for MINVNT during X1 are the same across all positions.
% % One-way anova shows that means for MINVNT during X5 are not the same for pos1 and pos3 but cannot reject for pos1 and pos2.
% % One-way anova shows that cannot reject the null that means for FRCM during X1 are the same across all positions.
% % One-way anova shows that cannot reject the null that means for FRCM during X1 are the same across all positions.
% % One-way anova shows that cannot reject the null that means for FRCM during X5 are the same across all positions.
% 
% %%
% % Params correlated with age or BMI?
% W = 3:7;
% X = 8:12;
% params  = {'MINVNTREL','FRCMREL'};
% posCond = X(5);
% allCond = W(5);
% param   = params{2};
% age = [22,41,22,20,38,28,27,21,21,21,19,21,21,21,21,22,44,21,42,23,25];
% bmi = [22.4,34.0,22.9,29.3,28.0,29.4,31.5,23.9,21.0,25.9,28.5,24.0,25.8,32.3,22.6,38.5,36.9,31.3,21.6,26.4,34.7];
% fn = fieldnames(E);
% for i=1:length(fn)
%    try
%        idx = str2double(fn{i}(2:end));
%        if ~isnan(idx) && idx~=14 && idx~=16
%            E.(fn{i}).age = age(idx);
%            E.(fn{i}).bmi = bmi(idx);
%        end
%    catch
%        continue
%    end
% end
% 
% %%
% params  = {'MINVNTREL','FRCMREL'};
% param   = params{2};
% xVars = {'age','bmi'};
% xVar = xVars{2};
% 
% x1 = []; x2 = [];
% y1 = []; y2 = [];
% for i=1:length(fn)
%    try
%        idx = str2double(fn{i}(2:end));
%        if ~isnan(idx) && idx~=14 && idx~=16
%            x1 = [x1, E.(fn{i}).age];
%            x2 = [x2, E.(fn{i}).bmi];
%            y1 = [y1, E.(fn{i}).(param).data(W(1))];
%            y2 = [y2, E.(fn{i}).(param).data(W(5))-E.(fn{i}).(param).data(W(1))];
%        end
%    catch
%        continue
%    end    
% end
% subplot(2,2,1);plot(x1,y1, 'o');xlabel(xVars{1});ylabel(param);title('W1');
% subplot(2,2,2);plot(x1,y2, 'o');xlabel(xVars{1});ylabel(param);title('W5-W1');
% subplot(2,2,3);plot(x2,y1, 'o');xlabel(xVars{2});ylabel(param);title('W1');
% subplot(2,2,4);plot(x2,y2, 'o');xlabel(xVars{2});ylabel(param);title('W5-W1');
% 
% fit1 = fitlm(x1, y1);
% fit2 = fitlm(x1, y2);
% fit3 = fitlm(x2, y1);
% fit4 = fitlm(x2, y2);
%%
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
function tbl = participant_data_to_table(Participant)
startingStruct = struct;
fn = fieldnames(Participant);
for i = 1:length(fn)
   variable = fn{i};
   variableFns = fieldnames(Participant.(variable));
   for j = 1:length(variableFns)
       if ~isempty(Participant.(variable).(variableFns{j}))
           startingStruct.(sprintf('%s_%s', variable, variableFns{j})) = Participant.(variable).(variableFns{j});
       end
   end
end
tbl = struct2table(startingStruct);
end % end function


function [pooled,pos] = pos_vs_all(E, param, posCond, allCond)
POS1    = 'P07,P08,P10,P11,P12'; % 5
POS2    = 'P04,P05,P06,P09,P15,P18'; % 6
POS3    = 'P01,P02,P03,P13,P17,P19'; % 6
fn      = fieldnames(E);
nPart   = length(split(sprintf('%s,%s,%s',POS1,POS2,POS3),','));
pooled  = nan(nPart,4);
pos     = nan(nPart,1);
pstn = 0;
a = 0;
b = 0;
p1idx=0; p2idx=0; p3idx=0; allIdx=0;
for i = 1:length(fn)
    do = true;
    p = fn{i};
    if contains(POS1, p)
        p1idx=p1idx+1;
        pooled(p1idx,1) = E.(p).(param).data(posCond);
        pstn = 1;
    elseif contains(POS2, p)
        p2idx=p2idx+1;
        pooled(p2idx,2) = E.(p).(param).data(posCond);
        pstn = 2;
    elseif contains(POS3, p)
        p3idx=p3idx+1;
        pooled(p3idx,3) = E.(p).(param).data(posCond);
        pstn = 3;
    else
        do = false; % not a subj field, dont add to all column.
    end
    if do
        allIdx=allIdx+1;
        pooled(allIdx,4) = E.(p).(param).data(allCond);
    end
    if pstn>0
       pos(allIdx)=pstn; 
    end
end
end

function isSig = test_significance( param, featTest )
% test significance
% 1.	W1 vs W5 for group average – is there a progressive change in parameter due to weight?
% 2.	X1 vs X5 for group average – is there a progressive change in parameter due to weight and exercise?
% 3.	W1 vs X1 for group average – is there a difference initially when weight is applied before and after exercise?
% 4.	W5 vs X5 for group average – Is the progressive effect of weight different before and after exercise?
% 5.	Is the W1-X1 difference significantly different than the X5-W5 difference?
phases  = {'U';'R';'W1';'W2';'W3';'W4';'W5';'X1';'X2';'X3';'X4';'X5';'P'};
featTest= split(featTest,':');
phase1  = featTest{1};
phase2  = featTest{2};

for i=1:length(phases)
   if strcmp(phase1, phases{i})
       idx1 = i;
   end
   if strcmp(phase2, phases{i})
       idx2 = i;
   end
end

if idx1==0 || idx2==0;keyboard;end

isSig   = cell(length(idx1),1);
alpha   = 0.05;
x       = param.data;
s       = param.sd;
n       = param.n;

if isempty(s) || isempty(n)
    return
end
mu1= x(idx1);   mu2= x(idx2);
s1 = s(idx1);   s2 = s(idx2);
n1 = n(idx1);   n2 = n(idx2);

if isnan(mu1) || isnan(mu2)
    p = 1;
elseif s1~=0 && s2~=0
    p = welch_t(mu1, s1, n1, mu2, s2, n2);
elseif (mu1==1&&s1==0) || (mu2==1&&s2==0) || (mu2==0&&s2==0) %testing against phase R
    df      = n1-1; % make sure df is the same type as X
    xmean   = mu1;
    sdpop   = s1;
    sqrtn   = sqrt(n1);
    xdiff   = (xmean - mu2);

    % Check for rounding issues causing spurious differences
    fix = (xdiff~=0) ...                                     % a difference
        & (abs(xdiff) < 100*sqrtn.*max(eps(xmean),eps(mu2)));  % but a small one
    if any(fix(:))
        % Fix any columns that are constant, even if computed difference is
        % non-zero but small
        constvalue = min(x,[],dim);
        fix = fix & all(x==constvalue | isnan(x),dim);
    end
    if any(fix(:))
        % Set difference and standard deviation to 0, and recompute mean
        xdiff(fix) = 0;
        sdpop(fix) = 0;
        xmean = xdiff+mu2;
    end
    ser = sdpop ./ sqrtn;
    tval = xdiff ./ ser;
    % Compute the correct p-value for the test, and confidence intervals
    % if requested.
    p = 2 * tcdf(-abs(tval), df);
else
    p=1;
end % end if

if p < alpha
    ast='*';
else
    ast='';
end % end if
    isSig = sprintf('%0.2f %s %0.2f%s(%.5f)', mu1,char(177),s1,ast,p);
%     isSig = sprintf('%0.2f %s %0.2f%s', mu1,char(177),s1,ast);
end % end function

% ----------------------------------------------------------------------- %

function outtable(D, tests, saveDir)
if ~strcmp(saveDir(end), '\')
    saveDir = horzcat(saveDir, '\');
end
fn      = fieldnames(D);
% RES     = D;
% col1    = {'U';'R';'W1';'W2';'W3';'W4';'W5';'X1';'X2';'X3';'X4';'X5';'P'};
% fn2     = {'TV','BF','FRCM','MINVNT'};
fn2     = {'FRCB','FRCM','FRCT'};
nrows   = length(tests);

for i = 1:length(fn2) % for each param
    tbl = cell(nrows, length(fn));
    for j = 1: length(fn) % for each participant
        if sum(contains( {'pos1','pos2','pos3','all'}, (fn{j}) )) == 1
            param = sprintf('%sREL',fn2{i});
        else
            param = fn2{i};
        end
        for r = 1:nrows
            featTest = tests{r};
            tbl{r,j} = test_significance( D.(fn{j}).(param), featTest );
        end
    end
    out = array2table( horzcat(tests,tbl), 'VariableNames', vertcat({'Condition'},fn));
    writetable(out, sprintf('%spaper_table_%s.csv', saveDir, fn2{i}));
end % end for i

end % end function

% ----------------------------------------------------------------------- % 

 function [D, header]= load_files(statdir)
 
cd(statdir);
files = ls;
D = struct;
have_header = false;
for i = 3: size(files, 1)
    f = files(i,:);
    file = strtrim(f);
%     tok= regexp(file, 'P._', 'match');
    tok = regexp(file, 'features', 'match');
    if length(tok)>0 && ~isempty(tok{1})
%         part_name=tok{1}(1:end-1);
        part_name = f(strfind(f,'P') : strfind(f,'P')+2);
    else
        continue
    end % end if
    if contains(file, '.csv')
        if ~have_header
            fid = fopen(file, 'r');
            header = textscan(fid, '%s', 1); 
            header = header{1};
            header = textscan(header{1}, '%s', 'Delimiter', ',');
            fclose(fid);
            header = header{1};
        end % end if
        D.(part_name)= readmatrix(file);
    else
        disp("Unrecognized file name!")
    end % end if
end % end for i

end % end function

% ----------------------------------------------------------------------- %

function E = get_position_data(D, header, params)
POS1    = 'P07,P08,P10,P11,P12'; % 5
% POS2    = 'P04,P05,P06,P09,P15,P16,P18'; % 7
POS2    = 'P04,P05,P06,P09,P15,P18'; % 6
% POS3    = 'P01,P02,P03,P13,P14,P17,P19'; % 7
POS3    = 'P01,P02,P03,P13,P17,P19'; % 6
par     = fieldnames(D);
E       = struct;
if nargin == 2
    params = header;
end
for i = 1:length(params)
    pos1Mu  = [];
    pos2Mu  = [];
    pos3Mu  = [];
    param = params{i};
    for j = 1:length(par)
        p = par{j};
        try
            E.(p).(param).data  = D.(p)(:,strcmp(header,param));
            E.(p).(param).sd    = D.(p)(:,strcmp(header,sprintf('%sSD',param)));
            E.(p).(param).n     = D.(p)(:,strcmp(header,sprintf('%sN',param)));
        catch
        end
        if contains(POS1, p)
            pos1Mu = [pos1Mu, D.(p)(:,strcmp(header,param))];
        end
        if contains(POS2, p)
            pos2Mu = [pos2Mu, D.(p)(:,strcmp(header,param))];
        end
        if contains(POS3, p)
            pos3Mu = [pos3Mu, D.(p)(:,strcmp(header,param))];
        end
    end
    E.pos1.(param).data = mean(pos1Mu, 2,'omitnan');
    E.pos1.(param).sd   = std( pos1Mu, [], 2,'omitnan');
    E.pos1.(param).n    = sum(~isnan(pos1Mu),2);
    
    E.pos2.(param).data = mean(pos2Mu, 2,'omitnan');
    E.pos2.(param).sd   = std( pos2Mu, [], 2,'omitnan');
    E.pos2.(param).n    = sum(~isnan(pos2Mu),2);
    
    E.pos3.(param).data = mean(pos3Mu, 2,'omitnan');
    E.pos3.(param).sd   = std( pos3Mu, [], 2,'omitnan');
    E.pos3.(param).n    = sum(~isnan(pos3Mu),2);
    
    E.all.(param).data  = mean([pos1Mu,pos2Mu,pos3Mu], 2,'omitnan');
    E.all.(param).sd    = std( [pos1Mu,pos2Mu,pos3Mu], [], 2,'omitnan');
    E.all.(param).n     = sum(~isnan([pos1Mu,pos2Mu,pos3Mu]),2);
end % end for
end % end function

% -------------------------------------------------------------------------

function tprob = welch_t(x1, s1, n1, x2, s2, n2)

tval    = (x1 - x2) / sqrt( s1^2 / n1 + s2^2 / n2 );
v       = (s1^2 / n1 + s2^2 / n2)^2 / ( s1^4/ (n1^2 * n1-1) + s2^4/ (n2^2 * n2-1) );
tdist2T = @(t,v) (1-betainc(v./(v+t.^2),v./2,0.5));    % 2-tailed t-distribution
tprob   = 1- tdist2T(tval,v);
% tdist1T = @(t,v) 1-(1-tdist2T(t,v))/2;              % 1-tailed t-distribution
% tprob   = 1- tdist1T(tval,v);
end % end function

% -------------------------------------------------------------------------

function results_figs(D)
SAVEDIR = 'C:\Users\Mark\Documents\GraduateStudies\LAB\EIT-restraint\zzMC\figures\paper\';
pp      = {'all','pos1','pos2','pos3'};
% relFeats= {'TVREL','BFREL','FRCMREL','MINVNTREL'};
% lab     = {'TV','BF','FRCM','MINVNT'};
relFeats= {'FRCBREL','FRCMREL','FRCTREL'};
lab     = {'FRCB','FRCM','FRCT'};
colours = {'#0072BD', '#D95319', '#EDB120', '#7E2F8E'};
markers = {'-o','--^',':s','-.d'};
% titles={'Normalized respiration rate', 'Normalized tidal volume', 'Delta FRC Middle Plane change from R as a ratio of mean R tidal volume'};
ylabels = {'Relative delta FRC (bottom plane)','Relative delta FRC (middle plane)','Relative delta FRC (top plane)'};
% ylabels = {'Relative tidal volume (RU * breath^{-1})','Relative respiration rate (breaths * min^{-1})','Relative delta FRC (middle plane)','Relative minute ventilation'};
lw      = 2;
cpsz    = 12;
mrkrsz  = 10;

for h = 1:length(relFeats)
    legEntries = [];
    fig = figure(h);
    fig.Units='normalized';fig.OuterPosition=[0 0 1 1];clf(); hold on;fig.PaperOrientation='landscape';
    ax1 = fig.Children;
    for i = 1:length(pp)
        mrkr    = markers{i};
        mrkrclr = colours{i};
        err     = D.(pp{i}).(relFeats{h}).sd ./ sqrt(D.(pp{i}).(relFeats{h}).n);
        if strcmp(pp{i}, 'all')
            errorbar(1,     D.(pp{i}).(relFeats{h}).data(1),     err(1),     mrkr, 'linewidth', lw, 'Color', mrkrclr, 'CapSize', cpsz, 'MarkerSize', mrkrsz);
            errorbar(2,     D.(pp{i}).(relFeats{h}).data(2),     err(2),     mrkr, 'linewidth', lw, 'Color', mrkrclr, 'CapSize', cpsz, 'MarkerSize', mrkrsz);
            errorbar(3:7,   D.(pp{i}).(relFeats{h}).data(3:7),   err(3:7),   mrkr, 'linewidth', lw, 'Color', mrkrclr, 'CapSize', cpsz, 'MarkerSize', mrkrsz);
%             plot(8:12,  D.(pp{i}).(relFeats{h}).data(8:12), mrkr, 'linewidth', lw, 'Color', mrkrclr, 'MarkerSize', mrkrsz);
            entry = errorbar(13,D.(pp{i}).(relFeats{h}).data(13),    err(13),    mrkr, 'linewidth', lw, 'Color', mrkrclr, 'CapSize', cpsz, 'MarkerSize', mrkrsz);
            txt1 = sprintf('n=%.0f',D.(pp{i}).(relFeats{h}).n(1));
            t1 = text(1,D.(pp{i}).(relFeats{h}).data(1)+err(1),txt1,'FontSize',16);
            txt2 = sprintf('n=%.0f',D.(pp{i}).(relFeats{h}).n(13));
            t2 = text(13,D.(pp{i}).(relFeats{h}).data(13)+err(13),txt2,'FontSize',16);
        else
            entry = errorbar(8:12,  D.(pp{i}).(relFeats{h}).data(8:12),  err(8:12),  mrkr, 'linewidth', lw, 'Color', mrkrclr, 'CapSize', cpsz, 'MarkerSize', mrkrsz);
        end
        legEntries = [legEntries, entry];
    end % end for i
    leg = legend(legEntries,...  
                '\begin{tabular}{p{.05cm}r}&group average\end{tabular}',...
                '\begin{tabular}{p{.05cm}r}&control \hspace{0.52cm}(arms at side)\end{tabular}',...
                '\begin{tabular}{p{.05cm}r}&restraint 1 (hands behind back)\end{tabular}',...
                '\begin{tabular}{p{.05cm}r}&restraint 2 (hands behind head)\end{tabular}',... 
                'FontSize', 15, 'Box', 'on', 'Location', 'best');
    set(leg,'interpreter','latex'); % set interpreter
        
    ax1.XTickLabels     = {'U', 'R', 'W_1','W_2','W_3','W_4','W_5','X_1','X_2','X_3','X_4','X_5','P'};
%     ax1.XTickLabels     = {'R', 'W_1','W_2','W_3','W_4','W_5','X_1','X_2','X_3','X_4','X_5'};
    ax1.XLim            = [0.5 length(ax1.XTickLabels)+.5];
    ax1.XTick           = 1:length(ax1.XTickLabels);
    ax1.XLabel.FontSize = 60;
    ax1.YLabel.FontSize = 60;
    set(gca,'XTickLabel', ax1.XTickLabels, 'fontsize', 15);
    ylabel(ylabels{h});
    xlabel('Experimental phase');
    
    t1.Units = 'normalized';
    t1.Position(1) = t1.Position(1)-t1.Extent(3)/2;
    t1.Position(2) = t1.Position(2)+t1.Extent(4);
    t2.Units = 'normalized';
    t2.Position(1) = t2.Position(1)-t2.Extent(3)/2;
    t2.Position(2) = t2.Position(2)+t2.Extent(4);
    printPDF(sprintf('%spaper_%s.pdf', SAVEDIR, lab{h}));
%     saveas(gcf, sprintf('%spaper_%s.svg', SAVEDIR, lab{h}));
end % end for h

end % end function

% -------------------------------------------------------------------------