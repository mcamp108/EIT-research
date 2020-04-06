% load data
bsln= 2;
statdir = 'C:\Users\Mark\Documents\GraduateStudies\LAB\EIT-restraint\zzMC\data\Mali Weighted Restraint\features\paper';
[D, header]= load_files(statdir);
E = get_position_data(D, header);
fg1= figure('units','normalized','outerposition',[0 0 1 1]);
results_figs(E);
% recordings= {'sref', 'pref', 'wref', 'wepos', 'epos'};
% participants= fieldnames(D);

% outtable(D, header);
%%
results_figs(D);
% %%
% figure;hold on;
% plot(D.all.TV.data);plot(D.pos2.TV.data);plot(D.pos3.TV.data);plot(D.P7.TV.data);
% plot(D.all.TV.data)
% plot(D.all.BF.data)
% 
% fn=fieldnames(D.pos3);
% for i=1:length(fn)
%     disp (fn{i});
%     disp(D.pos2.(fn{i}).data);disp(D.pos2.(fn{i}).sd);
% end
% 
% %% 1. for PW position, is there ate least one statistically significant
% % decrease in FRC relative to PR after the addition of weight?
% pos= 'p1';
% param= 'TVABS';
% n=3;
% v = repmat(2*n-2,13,1);
% sd0= repmat(D.(pos).(param).sd(2), 13, 1);
% mu0= repmat(D.(pos).(param).data(2), 13,1); % test different than prone rest.
% sdbar= D.(pos).(param).sd;
% mubar= D.(pos).(param).data;
% 
% tval= (mubar - mu0) ./ sqrt((sd0.^2 + sdbar.^2)./n);
% 
% tdist2T = @(t,v) (1-betainc(v./(v+t.^2),v./2,0.5));    % 2-tailed t-distribution
% tprob = 1- tdist2T(tval,v);
% disp(tprob);
% 
% %% 2. Is this decrease in FRC accompanied by decrease in TV?
% pos= 'pos2';
% param= 'TV';
% n=3;
% sd= D.(pos).(param).sd;
% mu0=ones(13,1);
% mubar= D.(pos).(param).data;
% 
% z= (mubar - mu0) ./ (sd ./ sqrt(n));
% disp(tcdf(z,n));
% 
% %% 3. Do they compensate by increasing BF?
% pos= 'p1';
% param= 'TVABS';
% n=D.(pos).(param).n;
% sd= D.(pos).(param).sd;
% mu0=repmat(D.(pos).(param).data(2),13,1);
% mubar= D.(pos).(param).data;
% 
% z= (mubar - mu0) ./ (sd ./ sqrt(n));
% disp(1-tcdf(z,n));
% 
% %%
% par2= {'p1','p2','p3','p4','p5','p6','p7'};
% param= 'FRCABS';
% allprob=zeros(13,7);
% alpha=0.05;
% txt= cell(13,7);
% for i=1:length(par2)
%     pos= par2{i};
%     
%     x1= D.(pos).(param).data;
%     s1= D.(pos).(param).sd;
%     n1= D.(pos).(param).n;
% 
%     x2= repmat(D.(pos).(param).data(2),13,1);
%     s2= repmat(D.(pos).(param).sd(2),13,1);
%     n2= repmat(D.(pos).(param).n(2),13,1);
% 
%     for j=1:13
%         if s1(j)~=0
%             allprob(j,i)= welch_t(x1(j), s1(j), n1(j), x2(j), s2(j), n2(j));
%         else
%             allprob(j,i)=1;
%         end % end if
%         txt{j,i}= horzcat(num2str(round(x1(j),2)),'(', num2str(round(s1(j),2)),')');
%         if allprob(j,i)<alpha
%             txt{j,i}=horzcat(txt{j,i},'*');
%         end % end if
%     end % end for j
%     
% end % end for i
% disp(txt);
% %% does adding weight decrease FRC
% for i=1:7
%     x1= D.(participants{i})(2,4);
%     s1= D.(participants{i})(2,5);
%     n1= D.(participants{i})(2,6);
% 
%     x2= D.(participants{i})(6,4);
%     s2= D.(participants{i})(6,5);
%     n2= D.(participants{i})(6,6);
% 
%     p= welch_t(x1, s1, n1, x2, s2, n2);
%     disp(p);
% end
% %% stats table
% results=D;
% var_names={};
% results_tbl=[];
% count=1;
% for p=1:length(participants)
%     for q=1:length(recordings)
%         record=recordings{q};
%         x= D.(participants{p}).(record)(:,1);
%         results.(participants{p}).(record)=zeros(6,2);
%         len_x=length(x);
%         current_res= [];
%         for r= 2:length(header)
%             y= D.(participants{p}).(record)(:,r);
%             P= polyfit(x,y,1);
%             yfit = P(1)*x+P(2);
%             change= round(100* (yfit(end)- yfit(1)) / yfit(end), 2);
%             R= corrcoef(yfit, y);            
%             t=R(1,2)* sqrt((len_x-2)/ (1-R(1,2)^2));
%             Pv=round(2*(1- tcdf(t, len_x- 2)), 2);
%             current_res(r-1,:)= [change, Pv];
%         end
%         results_tbl= [results_tbl,current_res];
%         var_names{count}= horzcat(participants{p},'_change_',record);
%         count= count+1;
%         var_names{count}= horzcat('p_',num2str(count));
%         count= count+1;
%     end % end for q
% end % end for p
% table_out = array2table(results_tbl, 'VariableNames', var_names);
% save_name= horzcat('pooled_results.csv');
% writetable(table_out,save_name);
% %%
% % for each person, how does parameter change with each recording?
% part_pos= fieldnames(D);
% for i= 1: numel(part_pos)
%     part= part_pos{i};
%     for j= 2:length(header)
%         parameter_name= header{j};
%         data_= [];
%         label= [];
%         for k= 1:length(recordings)
%             record= recordings{k};
%             data= D.(part).(record);
%             data_= [data_; data(:, j)];
%             m= repmat({record}, size(data, 1), 1);
%             label= [label; m];
%         end % end for k
%         
%         fig_title= remove_underscores(horzcat(part, ' feature - ', parameter_name));
%         boxplot(data_, label);
%         title(fig_title);
%         cd 'C:\Users\Mark\Documents\GraduateStudies\LAB\EIT-restraint\zzMC\figures\feature_comparison';
%         saveas(gcf, horzcat(fig_title, '.svg'));
%     end % end for j
% end % end for i

% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %


function isSig = test_significance( participant, header, param )

% test significance
% U vs R
% end X vs P
% each W time point vs R
% each X time point vs end W
isSig=cell(13,1);
alpha=0.05;

x1 = participant(:, strcmp(header, horzcat(param,'ABS')));
if strcmp(param, 'FIT')
    x1 = participant(:, strcmp(header, param));
end % end if
s1 = participant(:, strcmp(header, horzcat(param,'SD')));
n1 = participant(:, strcmp(header, horzcat(param,'N')));

x2= [repmat(x1(2),7,1) ; repmat(x1(7),5,1) ; x1(12)];
s2= [repmat(s1(2),7,1) ; repmat(s1(7),5,1) ; s1(12)];
n2= [repmat(n1(2),7,1) ; repmat(n1(7),5,1) ; n1(12)];

for j=1:13
    if s1(j)~=0 && s2(j)~=0
        p = welch_t(x1(j), s1(j), n1(j), x2(j), s2(j), n2(j));
    else
        p = 1;
    end % end if
    if p < alpha
        isSig{j}='*';
    else
        isSig{j}='';
    end % end if
end % end for j
isSig = {isSig};
end % end function

% ----------------------------------------------------------------------- %

function outtable(D, header)

fn = fieldnames(D);
RES = struct;
check = {'ABS', 'SD'};
col1 = {'U';'R';'W1';'W2';'W3';'W4';'W5';'X1';'X2';'X3';'X4';'X5';'P'};
for NPART=1:length(fn)
    for i=1:length(header)
        name = strtrim(header{i});
        if strcmp(name, 'FIT'); name = 'FITABS';end
        for j= 1:length(check)
            tok= regexp(name, check{j});
            if ~isempty(tok)
                param = name(1:tok-1);
                RES.(fn{NPART}).(param).(check{j}) = D.(fn{NPART})(:,i);
            end % end if
        end % end for j
    end % end for i
end % end for h

fn2 = fieldnames(RES.(fn{NPART}));
nrows = size(D.(fn{NPART})(:,i), 1);
for i=1:length(fn2) % for each param
    tbl = cell(nrows, NPART);
    for j=1:NPART % for each participant
        RES.(fn{j}).(fn2{i}).sig = test_significance( D.(fn{j}), header, fn2{i} );
        for r=1:nrows
            current = RES.(fn{j}).(fn2{i});
            tbl{r,j} = sprintf( '%0.2f(%0.2f)%s', current.(check{1})(r), current.(check{2})(r), current.sig{1}{r} );
        end % end for r
    end
    out = array2table( horzcat(col1,tbl), 'VariableNames', vertcat({'Condition'},fn));
    writetable(out, sprintf('paper_table_%s.csv', fn2{i}));
end % end for i

end % end function

% ----------------------------------------------------------------------- % 

 function [D, header]= load_files(statdir)
 
cd(statdir);
files= ls;
D= struct;
have_header=false;
for i= 3: size(files, 1)
    f=files(i,:);
    file= strtrim(f);
%     tok= regexp(file, 'P._', 'match');
    tok= regexp(file, 'features', 'match');
    if length(tok)>0 && ~isempty(tok{1})
%         part_name=tok{1}(1:end-1);
        part_name=f(1:2);
    else
        continue
    end % end if
    if contains(file, '.csv')
        if ~have_header
            fid = fopen(file, 'r');
            header = textscan(fid, '%s', 1); 
            header=header{1};
            header = textscan(header{1}, '%s', 'Delimiter', ',');
            fclose(fid);
            header=header{1};
        end % end if
        D.(part_name)= readmatrix(file);
    else
        disp("Unrecognized file name!")
    end % end if
end % end for i

end % end function

% ----------------------------------------------------------------------- %

function D= get_position_data(D, header)

    for i=1:length(header)
        param=header{i};
        D.pos3.(param).data= mean([ D.P1(:,i), D.P2(:,i), D.P3(:,i)], 2);
        D.pos3.(param).sd= std([    D.P1(:,i), D.P2(:,i), D.P3(:,i)], [],2);
        
        D.pos2.(param).data= mean([ D.P4(:,i), D.P5(:,i), D.P6(:,i)], 2);
        D.pos2.(param).sd= std([    D.P4(:,i), D.P5(:,i), D.P6(:,i)], [],2);
        
        D.pos1.(param).data=        D.P7(:,i);
        D.pos1.(param).sd= std(     D.P7(:,i), [], 2);
        
        D.all.(param).data= mean([  D.P1(:,i), D.P2(:,i), D.P3(:,i), D.P4(:,i), D.P5(:,i), D.P6(:,i), D.P7(:,i)], 2);
        D.all.(param).sd= std([     D.P1(:,i), D.P2(:,i), D.P3(:,i), D.P4(:,i), D.P5(:,i), D.P6(:,i), D.P7(:,i)], [],2);
    end % end for
par= {'P1','P2','P3','P4','P5','P6','P7'};
par2= {'p1','p2','p3','p4','p5','p6','p7'};
for i= 1:length(par)
    D.(par2{i}).BFABS.data=D.(par{i})(:,1);
    D.(par2{i}).BFABS.sd= D.(par{i})(:,2);
    D.(par2{i}).BFABS.n= D.(par{i})(:,3);
    
    D.(par2{i}).FRCABS.data= D.(par{i})(:,5);
    D.(par2{i}).FRCABS.sd= D.(par{i})(:,6);
    D.(par2{i}).FRCABS.n= D.(par{i})(:,7);
    
    D.(par2{i}).TVABS.data= D.(par{i})(:,9);
    D.(par2{i}).TVABS.sd= D.(par{i})(:,10);
    D.(par2{i}).TVABS.n= D.(par{i})(:,11);
end
end % end function

% -------------------------------------------------------------------------

function plot_parameters(data, variablenames, select)

recordings= {'sref', 'pref', 'wref', 'wepos', 'epos'};
timepoints= [0 1 2 3 4]*60;         n_timepoints=length(timepoints);
x_ticklabels=cell(1,13);            rows=[1,2,3,8,13];
n_variables=length(select);         colors= ['b','r','m'];
y=zeros(13,n_variables);            clf;                   

for i= 1:length(recordings)
    record= recordings{i};
    xax= data.(record)(:,1);
    hold on
    
    for j=1:n_variables
        idx=select(j);
        
        if i==1 || i==2 || i==5 %sref pref wepos
            y(rows(i),j)= mean(data.(record)(:,idx));
%             y(rows(i),j)= mean( data.(record)(:,idx) ./ data.(record)(:,4), 1);
            x_ticklabels{rows(i)}= record;
        else
            something= sum(xax>timepoints, 2);
            
            for k=0:n_timepoints-1
                row=rows(i)+k;
                y(row,j)= mean(data.(record)(something==k+1,idx));
%                 y(row,j)= mean( data.(record)(something==k+1,idx) ./ data.(record)(something==k+1,4), 1);
                x_ticklabels{row}= horzcat(record, ' ', num2str(k+1));
            end % end for k
            
        end % end if
        
    end % end for j
    
end % end for i
y= y./ max(y, [], 1);
for i=1:n_variables
    hold on    
    if length(select)== 2 && select== [2, 3]
        subplot(n_variables,1,1);plot(y(:,i),colors(i));
        xticks(1:13); xticklabels(x_ticklabels); xtickangle(45);
        legend('Inspiratory Time', 'Expiratory Time');
        xlabel('Position');ylabel('Proportion of Total breath time (%)');
    else
        subplot(n_variables,1,i);plot(y(:,i),colors(i));
        xticks(1:13); xticklabels(x_ticklabels); xtickangle(45);
        legend(horzcat('Normalized ', remove_underscores(variablenames{select(i)})));
        xlabel('Position');ylabel('Arbitrary Units');
    end
end % end for i

end % end function

% -------------------------------------------------------------------------

function tprob= welch_t(x1, s1, n1, x2, s2, n2)

tval= (x1 - x2) / sqrt( s1^2 / n1 + s2^2 / n2 );
v= (s1^2 / n1 + s2^2 / n2)^2 / ( s1^4/ (n1^2 * n1-1) + s2^4/ (n2^2 * n2-1) );

tdist2T = @(t,v) (1-betainc(v./(v+t.^2),v./2,0.5));    % 2-tailed t-distribution
tprob = 1- tdist2T(tval,v);

end % end function

% -------------------------------------------------------------------------

function results_figs(D)
cd 'C:\Users\Mark\Documents\GraduateStudies\LAB\EIT-restraint\zzMC\figures\paper';
pp={'all','pos1','pos2','pos3'};
abs={'BF','TV','FRCM'};
lab={'BF','TV','FRCM'};
titles={'Normalized respiration rate', 'Normalized tidal volume', 'Delta FRC change from R as a percentage of mean R tidal volume'};
ylabels={'Normalized respiration rate (breaths/min)', 'Normalized tidal volume (ml/breath)', 'Delta FRC (% change)'};
for h= 1:length(abs)
    fig=figure('units','normalized','outerposition',[0 0 1 1]);clf; hold on;
    fig.PaperOrientation='landscape';
    ax1=fig.Children;
    for i=1:length(pp)
        errorbar(1:length(D.(pp{i}).(abs{h}).data), D.(pp{i}).(abs{h}).data, -D.(pp{i}).(abs{h}).sd, D.(pp{i}).(abs{h}).sd, 'linewidth', 4,'CapSize',12);
    end % end for i
    legend('all','position 1', 'position 2','position 3');
    title(titles{h}, 'fontsize', 20);
    ax1.XLim=[0.5 13.5];    
    ax1.XTick=1:13;     
    ax1.XLabel.FontSize= 60;
    ax1.YLabel.FontSize= 60;
    ax1.XTickLabels={'U', 'R', 'W_1','W_2','W_3','W_4','W_5','X_1','X_2','X_3','X_4','X_5','P'};
    xtickangle(30);
    ylabel(ylabels{h});
    saveas(gcf, horzcat('paper_',lab{h}, '.svg'));
end % end for h
end % end function

% -------------------------------------------------------------------------