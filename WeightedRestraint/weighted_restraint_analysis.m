% load data
D= load_files();
% fg1= figure('units','normalized','outerposition',[0 0 1 1]);
% header= {'BF','BFsd','BFn', 'FRC','FRCsd','FRCn', 'TV','TVsd','TVn','Ti_by_Tt','Ti_by_Ttsd','Ti_by_Ttn'};
header= {'Ti_by_Tt','BF','FRC', 'TV'};
recordings= {'sref', 'pref', 'wref', 'wepos', 'epos'};
participants= fieldnames(D);

for i=1:length(header)
    param=header{i};
    D.pos1.(param).data= mean([ D.P1(:,i), D.P2(:,i), D.P3(:,i)], 2);
    D.pos1.(param).sd= std([    D.P1(:,i), D.P2(:,i), D.P3(:,i)], [],2);
    D.pos2.(param).data= mean([ D.P4(:,i), D.P5(:,i), D.P6(:,i)], 2);
    D.pos2.(param).sd= std([    D.P4(:,i), D.P5(:,i), D.P6(:,i)], [],2);
    D.all.(param).data= mean([  D.P1(:,i), D.P2(:,i), D.P3(:,i), D.P4(:,i), D.P5(:,i), D.P6(:,i)], 2);
    D.all.(param).sd= std([     D.P1(:,i), D.P2(:,i), D.P3(:,i), D.P4(:,i), D.P5(:,i), D.P6(:,i)], [],2);
end

%% 1. for PW position, is there ate least one statistically significant
% decrease in FRC relative to PR after the addition of weight?
pos= 'pos1';
param= 'TV';
n=3;
sd= D.(pos).(param).sd;
mu0=ones(13,1);
mubar= D.(pos).(param).data;

z= (mubar - mu0) ./ (sd ./ sqrt(n));
disp(tcdf(z,n));
%% 2. Is this decrease in FRC accompanied by decrease in TV?
pos= 'pos2';
param= 'TV';
n=3;
sd= D.(pos).(param).sd;
mu0=ones(13,1);
mubar= D.(pos).(param).data;

z= (mubar - mu0) ./ (sd ./ sqrt(n));
disp(tcdf(z,n));

%% 3. Do they compensate by increasing BF?
pos= 'pos2';
param= 'BF';
n=3;
sd= D.(pos).(param).sd;
mu0=ones(13,1);
mubar= D.(pos).(param).data;

z= (mubar - mu0) ./ (sd ./ sqrt(n));
disp(1-tcdf(z,n));

%% does adding weight decrease FRC
for i=1:7
    x1= D.(participants{i})(2,4);
    s1= D.(participants{i})(2,5);
    n1= D.(participants{i})(2,6);

    x2= D.(participants{i})(6,4);
    s2= D.(participants{i})(6,5);
    n2= D.(participants{i})(6,6);

    p= welch_t(x1, s1, n1, x2, s2, n2);
    disp(p);
end
%% stats table
results=D;
var_names={};
results_tbl=[];
count=1;
for p=1:length(participants)
    for q=1:length(recordings)
        record=recordings{q};
        x= D.(participants{p}).(record)(:,1);
        results.(participants{p}).(record)=zeros(6,2);
        len_x=length(x);
        current_res= [];
        for r= 2:length(header)
            y= D.(participants{p}).(record)(:,r);
            P= polyfit(x,y,1);
            yfit = P(1)*x+P(2);
            change= round(100* (yfit(end)- yfit(1)) / yfit(end), 2);
            R= corrcoef(yfit, y);            
            t=R(1,2)* sqrt((len_x-2)/ (1-R(1,2)^2));
            Pv=round(2*(1- tcdf(t, len_x- 2)), 2);
            current_res(r-1,:)= [change, Pv];
        end
        results_tbl= [results_tbl,current_res];
        var_names{count}= horzcat(participants{p},'_change_',record);
        count= count+1;
        var_names{count}= horzcat('p_',num2str(count));
        count= count+1;
    end % end for q
end % end for p
table_out = array2table(results_tbl, 'VariableNames', var_names);
save_name= horzcat('pooled_results.csv');
writetable(table_out,save_name);
%%
% % normalize values to pref
% for p=1:length(participants)
%     norm_to= mean(D.(participants{p}).pref(:,2:end), 1);
%     for q=1:length(recordings)
%         record=recordings{q};
%         D.(participants{p}).(record)(:,2:end)= D.(participants{p}).(record)(:,2:end) ./ norm_to;
%     end % end for q
% end % end for p



%%
% for each person, how does parameter change with each recording?
part_pos= fieldnames(D);
for i= 1: numel(part_pos)
    part= part_pos{i};
    for j= 2:length(header)
        parameter_name= header{j};
        data_= [];
        label= [];
        for k= 1:length(recordings)
            record= recordings{k};
            data= D.(part).(record);
            data_= [data_; data(:, j)];
            m= repmat({record}, size(data, 1), 1);
            label= [label; m];
        end % end for k
        
        fig_title= remove_underscores(horzcat(part, ' feature - ', parameter_name));
        boxplot(data_, label);
        title(fig_title);
        cd 'C:\Users\Mark\Documents\GraduateStudies\LAB\EIT-restraint\zzMC\figures\feature_comparison';
        saveas(gcf, horzcat(fig_title, '.svg'));
    end % end for j
end % end for i
%%
for k= 1:length(recordings)
    record= recordings{k};
    
    for j= 2:length(header)
        parameter_name= header{j};
        data_= [];
        label= [];
        
        for i= 1: numel(part_pos)
            part= part_pos{i};
            data= D.(part).(record);
            data_= [data_; data(:, j)];
            m= repmat({part}, size(data, 1), 1);
            label= [label; m];
        end % end for k
        
        fig_title= remove_underscores(horzcat(record, ' feature - ', parameter_name));
        boxplot(data_, label);
        title(fig_title);
        saveas(fg1, horzcat(fig_title, '.svg'));
    end % end for j
end % end for i
 %%
% header= {'Time_s', 'Ti', 'Te', 'Ttotal', 'exp_to_exp', 'FRC', 'tidal_volume', 'breath_frequency'};
% recordings= {'sref', 'pref', 'wref', 'wepos', 'epos'};
select=(6:8);
% cd 'C:\Users\Mark\Documents\GraduateStudies\LAB\EIT-restraint\zzMC\figures\feature_comparison';
cd 'C:\Users\Mark\Documents\GraduateStudies\LAB\EIT-restraint\zzMC\figures\feature_comparison\ventilation plots';
for i=1:length(participants)
    plot_parameters(D.(participants{i}), header, select);
    sgtitle(participants{i});
    saveas(gcf,horzcat(participants{i}, ' vent_plots.svg'));
end % end for i
 

 function D= load_files()
 
cd 'C:\Users\Mark\Documents\GraduateStudies\LAB\EIT-restraint\zzMC\data\Mali Weighted Restraint\features';
files= ls;
D= struct;
for i= 3: size(files, 1)
    f=files(i,:);
    file= strtrim(f);
%     tok= regexp(file, 'P._', 'match');
    tok= regexp(file, 'features (', 'match');
    if length(tok)>0 && ~isempty(tok{1})
%         part_name=tok{1}(1:end-1);
        part_name=horzcat('P',f(end-5));
    else
        continue
    end % end if
    if contains(file, '.csv')
        D.(part_name)= readmatrix(file);
    else
        disp("Unrecognized file name!")
    end % end if
end % end for i

end % end function


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

function p= welch_t(x1, s1, n1, x2, s2, n2)

t= (x1 - x2) / sqrt( s1^2 / n1 + s2^2 / n2 );
v= (s1^2 / n1 + s2^2 / n2)^2 / ( s1^4/ (n1^2 * n1-1) + s2^4/ (n2^2 * n2-1) );
p = 1-tcdf(t,v);

end % end function

