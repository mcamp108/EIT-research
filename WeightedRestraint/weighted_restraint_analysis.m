% load data
cd 'C:\Users\Mark\Documents\GraduateStudies\LAB\WeightedRestraint\data\Mali Weighted Restraint';
dirs= ls;
D= struct;
for i= 3: size(ls, 1)
    folder= strtrim(dirs(i, :));

    if ~strcmp(folder, 'Other')
        cd(folder);
    else
        continue
    end % end for
    
    participants= folder(end-1:end);
    hdr= length(horzcat(folder, '_'))+ 1;
    files= ls;
    for f= 3: size(files, 1)
        file= strtrim(files(f, :));
        if contains(file, '_features.csv')
            if contains(file, 'proneRef_')
                D.(participants).pref= readmatrix(file);
            elseif contains(file, 'proneRefWeight')
                D.(participants).wref= readmatrix(file);
            elseif contains(file, 'standingReference')
                D.(participants).sref= readmatrix(file);
            elseif contains(file, 'proneWeightExercise')
                D.(participants).wepos= readmatrix(file);
            elseif contains(file, 'proneNoWeightExercise')
                D.(participants).epos= readmatrix(file);
            else
                disp("Unrecognized file name!")
            end % end if
        end % end if
    end % end for f
    cd ../
end % end for i
%%
participants= fieldnames(D);
header= {   'Time_s', 'ins_to_ins',               'exp_to_exp',               'ins_to_exp',...
                    'del_z_lung_max',           'del_z_lung_min',           'left_lung_filling',...          
                    'right_lung_filling',       'left_lung_median_value',	'right_lung_median_value',...  
                    'global_inhomogeneity',     'center_of_ventilation'};
recordings= {'sref', 'pref', 'wref', 'wepos', 'epos'};

% sort by posture
for i= 1:length(recordings)
    record= recordings{i};
    D.pos1.(record)= [];
    D.pos2.(record)= [];
    D.pos3.(record)= [];
    
    for j= 1: numel(participants)
        data= D.(participants{j}).(record);
        if j< 4
            D.pos3.(record)= [D.pos3.(record); data];
        elseif j< 7
            D.pos2.(record)= [D.pos2.(record); data];
        else
            D.pos1.(record)= [D.pos1.(record); data];
        end % end if
    end % end for j
end % end for i

%%
% for each person, how does parameter change with each recording?
part_pos= fieldnames(D);
figure('units','normalized','outerposition',[0 0 1 1]);
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
        
        fig_title= remove_underscores(horzcat(part, ' feature- ', parameter_name));
        boxplot(data_, label);
        title(fig_title);
        print_convert(horzcat(fig_title, '.png'));
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
        
        fig_title= remove_underscores(horzcat(record, ' feature- ', parameter_name));
        boxplot(data_, label);
        title(fig_title);
        print_convert(horzcat(fig_title, '.png'));
    end % end for j
end % end for i
 