% load data
cd 'C:\Users\Mark\Documents\GraduateStudies\LAB\WeightedRestraint\data\Mali Weighted Restraint';
dirs= ls;
D= struct;
for i= 3: size(ls, 1)
    folder= dirs(i, :);
    while strcmp(folder(end), ' ')
        folder= folder(1:end-1);
    end % end while
    if ~strcmp(folder, 'Other')
        cd(folder);
    else
        continue
    end % end for
    participant= folder(end-1:end);
    hdr= length(horzcat(folder, '_'))+ 1;
    files= ls;
    for f= 3: size(files, 1)
        file= files(f, :);
        while strcmp(file(end), ' ')
            file= file(1:end-1);
        end % end while
        if strcmp(file(end-3:end), '.csv')
            if strcmp(file(hdr:hdr+8), 'proneRef_')
                D.(participant).pref= readmatrix(file);
            elseif strcmp(file(hdr:hdr+13), 'proneRefWeight')
                D.(participant).wref= readmatrix(file);
            elseif strcmp(file(hdr:hdr+16), 'standingReference')
                D.(participant).sref= readmatrix(file);
            elseif strcmp(file(hdr:hdr+18), 'proneWeightExercise')
                D.(participant).wepos= readmatrix(file);
            elseif strcmp(file(hdr:hdr+20), 'proneNoWeightExercise')
                D.(participant).epos= readmatrix(file);
            else
                disp("Unrecognized file name!")
            end % end if
        end % end if
    end % end for
    cd ../
end % end for
%%
fn= fieldnames(D);
header= {   'ins_to_ins',               'exp_to_exp',               'ins_to_exp',...
            'del_z_lung_max',           'del_z_lung_min',           'left_lung_area',...
            'right_lung_area',          'left_lung_filling',        'right_lung_filling',...
            'left_lung_median_value',	'right_lung_median_value',  'global_inhomogeneity',...
            'center_of_ventilation'};
recordings= {'pref', 'wref', 'sref', 'wepos', 'epos'};


for i= 1:length(recordings)
    pos1= [];
    pos2= [];
    pos3= [];
    record= recordings{i};
    figure(i);clf;
    sgtitle(record);
    for j= 1:numel(fn)
        data= D.(fn{j}).(record);
        if j< 4
            pos3= [pos3; data];
        elseif j< 7
            pos2= [pos2; data];
        else
            pos1= [pos1; data];
        end % end if
    end % end for
    
    pos1(isnan(pos1))=0;
    pos2(isnan(pos2))=0;
    pos3(isnan(pos3))=0;
    
    p1_mean= mean(pos1, 1);
    p1_std= std(pos1, 1);
    p2_mean= mean(pos2, 1);
    p2_std= std(pos2, 1);
    p3_mean= mean(pos3, 1);
    p3_std= std(pos3, 1);
    for k= 1:length(p1_mean)
        if k==11
            keyboard;
        end
        subplot(4,4,k);
        p1_pd = makedist('Normal','mu',p1_mean(k),'sigma',p1_std(k));
        x= sort(pos1(:,k));
        p1_pdf= pdf(p1_pd, x);
        plot(x, p1_pdf);
        hold on
        p2_pd = makedist('Normal','mu',p2_mean(k),'sigma',p2_std(k));
        x= sort(pos2(:,k));
        p2_pdf= pdf(p2_pd, x);
        plot(x, p2_pdf);

        p3_pd = makedist('Normal','mu',p3_mean(k),'sigma',p3_std(k));
        x= sort(pos3(:,k));
        p3_pdf= pdf(p3_pd, x);
        plot(x, p3_pdf);
        xlabel(header{k});
%         legend('Posture 1', 'Posture 2', 'Posture 3');
        hold off
    end % end for
    
end % end for

%%

