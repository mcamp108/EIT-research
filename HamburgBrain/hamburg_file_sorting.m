
%% 1st layer
cd 'C:\Users\Mark\Documents\GraduateStudies\LAB\HamburgBrain\Data\12.2\DICOM 12.2\12-2_MRI';
files= ls;
files= files(3:end, :);
folder_names= files(:, 4:7);

last_fol= '';
for i= 1:size(folder_names, 1)
    fol= folder_names(i, :);
    if ~strcmp(last_fol, fol)
        mkdir(fol);
        last_fol= fol;
        for f= 1: size(files, 1)
            file= files(f, :);
            if strcmp(last_fol, file(4:7))
                movefile(file, last_fol);
            end
        end
    end
end
%% second layer

cd 'C:\Users\Mark\Documents\GraduateStudies\LAB\HamburgBrain\Data\12.2\DICOM 12.2\12-2_MRI';
folders= ls;
folders= folders(3:end, :);
for j= 1:size(folders, 1)
    k= folders(j,:);
    cd(k);
    
    % now in subfolder
    files= ls;
    files= files(3:end, :);
    if (size(files, 2)<= 16) || (strcmp('0001', files(end, 9:12)))
    else
        folder_names= files(:, 9:12);
        last_fol= '';
        for i= 1:size(folder_names, 1)
            fol= folder_names(i, :);
            if ~strcmp(last_fol, fol)
                mkdir(fol);
                last_fol= fol;
                files_to_move=[];
                for f= 1: size(files, 1)
                    file= files(f, :);
                    if strcmp(last_fol, file(9:12))
                        movefile(file, last_fol);
                    end % end if
                end % end for
            end % end if
        end % end for
    end % end if
    cd ../
end



%%
% sorting subdirectory

cd 'C:\Users\Mark\Documents\GraduateStudies\LAB\HamburgBrain\Data\8.2\DICOM 8.2\8-2_MRI\11_brainMRgrainy';
files= ls;
files= files(3:end, :);
folder_names= files(:, 9:12);

last_fol= '';
for i= 1:size(folder_names, 1)
    fol= folder_names(i, :);
    if ~strcmp(last_fol, fol)
        dirname= "set " + fol;
        mkdir(dirname);
        last_fol= fol;
        for f= 1: size(files, 1)
            file= files(f, :);
            if strcmp(fol, file(14:17))
                movefile(file, dirname);
            end
        end
    end
end



