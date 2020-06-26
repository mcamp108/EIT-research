dir = 'E:\University\Graduate Studies\Hamburg Data\DICOM 9.2.-20200522T202320Z-001\DICOM 9.2\Ta_Airchill_7_19990101';
% 1st layer
cd(dir);
files = ls;
files = files(3:end, :);
folder_names = files(:, 4:7);
done = false;
last_fol = '';
while ~done
    key = folder_names(1, :);
    fol = num2str( str2double(key) );
    if ~strcmp(last_fol, fol)
        mkdir( fol );
        last_fol = fol;
        for f = 1: size(files, 1)
            file = files(f, :);
            if strcmp(key, file(4:7))
                movefile(file, last_fol);
            else
                files = files(f:end, :);
                folder_names = folder_names(f:end, :);
                break
            end % end if
        end % end for
        
        if f == size(folder_names, 1)
            done = true;
        end % end if
    
    end % end if

end % end while
%%
% sorting subdirectory
cd(dir);
folders = ls;
folders = folders(3:end, :);

for i = 1:size(folders, 1)
    cd(folders(i, :));
    files = ls;
    files = files(3:end, :);
    if length(files(1,:)) == 21
        folder_names = files(:, 14:17);
        for j = 1:size(folder_names, 1)
            fol = folder_names(j, :);
            file = files(j, :);
            if exist(horzcat(pwd,'\',fol), 'dir') ~= 7
                mkdir(fol);
            end % end if
            if strcmp(fol, file(14:17))
                movefile(file, fol);
            end % end if   
        end % end for
    end % end if
    cd(dir);
end
