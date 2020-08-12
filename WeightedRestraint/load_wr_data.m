function files= load_wr_data(folder)
    
while strcmp(folder(end), ' ')
    folder= folder(1:end-1);
end % end while

hdr= length(horzcat(folder, '_'))+ 1;
cd(folder);
files= ls;
for f= 3: size(files, 1)
    file= files(f, :);
    while strcmp(file(end), ' ')
        file= file(1:end-1);
    end % end while
    if strcmp(file(end-3:end), '.eit')
        if strcmp(file(hdr:hdr+8), 'proneRef.')
            pref= file;
        elseif strcmp(file(hdr:hdr+13), 'proneRefWeight')
            wref= file;
        elseif strcmp(file(hdr:hdr+16), 'standingReference')
            sref= file;
        elseif strcmp(file(hdr:hdr+18), 'proneWeightExercise')
            wepos= file;
        elseif strcmp(file(hdr:hdr+20), 'proneNoWeightExercise')
            epos= file;
        else
            disp("Unrecognized file name!")
        end % end if
    end % end if
end % end for f

files= {sref, pref, wref, wepos, epos};
% cd ../.
end % end function