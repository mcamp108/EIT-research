function filesOut = load_wr_data(path)
    thing   = strfind(path, '\');
    folder  = path(thing(end)+1 : end);
    hdr     = length(horzcat(folder, '_'))+ 1;
    files   = ls(path);
    filesOut = cell(5,1);
    for f = 3: size(files, 1)
        file = files(f, :);
        FILE = lower(file);
        while strcmp(file(end), ' ')
            file = file(1:end-1);
        end % end while
        if contains(FILE, 'standingreference') || contains(FILE, 'minstanding.eit')
            filesOut{1} = file;
        elseif contains(FILE, 'proneref.eit') || contains(FILE, 'minprone.eit')
            filesOut{2} = file;
        elseif contains(FILE, 'pronerefweight') || contains(FILE, 'minwithweight.eit')
            filesOut{3} = file;
        elseif contains(FILE, 'proneweightexercise') || contains(FILE, 'minwithwafterexercise.eit')
            filesOut{4} = file;
        elseif contains(FILE, 'pronenoweightexercise')
            filesOut{5} = file;
        else
            fprintf('\nUnrecognized: %s', file);
        end % end if
    end % end for f
end % end function