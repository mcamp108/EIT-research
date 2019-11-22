function title_name= remove_underscores(file)
% Replace underscores with '-' and remove .eit suffix from file name
if strcmp(file(end-3:end), '.eit')
    file= file(1:end-4);
end % end if

title_name= "";
for c= file(1: length(file))
    if c== "_"
        title_name= title_name+ " ";
    else
        title_name= title_name+ c;
    end % end if
end % end for
title_name= char(title_name);
end % end function