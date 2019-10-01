function title_name= remove_underscores(file)
% Replace underscores with '-' and remove .eit suffix from file name
title_name= "";
for c= file(1: length(file)-4)
    if c== "_"
        title_name= title_name+ " ";
    else
        title_name= title_name+ c;
    end % end if
end % end for
end % end function