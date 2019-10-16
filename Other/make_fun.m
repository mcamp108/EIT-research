function [] = make_fun(V)

fun_name = horzcat('function out= ', V, '(in)');
end_function= 'end % end function';
str_date=horzcat('%   ', datestr(date,'dd.mmm.yyyy'));
fid = fopen('C:\Users\Mark\Documents\GraduateStudies\LAB\Scripts\EIT-reserach\Other\my_template.m');
F = fread(fid, '*char')';
fclose(fid);
text= [fun_name newline F newline str_date newline newline newline end_function];
fid = fopen([V '.m'], 'w');
fwrite(fid, text);
fclose(fid);
open([V '.m'])

end % end function