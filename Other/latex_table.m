function latex_table(file)
% -------------------------------------------------------------------------
% DESCRIPTION:
%   print table in latex format from CSV file.
% -------------------------------------------------------------------------
% PARAMETERS:
%   text or csv file.
% -------------------------------------------------------------------------   
% RETURNS:
%   prints table cells in latex format.
% -------------------------------------------------------------------------   
% AUTHOR:
%   Mark Campbell
%   Carleton University
%   markacampbell@cmail.carleton.ca
% -------------------------------------------------------------------------
% VERSION:
%   1.0.0
% -------------------------------------------------------------------------
%%
f = tdfread(file, ',');
cols = fieldnames(f);
n_cols = length(cols);
col_spec = [];

for i = 1:n_cols
    col_spec = [col_spec 'c'];
end

fprintf('\\begin{table}[H] \\footnotesize\n\\centering\n\\caption{}\n\\begin{tabular}{%s}\n', col_spec);

for i = 1:size(f.(cols{1}),1) % for each line
    line = '';
    for j = 1:n_cols - 1
        line = horzcat(line, strtrim(f.(cols{j})(i, :)), ' & ');
    end
    line = horzcat(line, strtrim(f.(cols{j})(i, :)), ' \\ \hline');

    fprintf('%s\n', line);
end % end for
fprintf('\\end{tabular}\n\\caption*{}\n\\label{table:}\n\\end{table}\n');


end % end function