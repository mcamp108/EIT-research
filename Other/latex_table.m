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
%
f = tdfread(file, ',');
cols = fieldnames(f);
n_cols = length(cols);
col_spec = [];

for i = 1:n_cols
    col_spec = [col_spec 'c'];
end

fprintf('\\begin{table}[!ht] \\footnotesize\n\\centering\n\\begin{tabular}{%s}\n', col_spec);
% first row
line = '';
for i = 1:n_cols % for each line
    characters = cols{i};
    if isnumeric(characters)
        characters = num2str(round(characters,4));
    end
    if i < n_cols
        line = horzcat(line, strtrim(characters), ' & ');
    else
%         line = horzcat(line, strtrim(characters), ' \\ \hline');
        line = horzcat(line, strtrim(characters), ' \\');
    end
end

fprintf('%s\n', line);

% subsequent rows
for i = 1:size(f.(cols{1}),1) % for each line
    line = '';
    for j = 1:n_cols
        characters = f.(cols{j})(i, :);
        try characters = str2double(characters);catch;end
        if isnan(characters);characters = f.(cols{j})(i, :);end;
        if isnumeric(characters)
            characters = num2str(round(characters,4));
        end
        if j < n_cols
            line = horzcat(line, strtrim(characters), ' & ');
        else
%             line = horzcat(line, strtrim(characters), ' \\ \hline');
            line = horzcat(line, strtrim(characters), ' \\');
        end
    end
    fprintf('%s\n', line);
end % end for
fprintf('\\end{tabular}\n\\caption{}\n\\label{table:}\n\\end{table}\n');

end % end function