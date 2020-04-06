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
fid = fopen(file, 'r');
txt = textscan(fid, '%s');
% txt{1}{1} = txt{1}{1}(4:end);
n_cols= length(regexp(txt{1}{1}, ',')) + 1;
col_spec=[];
for i=1:n_cols
    col_spec=[col_spec 'c'];
end

fprintf('\\begin{table}[H] \\footnotesize\n\\centering\n\\caption{}\n\\begin{tabular}{%s}\n', col_spec);
for i=1:size(txt{:},1)
    line = horzcat(txt{1}{i}, ' \\ \hline');
    rmchar= regexp(line, 'Â');
    if ~isempty(rmchar)
        line(rmchar)='';
    end % end if
%     j= length(line);
%     done=false; 
%     while ~done
    splitline = strsplit(line, ',');
    if length(splitline) > 1
        for k=1:length(splitline)-1
            fprintf('%s\t&\t', splitline{k});
        end % end for k
    end % end if
    fprintf('%s\n', splitline{k+1});
%         commaloc= regexp(line, ',', 'once');
%         if isempty(commaloc)
%             done=true;
%         else
%             if commaloc== 1
%                 line=fprintf('\t&\t%s\n', line(commaloc+1:end));
%             elseif commaloc== j
%                 line=fprintf('%s\t&\t', line(1:commaloc-1));
%             else
%                 line=fprintf('%s\t&\t%s', line(1:commaloc-1), line(commaloc+1:end));
%             end % end if
%             j=j+2;
%         end % end if
%     end % end while
%     disp(line);
end % end for
fprintf('\\end{tabular}\n\\caption*{}\n\\label{table:}\n\\end{table}\n');


end % end function