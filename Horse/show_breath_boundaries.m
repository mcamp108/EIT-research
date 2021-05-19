function show_breath_boundaries(data, triads)
% -------------------------------------------------------------------------
% DESCRIPTION:
%
% -------------------------------------------------------------------------
% PARAMETERS:
% 
% -------------------------------------------------------------------------   
% RETURNS:
% 
% -------------------------------------------------------------------------   
% AUTHOR:
%   Mark Campbell
%   Carleton University
%   markacampbell@sce.carleton.ca
% -------------------------------------------------------------------------
% VERSION:
%   1.0.0
% -------------------------------------------------------------------------
% (C) 2019-2021 Mark Campbell. License: GPL version 2 or version 3
% -------------------------------------------------------------------------

% show seconds
gsig = sum(data,1);
textY = min(gsig);
plot(gsig); hold on;
for i = 1:size(triads,1)
   plot(triads(i,:), gsig(triads(i,:)), 'or');
   text(triads(i,2), textY, num2str(i), 'FontSize',12);
end

end % end function