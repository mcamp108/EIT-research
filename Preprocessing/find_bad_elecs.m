function scores= find_bad_elecs(data, imdl)
% -------------------------------------------------------------------------
% DESCRIPTION:
%   clip_scores= find_bad_elecs(data, imdl)
% give electrodes a score based on number of clipped measurements occuring
% while that electrode was measuring. Divided by total number of frames in
% data. A score of > 1 is considered bad.
% -------------------------------------------------------------------------
% PARAMETERS:
% 
%
%
% -------------------------------------------------------------------------   
% RETURNS:
% 
%
%
% -------------------------------------------------------------------------   
% AUTHOR:
%   Mark Campbell
%   Carleton University
%   markacampbell@cmail.carleton.ca
% -------------------------------------------------------------------------
% VERSION:
%   1.0.0
% -------------------------------------------------------------------------

scores= elec_clip_scores(data, imdl);

end % end function

function clip_scores= elec_clip_scores(data, imdl)

msel= imdl.fwd_model.meas_select;
% find electrodes responsible for each selected measurement
mm = find(msel);
elec_used = zeros(length(mm),2);
count=1;
for i=1:32
    [r,c]= find(imdl.fwd_model.stimulation(i).meas_pattern);
    for j=1:max(r)
       elec_used(count,:)= c(r==j)';
       count= count+ 1;
    end
end % end for

% find total clipped data points for each selected measurement
data_= data(mm,:);
clipped= (abs(imag(data_)) >= 1.7e-3) + (abs(data_) >= 1.7e-3);
clipped(clipped>1) = 1;
tlt_clipped= sum(clipped, 2);

% give electrodes a score based on number of clipped measurements occuring
% while that electrode was measuring. Divided by total number of frames in
% data. A score of > 1 is considered bad.

score= 1;
clip_scores= zeros(32,1);
while score> 0
    score= max(tlt_clipped);
    worst= find(tlt_clipped==score);
    tlt_clipped(worst)= 0;
    for i=1:length(worst)
        elecs= elec_used(worst(i),:);
        if ~isempty(elecs)
            clip_scores(elecs)= clip_scores(elecs)+ score;
        end % end if
    end % end for
end % end while

clip_scores= clip_scores./size(data,2);

end % end function