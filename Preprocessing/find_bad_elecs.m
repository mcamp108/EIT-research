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

% find electrodes responsible for each selected measurement
msel = imdl.fwd_model.meas_select;
nElec = size(imdl.fwd_model.electrode, 2);
mm = find(msel);
elecUsed = zeros(length(mm),2);
count = 1;
nMeasWithElec = zeros(32, 1);

for i = 1:nElec
    [r,c] = find(imdl.fwd_model.stimulation(i).meas_pattern);
    for j = 1:max(r)
       elecUsed(count,:) = c(r == j)';
       count = count+ 1;
    end
end % end for

for i = 1:nElec
    nMeasWithElec(i) = sum(elecUsed == i, 'all');
end % end for

data_= data(mm,:);
clipped = find_clipped_agresive(data_);
tltClipped = sum(clipped, 2);

% give electrodes a score based on number of clipped measurements occuring
% while that electrode was measuring. Divided by total number of frames in
% data. A score of > 1 is considered bad.
score = 1;
clip_scores = zeros(32,1);

while score > 0
    score = max(tltClipped);
    worst = find(tltClipped == score);
    tltClipped(worst) = 0;
    
    for i=1:length(worst)
        elecs = elecUsed(worst(i),:);
    
        if ~isempty(elecs)
            clip_scores(elecs) = clip_scores(elecs) + score;
        end % end if
    
    end % end for

end % end while

clip_scores = clip_scores ./ (size(data,2) * nMeasWithElec);

end % end function

function clipped = find_clipped(data)
% find total clipped data points for each selected measurement Draw a
% circle and slowly decrease the radius. Tally how many points are outside
% of that circle. When the number of points outside the circle drops
% sharply, we can be confident that all clipped data points are outside of
% this circle. 
% TODO: Needs to be tested on data without clipping.

ALPHA = 1.0;
DELTA = 0.01;
ITERS = 90;

outsideCirc = zeros(ITERS, 1);
distFromC = sqrt( real(data).^2 + imag(data).^2 );
maxDist = max(distFromC, [], 'all');

for i=1:ITERS
    ALPHA = ALPHA - DELTA;
    outsideCirc(i) = sum(distFromC >= maxDist * ALPHA, 'all');
end % end for

foundInIter = diff(outsideCirc);
% stoppingIter = find(foundInIter == max(foundInIter)) + 2;
stoppingIter = find(diff(foundInIter) > 0, 1, 'first');
clipDist = maxDist * (1 - (stoppingIter * DELTA) );
clipped = distFromC >= clipDist;

end % end function


function clipped = find_clipped_agresive(data)
% find total clipped data points for each selected measurement Draw a
% circle and slowly decrease the radius. Tally how many points are outside
% of that circle. When the number of points outside the circle drops
% sharply, we can be confident that all clipped data points are outside of
% this circle. 
% TODO: Needs to be tested on data without clipping.

ALPHA = 0;
DELTA = 0.05;
ITERS = 21;

insideCirc = zeros(ITERS, 1);
distFromC = sqrt( real(data).^2 + imag(data).^2 );
maxDist = max(distFromC, [], 'all');

for i = 1:ITERS
    ALPHA = ALPHA + DELTA;
    insideCirc(i) = sum(distFromC < maxDist * ALPHA, 'all');
end % end for

foundInIter = diff(insideCirc);
% stoppingIter = find(foundInIter == max(foundInIter)) + 2;
stoppingIter = find(diff(foundInIter) > 0, 1, 'last');
clipDist = maxDist * (stoppingIter * DELTA);
clipped = distFromC >= clipDist;

end % end function