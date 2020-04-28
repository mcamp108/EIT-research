function [wElecs, scores, mmScores, stimScores] = worst_n_elecs(D, imdl, n)
% -------------------------------------------------------------------------
% DESCRIPTION:
%   [wElecs, scores, mmScores, stimScores] = worst_n_elecs(D, imdl, n)
%
% Finds the noisiest n electrodes across all data files in D based on
% number of clipped measurements occuring while that electrode was
% measuring. The scores represent the proportion of total measurements made
% using that electrode that were clipped or faulty. Scores range from 0 to
% 1.
% -------------------------------------------------------------------------
% PARAMETERS:
%   D:
%       Either a struct where each field is an EIT data matrix or an array
%       consisting of EIT data from a single file.
% -------------------------------------------------------------------------   
% RETURNS:
%   wElecs:
%       An array containing the worst n electrodes by number, sorted from
%       worst to best.
%   scores:
%       The score associated with each electrode in wElecs, in the same
%       order as wElecs.
% -------------------------------------------------------------------------   
% AUTHOR:
%   Mark Campbell
%   Carleton University
%   markacampbell@cmail.carleton.ca
% -------------------------------------------------------------------------
% VERSION:
%   1.1.0
% -------------------------------------------------------------------------

% find worst n electrodes across all sequences
data = {};
if isstruct(D)
    fn = fieldnames(D);
    n_files = length(fn);
    for i= 1:n_files
       data{i} = D.(fn{i}).eit.data;
    end % end for i
elseif iscell(D)
    n_files = length(D);
    data = D;
else
    n_files = 1;
    data = D;
end % end if

elec_scores = zeros(32, n_files);
mmScores = zeros( length(find(imdl.fwd_model.meas_select)), n_files );
stimScores = zeros(32, n_files);

for i= 1:n_files
    if n_files == 1
        [elecScores, mmScores(:,i), stimScores(:,i)] = elec_clip_scores(data, imdl);
    else
        [elecScores, mmScores(:,i), stimScores(:,i)] = elec_clip_scores(data{i}, imdl);
    end % end if
    elec_scores(:,i) = resolve_scores(elecScores, imdl);
%         plot(elecScores); hold on; 
%         plot(stimScores(:,i));
%         plot(elec_scores(:,i));
%         keyboard;
end % end for

if n_files > 1
    elec_scores= mean(elec_scores, 2);
    mmScores= mean(mmScores, 2);
    stimScores= mean(stimScores, 2);
end % end if

[hi_lo_scores, elecs]= sort(elec_scores, 'descend');
wElecs = elecs(hi_lo_scores > 0);

if length(wElecs) > n
    wElecs = wElecs(1:n);
end % end if

if nargout > 1
    scores = hi_lo_scores(1:n);
end % end if

end % end function


function resolved = resolve_scores(elecScores, imdl)
% If one electrode is fauly, its partners will appear faulty as well. This
% function will detect if one electrode is making its parterns look bad,
% and resolve its partner's scores
resolved = elecScores;
nElec = size(imdl.fwd_model.electrode, 2);
stimPairs = zeros(nElec, 2);
for i = 1:nElec
    stimPairs(i, 1) = find(imdl.fwd_model.stimulation(i).stim_pattern == -1);
    stimPairs(i, 2) = find(imdl.fwd_model.stimulation(i).stim_pattern == 1);    
end % end for

for i = 1:nElec
    [r,c] = find(stimPairs == i);
    partners = unique( stimPairs(r,c) );
    partners = partners(partners ~= i);
    for j = 1:length(partners)
        p = partners(j);
        [r,c] = find(stimPairs == p);
        jPartners = unique( stimPairs(r,c) );
        jPartners = jPartners(jPartners ~= p);
        jPartner = jPartners(jPartners ~= i);
        if elecScores(i) >= ( elecScores(p) + elecScores(jPartner) )
            resolved(p) = elecScores(jPartner);
        end % end if
    end % end for
end % end for

end % end function


% ======================================================================= %

function [clip_scores, mmScores, stimScores] = elec_clip_scores(data, imdl)
% clip_scores:
%     give electrodes a score based on proportion of clipped measurements
%     occuring while that electrode was measuring. Score ranges from 0 to
%     1.
% mmScores:
%     For each selected measurement, score is given based on portion of total
%     measurements that were clipped (scores range from 0 to 1).

    msel = imdl.fwd_model.meas_select;
    nElec = size(imdl.fwd_model.electrode, 2);
    mm = find(msel);
    elecUsed = zeros(length(mm),2);
    count = 1;
    nMeasWithElec = zeros(nElec, 1);
    stimPairs = zeros(nElec, 2);
    
    for i = 1:nElec
        [r,c] = find(imdl.fwd_model.stimulation(i).meas_pattern);
        stimPairs(i, 1) = find(imdl.fwd_model.stimulation(i).stim_pattern == -1);
        stimPairs(i, 2) = find(imdl.fwd_model.stimulation(i).stim_pattern == 1);
        for j = 1:max(r)
           elecUsed(count,:) = c(r == j)';
           count = count+ 1;
        end % end for
    end % end for

    for i = 1:nElec
        nMeasWithElec(i) = sum(elecUsed == i, 'all');
    end % end for

    data_= data(mm,:);
    clipped = find_clipped_agresive(data_);
    tltClipped = sum(clipped, 2);
    mmScores = sum(clipped, 2) ./ size(data_, 2); % output
    measPerPair = length(mm) / nElec;

    score = 1;
    clip_scores = zeros(nElec,1);
    stimScores = zeros(nElec,1); % number of clipped measurements when this electrode is measuring
    
    for i = 1:nElec
        stop  = i * measPerPair;
        start = stop - measPerPair + 1;
        idx = stimPairs(i, :);
        stimScores( idx ) = stimScores( idx ) + sum( tltClipped(start:stop) );
    end % end for
    stimScores = stimScores ./ (measPerPair * size(data, 2));
    
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

% ======================================================================= %

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



