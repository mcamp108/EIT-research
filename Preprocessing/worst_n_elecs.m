function [worst, scores] = worst_n_elecs(data, imdl, params)
% -------------------------------------------------------------------------
% DESCRIPTION:
% [worst, scores] = worst_n_elecs(data, imdl, params)
%
% Finds the noisiest n electrodes across all data files in D based on
% number of clipped measurements occuring while that electrode was
% measuring. The scores represent the proportion of total measurements made
% using that electrode that were clipped or faulty. Scores range from 0 to
% 1.
% -------------------------------------------------------------------------
% PARAMETERS:
%   dataIn (cell):
%       contains one or more EIT data arrays
%   imdl (EIDORS inverse model object):
%       Inverse model object
%   params (struct):
%       required fields:
%           n (int):      
%               Number of electrodes to be excluded
%           type (string):   
%               'elec' to find the worst n electrodes, or 'meas' to find
%               the worst n measurements
% -------------------------------------------------------------------------   
% RETURNS:
%   worst:
%       An array containing the worst n electrodes or measurements
%       (specifiable) by number, sorted from worst to best.
%   scores:
%       The score associated with each electrode or measurement in worst,
%       in the same order as worst. Proportion of measurements that were
%       faulty when the electrode was measuring, or proportion of faulty
%       measurements for each individual measurement.
% -------------------------------------------------------------------------   
% AUTHOR:
%   Mark Campbell
%   Carleton University
%   markacampbell@cmail.carleton.ca
% -------------------------------------------------------------------------
% VERSION:
%   2.0.0
% -------------------------------------------------------------------------
% (C) 2019-2021 Mark Campbell. License: GPL version 2 or version 3
% -------------------------------------------------------------------------
params = parse_options(params, imdl);
n = params.n;
nFiles = length(data);
elecScores = zeros(params.nElecs, nFiles);
mmScores = zeros( params.nMeas, nFiles );

for i = 1: nFiles
    [elecScores, mmScores(:, i)] = elec_clip_scores(data{i}, imdl, params);
    elecScores(:, i) = resolve_scores(elecScores, imdl, params);
end % end for

if nFiles > 1
    elecScores = mean(elecScores, 2);
    mmScores = mean(mmScores, 2);
end % end if

if strcmp(params.type, 'elec')
    [hi_lo_scores, order] = sort(elecScores, 'descend');
elseif strcmp(params.type, 'meas')
    [hi_lo_scores, order] = sort(mmScores, 'descend');
end % end if

worst = order(hi_lo_scores > 0);
numFound = min( [length(worst), n] );

worst = worst(1: numFound);
scores = hi_lo_scores(1: numFound);

end % end function

% ======================================================================= %

function resolved = resolve_scores(elecScores, imdl, params)
    % If one electrode is fauly, its partners will appear faulty as well.
    % This function will detect if one electrode is making its parterns
    % look bad, and resolve its partners' scores
    resolved = elecScores;
    nElecs = params.nElecs;
    stimPairs = zeros(nElecs, 2);

    for i = 1: nElecs
        stimPairs(i, 1) = find(imdl.fwd_model.stimulation(i).stim_pattern < 0);
        stimPairs(i, 2) = find(imdl.fwd_model.stimulation(i).stim_pattern > 0);    
    end % end for

    for i = 1: nElecs
        [r, c] = find(stimPairs == i);
        partners = unique( stimPairs(r, c) );
        partners = partners(partners ~= i);
        for j = 1: length(partners)
            p = partners(j);
            [r, c] = find(stimPairs == p);
            jPartners = unique( stimPairs(r, c) );
            jPartners = jPartners(jPartners ~= p);
            jPartner = jPartners(jPartners ~= i);
            if elecScores(i) >= ( elecScores(p) + elecScores(jPartner) )
                resolved(p) = elecScores(jPartner);
            end % end if
        end % end for
    end % end for
end % end function

% ======================================================================= %

function [elec_scores, mmScores] = elec_clip_scores(data, imdl, params)
    % A score is the proportion of measurements that are clipped or have
    % negative in-phase component. Scores range from 0 to 1.
    %
    % elec_scores:
    %     Proportion of faulty measurements while the electrode was
    %     measuring
    % mmScores:
    %     The score for each selected measurement (specified by
    %     fwd_model.meas_select)
    
    nElecs = params.nElecs;
    mm = params.mm;
    elecUsed = zeros(params.nMeas, 2);
    nMeasWithElec = zeros(nElecs, 1);
    stimPairs = zeros(nElecs, 2);
    
    % Find total number of measurements per frame for each electrode
    count = 1;
    for i = 1: nElecs
        [r, c] = find(imdl.fwd_model.stimulation(i).meas_pattern);
        stimPairs(i, 1) = find(imdl.fwd_model.stimulation(i).stim_pattern < 0);
        stimPairs(i, 2) = find(imdl.fwd_model.stimulation(i).stim_pattern > 0);
        for j = 1:max(r)
           elecUsed(count, :) = c(r == j)';
           count = count + 1;
        end % end for
    end % end for

    for i = 1: nElecs
        nMeasWithElec(i) = sum(elecUsed == i, 'all');
    end % end for

    data_ = data(mm, :);
    
    % mark clipped measurements.
    clipped = find_clipped_agresive(data_);     
    
    % mark measurements with negative in-phase component.
    rDataNeg = real(data_) < 0;                 
    
%     plot(data_(~rDataNeg), '+g'); hold on;
%     plot(data_(rDataNeg), '+r');
%     axis equal;
%     xlabel('In-phase');
%     ylabel('Quadrature');
    
    % Calculate total score
    allMarked = sum( (clipped * 1 + rDataNeg * 1) > 0, 2);
    mmScores = allMarked ./ size(data_, 2);

    score = 1;
    elec_scores = zeros(nElecs, 1);
    while score > 0
        score = max(allMarked);
        worst = find(allMarked == score);
        allMarked(worst) = 0;
        for i = 1: length(worst)
            elecs = elecUsed(worst(i), :);
            if ~isempty(elecs)
                elec_scores(elecs) = elec_scores(elecs) + score;
            end % end if
        end % end for
    end % end while
    elec_scores = elec_scores ./ (size(data, 2) * nMeasWithElec);
end % end function

% ======================================================================= %

function clipped = find_clipped_agresive(data)
    % find total clipped data points for each selected measurement Draw a
    % circle and slowly decrease the radius. Tally how many points are
    % outside of that circle. When the number of points outside the circle
    % drops sharply, we can be confident that all clipped data points are
    % outside of this circle.

    ALPHA = 0;
    DELTA = 0.05;
    ITERS = 21;

    insideCirc = zeros(ITERS, 1);
    distFromC = sqrt( real(data) .^ 2 + imag(data) .^ 2 );
    maxDist = max(distFromC, [], 'all');

    for i = 1: ITERS
        ALPHA = ALPHA + DELTA;
        insideCirc(i) = sum(distFromC < maxDist * ALPHA, 'all');
    end % end for

    foundInIter = diff(insideCirc);
    stoppingIter = find(diff(foundInIter) > 0, 1, 'last');
    clipDist = maxDist * (stoppingIter * DELTA);
    clipped = distFromC >= clipDist;
    
%     % demo figures
%     ALPHA = 0;
%     for i = 1: ITERS-1
%         ALPHA = ALPHA + DELTA;
%         insideIdx = distFromC < maxDist * ALPHA;
%         inside = data(insideIdx);
%         outside = data(~insideIdx);
%         close all;
%         figure('units','normalized','outerposition',[-1.0042 0.2259 0.8083 0.7778]);
%         
%         subplot(2, 2, [1, 3]);
%         plot(real(outside), imag(outside), '+k');
%         hold on;
%         plot(real(inside), imag(inside), '+g');
%         axis equal;
%         xlabel('In-phase');
%         ylabel('Quadrature');
%         
%         subplot(2, 2, 2);
%         semilogy(insideCirc(1: i), '-ok');
%         xlabel('Iteration number');
%         xlim([0, ITERS + 1]);
%         ylabel('Points inside circle');
%         ylim([min(insideCirc)-10, max(insideCirc)+10]);
%         
%         subplot(2, 2, 4);
%         semilogy(foundInIter(1: i), '-ok');
%         xlabel('Iteration number');
%         xlim([0, ITERS + 1]);
%         ylabel('Points found this iteration');
%         ylim([0.1, max(foundInIter)]);
%         
%         ttl = sprintf('frame_%0.f.png', i);
%         
%         pause(1);
%         saveas(gcf(), horzcat('C:\Users\Mark\Documents\EIT\proj\2021\EIT2021_conf\eqadr\', ttl));
%         
%     end
%     
%     good = data(~clipped);
%     bad = data(clipped);
%     subplot(2, 2, [1, 3]);
%     cla();
%     plot(real(bad), imag(bad), '+r');
%     hold on;
%     plot(real(good), imag(good), '+g');
%     axis equal;
%     xlabel('In-phase');
%     ylabel('Quadrature');
% 
%     subplot(2, 2, 4);
%     xline(stoppingIter, '--');
%     ttl = sprintf('frame_%0.f.png', i+1);
%     saveas(gcf(), horzcat('C:\Users\Mark\Documents\EIT\proj\2021\EIT2021_conf\eqadr\', ttl));
    
end % end function


function params = parse_options(optIn, imdl)
    assert(isfield(imdl, 'solve_use_matrix'), "The inverse model must contain the 'solve_use_matrix' field. Re-generate the inverse model after setting the 'keep_intermediate_results' field in the inverse model options struct to true.");
    params = optIn;
    params.nElecs = length(imdl.fwd_model.electrode);
    params.mm = find(imdl.fwd_model.meas_select);
    params.nMeas = sum(imdl.fwd_model.meas_select);
    if ~isfield(optIn, 'thresh')
        params.thresh = 0.25;
    end % end if
    if ~isfield(optIn, 'type')
        params.type = 'elec';
    end % end if
    if ~isfield(optIn, 'n')
        params.n = inf;
    end % end if
    if strcmp(params.type, 'elec')
        params.n = min( [params.nElecs, params.n] );
    elseif strcmp(params.type, 'meas')
        params.n = min( [params.nMeas, params.n] );
    end % end if
end