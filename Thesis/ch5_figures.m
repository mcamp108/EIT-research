% Chapter 5 Figures
%--------------------------------------------------------------------------
figure('units','normalized','outerposition',[0 0 1 1]); clf;
SAVEDIR = 'C:\Users\Mark\Documents\GraduateStudies\LAB\Thesis\files\figures\Chap_5\';
pig= '11-2';
close all
ELECSCORETHRESH = 0.25;
[fmdl, imdl] = get_pig_mdl(pig);
eit_files= {    'EIT_11.2_Nativ_1_Schleuse.eit',                'EIT_11.2_Nativ_2_ZVK.eit',... 
                'EIT_Sequenz_3_Schleuse.eit',                   'EIT_Seqeunz4_zvk.eit',...
                'EIT_11.2_4h_nach_Stroke_Schleuse_Sequ5.eit',   'EIT_4h_nach_Stroke_ZVk_Sequ6.eit'};
perf_files= {   'EIT_11.2_Nativ_1_Schleuse.mat',                'EIT_11.2_Nativ_2_ZVK.mat',...
                'EIT_Sequenz_3_Schleuse.mat',                   'EIT_Seqeunz4_zvk.mat',...
                'EIT_11.2_4h_nach_Stroke_Schleuse_Sequ5.mat',   'EIT_4h_nach_Stroke_ZVk_Sequ6.mat'};  
sync1= [1225,   211,    261,    55,     495,    169]; sync2= [6369,   5694,   7008,   5814,   5322,   6666];
e_apn= [-10.2,  10.7,   25.1,   22.3,   8.1,    12.1]; e_inj= [11.5,   15.3,   32.2,   66,     18.2,   16.2];
e_vnt= [89.1,   91.9,   113.7,  104.4,  80.4,   85.8]; p_apn= [11690,  10714,  36930,  31092,  16000,  19419];
p_inj= [33376,  15269,  43994,  74955,  26149,  23516]; p_vnt= [111004, 91873,  125489, 113239, 88339,  93187];
D = hamburg_load_data(eit_files, perf_files, sync1, sync2, e_apn, e_inj, e_vnt, p_apn, p_inj, p_vnt, pig, imdl);
%--------------------------------------------------------------------------
% find worst 6 electrodes
[rmElecs, scores, mmScores, ~] = worst_n_elecs(D, imdl, 6);
rmMeas = mmScores >= ELECSCORETHRESH;
fprintf('measurements removed: %s \n', num2str(sum(rmMeas)));
imdl_comp = comp_RM_bad_elec(imdl, rmMeas, 'meas');
% Load data
ref = 'self';
fn = fieldnames(D);
lowpass_cutoff = 4;

for i=1
    dd = D.(fn{i}).eit.data;
    absZ = sqrt( real(dd).^2 + imag(dd).^2 );
    D.(fn{i}).eit.fdata = lowpass_iir(absZ', lowpass_cutoff, D.(fn{i}).eit.fs)';
    D.(fn{i}).eit.fdata = lowpass_iir(D.(fn{i}).eit.data', lowpass_cutoff, D.(fn{i}).eit.fs)';
end % end for
useData = D.seq1.eit.fdata;

% Set colour limits
imdl_comp_resolve = comp_RM_bad_elec(imdl, 1);
imdl_comp_noResolve = comp_RM_bad_elec(imdl, [1 6 28]);
imgr_nocomp = inv_solve(imdl, mean(useData,2), useData);
imgr_resolve = inv_solve(imdl_comp_resolve, mean(useData,2), useData);
imgr_noResolve = inv_solve(imdl_comp_noResolve, mean(useData,2), useData);

clim = -inf;
for i = [imgr_nocomp, imgr_resolve, imgr_noResolve]
    clim = max( clim, max( i.elem_data(:) ) );
end
for i = [imgr_nocomp, imgr_resolve, imgr_noResolve]
    i.calc_colours.clim = clim;
end

%--------------------------------------------------------------------------
% imgr_nocomp(isnan(imgr_nocomp))= 0;
% imgr_resolve(isnan(imgr_resolve))= 0;
% imgr_noResolve(isnan(imgr_noResolve))= 0;

%--------------------------------------------------------------------------        
% make figures
frm = 200;
imgr_compare = imgr_resolve;
imgr_compare.elem_data = [  imgr_nocomp.elem_data(:,frm), imgr_noResolve.elem_data(:,frm),  imgr_resolve.elem_data(:,frm), zeros(3172,1)];
slices_compare = calc_colours( calc_slices(imgr_compare) );
bkg = slices_compare(:,:,1)==1;
slices_compare(:,:,4) = slices_compare(:,:,3) - slices_compare(:,:,2) + 125;
for i=4
    temp = slices_compare(:,:,i); temp(bkg) = 1;
    slices_compare(:,:,i) = temp;
end
final_img = mk_mosaic(slices_compare,[1, 1],[], 4);
final_img(isnan(final_img)) = 0;
image(final_img);axis equal; axis image; axis tight;
colorbar;
fg = gcf();
ax = fg.Children(2);
ax.XTick = [1,3,5,7]*32+1;
ax.XTickLabel = {'No Compensation';'No Score Resolution';'Score Resolution';'Resolution Difference'};
ax.YTick = 0;
ax.XAxis.FontSize=20;
% ax.YAxis.FontSize=20; ax.YAxis.FontWeight='bold'; ax.YAxis.TickLabelRotation=-45;
printPDF(sprintf('%sResolveImgs',SAVEDIR));

% =========================================================================

function D = hamburg_load_data(eit_files, perf_files, sync1, sync2, e_apn, e_inj, e_vnt, p_apn, p_inj, p_vnt, pig, imdl)
    RMTHRESHOLD=0.25;
    % Field names
    fn = {'seq1', 'seq2', 'seq3', 'seq4', 'seq5', 'seq6'};
    n_files = length(eit_files);
    for i = 1:n_files
        % load data
        [vv,aux] = eidors_readdata(eit_files{i});
        eitfs = median(1e6/median(diff(aux.t_rel))); % framerate
%         [ dd, D.(field).aux ] = eidors_readdata(file);
        % assign struct fields
        D.(fn{i}).eit.data= vv;
        impedanceFactor =  2.048 / (2^12 * 0.173 * 0.003) / 2^15; % = 0.9633 / 2^15;
        D.(fn{i}).eit.elec_impedance= impedanceFactor* aux.elec_impedance;
        D.(fn{i}).name = horzcat( num2str(i), '-', char(remove_underscores(eit_files{i})) );
        D.(fn{i}).pig = pig;
        D.(fn{i}).eit.fs = eitfs;
        D.(fn{i}).eit.sync1 = sync1(i);
        D.(fn{i}).eit.sync2 = sync2(i);
        D.(fn{i}).eit.apn = D.(fn{i}).eit.sync1 + round(e_apn(i)* eitfs);
        D.(fn{i}).eit.inj = D.(fn{i}).eit.sync1 + round(e_inj(i)* eitfs);
        D.(fn{i}).eit.vnt = D.(fn{i}).eit.sync1 + round(e_vnt(i)* eitfs);
        D.(fn{i}).perf = load( perf_files{i} );
        D.(fn{i}).perf.apn = p_apn(i);
        D.(fn{i}).perf.inj = p_inj(i);
        D.(fn{i}).perf.vnt = p_vnt(i);

        % allign eit and perfusion sequences, trim syncronization spikes
        D.(fn{i}) = allign_eit_and_perf(D.(fn{i}));
        D.(fn{i}) = find_perf_landmarks(D.(fn{i}));
        
        inspect_eit_elec_and_data({D.(fn{i}).eit.data, aux}, imdl, RMTHRESHOLD);
%         printPDF('C:\Users\Mark\Documents\GraduateStudies\LAB\Thesis\files\figures\Chap_5\noResolve.pdf');
    end % end for

end % end function

% ======================================================================= %

function seq= allign_eit_and_perf(seq)
    % Find optimal alignment between eit and perfusion data
    % Annotate peaks and troughs of perfusion signal and transfer annotations
    % to EIT data
    trim = 0;
    switch seq.name
        case '4-EIT nach Perfusionsminderun 2 9.2'; trim = 15;
        case '3-EIT nach Embolisation 1 9.2'; trim = 6;
        case '2-EIT 9.2 Nativ 2'; trim = 5;
        case '6-EIT nach 6h Perfusionsminderung 10.02 ZVK Sequ6'; trim = 6;
        case '6-EIT 4h nach Stroke ZVk Sequ6'; trim = 2;
        case '2-EIT 12.2 NATIV ZVK Sequ 2'; trim = 5;
%         case
    end % end switch
    estart = seq.eit.sync1 + seq.eit.fs; % desired starting frame is 1 s after first sync
    eend = seq.eit.sync2 - seq.eit.fs; % desired end frame is 1 s before last sync
    
    % downsample perf data to 100 samples per second
    seq.perf.data = seq.perf.data( 1:10:length(seq.perf.data) );
    seq.perf.tickrate = seq.perf.tickrate / 10;
    seq.perf.apn = round( seq.perf.apn / 10 );
    seq.perf.inj = round( seq.perf.inj / 10 );
    seq.perf.vnt = round( seq.perf.vnt / 10 );
    
    pstart= round( (estart/ seq.eit.fs)* seq.perf.tickrate );
    pend=   round( (eend/ seq.eit.fs)* seq.perf.tickrate );

    eapn= seq.eit.apn/ seq.eit.fs;
    einj= seq.eit.inj/ seq.eit.fs;
    evnt= seq.eit.vnt/ seq.eit.fs;
    papn= seq.perf.apn/ seq.perf.tickrate;
    pinj= seq.perf.inj/ seq.perf.tickrate;
    pvnt= seq.perf.vnt/ seq.perf.tickrate;

    t1= eapn- papn;
    t2= einj- pinj;
    t3= evnt- pvnt;
    shift= mean( [t1, t2, t3] );

    if shift> 0 % EIT sequence is ahead of perfusion sequence
        ecut= round( estart + shift * seq.eit.fs);
        pcut= pstart;   
    elseif shift< 0 % Perfusion sequence is ahead of EIT sequence
        ecut= estart;
        pcut= pstart+ abs(round( shift*  seq.perf.tickrate));
    end % end if
    eend = round(eend);
    ecut = round( ecut + trim * seq.eit.fs );
    pcut = round( pcut + trim * seq.perf.tickrate);

    seq.eit.data=   seq.eit.data(:, ecut+1: eend);

    seq.eit.sync1 = seq.eit.sync1- ecut;
    seq.eit.sync2 = seq.eit.sync2- ecut;
    seq.eit.apn =   seq.eit.apn- ecut;
    seq.eit.inj =   seq.eit.inj- ecut;
    seq.eit.vnt =   seq.eit.vnt- ecut;

    seq.perf.data = movmean( seq.perf.data(:, pcut + 1: pend + pcut), 15);
    seq.perf.apn =  seq.perf.apn - pcut;
    seq.perf.inj =  seq.perf.inj - pcut;
    seq.perf.vnt =  seq.perf.vnt - pcut;
    
end % end function

% ======================================================================= %

function imgr= get_imgr(fdata, ref, imdl)
    
    vv= real(fdata);
    imgr= inv_solve(imdl, ref, vv);
    imgr.calc_colours.ref_level= 0;

end % end function

% ======================================================================= %

% function [w_elecs, scores] = worst_n_elecs(D, imdl, n)
% 
%     % find worst n electrodes across all sequences
%     fn= fieldnames(D);
%     n_files=length(fn);
%     elec_scores= zeros(32, n_files);
% 
%     for i= 1:n_files
%         elec_scores(:,i)= find_bad_elecs( D.(fn{i}).eit.data, imdl );
%     end % end for
% 
%     elec_scores= mean(elec_scores, 2);
%     [hi_lo_scores, elecs]= sort(elec_scores, 'descend');
%     w_elecs = elecs(hi_lo_scores > 0);
% 
%     if length(w_elecs) > n
%         w_elecs = w_elecs(1:n);
%     end % end if
% 
%     if nargout == 2
%         scores = hi_lo_scores(1:n);
%     end % end if
% 
% end % end function

% ======================================================================= %

function seq= find_perf_landmarks(seq)
    
    FS = seq.perf.tickrate;
    
    % Find peaks and valleys in perfusion Signal
    perf_data = seq.perf.data;
    perf_data = perf_data- mean(perf_data); % zero data
    
    if strcmp(seq.name, '4-EIT nach Perfusionsminderun 2 9.2')
        v_search = peakfinder(perf_data, 5.5, [], -1);
    else
        v_search = peakfinder(perf_data, 8, [], -1);
    end % end if
    
    if v_search(1) == 1
        v_search = v_search(2:end);
    end % end if
    
    % Make sure we start at a true perfusion valley
    while perf_data(v_search(1)) - perf_data(v_search(2)) > 10
        v_search = v_search(2:end);
    end % end if
    
    v_idx = zeros(1, length(v_search));
    v_idx(1) = v_search(1);
    
    if strcmp(seq.name, '4-EIT nach Perfusionsminderun 2 9.2')
        takeIdx = v_search >= 130;
        v_search = v_search(takeIdx == 1);
        v_search = v_search(diff(v_search) > 30);
    end % end if
    
    j = 1;
    alpha = 5;
    if strcmp(seq.pig, '9-2')
        BETA = 10;
        THETA = 60;
    else
        BETA = 20;
    end % end if
    
    % Avoid adding higher false valleys
    for i = 2: length(v_search)
        if abs( perf_data(v_idx(j)) - perf_data(v_search(i)) ) < 10
            j = j + 1;
            v_idx(j) = v_search(i);
        % sometime signal wanders downwards. If the height between is not
        % too large and the timing is within acceptable bounds of a cardiac
        % cycle, accept it.

        elseif ( j > 3 ) && ( abs( perf_data(v_idx(j)) - perf_data(v_search(i)) ) < BETA ) && ( v_search(i) - v_idx(j) >= min(diff(v_idx(v_idx~=0))) - alpha )
            j = j + 1;
            v_idx(j) = v_search(i);
        end % end if
    end % end for
    
    v_idx = v_idx(v_idx ~= 0);
    
    % Remove valleys within the synchronization flush
    invalids = [seq.perf.apn - FS, seq.perf.apn + FS;...
                seq.perf.inj - FS, seq.perf.inj + FS;...
                seq.perf.vnt - FS, seq.perf.vnt + FS];
    for i = 1:size(invalids, 1)
        rm = find( (v_idx > invalids(i,1)) + (v_idx < invalids(i,2)) == 2 );
        if ~isempty(rm)
            v_idx(rm) = 0;
        end % end if
    end % end if
    v_idx = v_idx(v_idx ~= 0);
        
    % Remove valleys that are too far apart
    valdif = diff(v_idx);
    maxdif = mean(valdif) + std(valdif);
    if std(valdif) > mean(valdif)
        maxdif = quantile(valdif,4);
        maxdif = maxdif(end);
    end % end if
    v_idx(find(valdif > maxdif) + 1) = 0;
    v_idx = v_idx(v_idx ~= 0);
    
    % prepare ouput
    valdif = diff(v_idx);
    seq.perf.peaks = zeros(1, length(v_idx) - 1);
    seq.perf.vals = zeros(2,length(v_idx) - 1);
    for i = 1:length(v_idx) - 1
        % exclude peaks that would be within perfusion flushes
        if valdif (i) < maxdif
            if strcmp(seq.name, '4-EIT nach Perfusionsminderun 2 9.2')
                dd = perf_data(v_idx(i)+ THETA: v_idx(i+1) );
                idx = find(dd == max(dd), 1) + v_idx(i)+THETA - 1;
                if length(idx) > 1
                    idx = idx(2);
                elseif isempty(idx)
                    dd = perf_data(v_idx(i): v_idx(i+1) );
                    idx = find(dd == max(dd), 1) + v_idx(i) - 1;
                end % end if
            else
                dd = perf_data(v_idx(i): v_idx(i+1));
                idx = find(dd == max(dd), 1) + v_idx(i) - 1;
            end % end if
            seq.perf.vals(:,i) = [v_idx(i), v_idx(i+1)];
            seq.perf.peaks(i) = idx;
        end
    end % end for i
    seq.perf.vals = seq.perf.vals(:, seq.perf.vals(1,:) ~= 0);
    seq.perf.peaks = seq.perf.peaks(:, seq.perf.peaks(1,:) ~= 0);
    
    % Transfer these annotations to the EIT data
    seq.eit.peaks=  round( (seq.perf.peaks ./ FS).* seq.eit.fs );
    seq.eit.vals=   round( (seq.perf.vals ./ FS).* seq.eit.fs );
    
    if seq.eit.peaks(1) < 1
        seq.eit.peaks(1)= 1;
    end % end if
    
    if seq.eit.vals(1) < 1
        seq.eit.vals(1)= 1;
    end % end if
% % testing
% figure; plot(perf_data); hold on; plot(seq.perf.peaks, perf_data(seq.perf.peaks), 'o');
end % end function

% ======================================================================= %

function glbl_ref = get_glbl_ref(D, pig)
    
    switch pig
        
        case '8-2'
            glbl_ref = mean( D.seq1.eit.fdata( :, D.seq1.eit.apn: D.seq1.eit.inj ), 2 );
        case '9-2'
            glbl_ref = mean( D.seq2.eit.fdata( :, D.seq1.eit.apn: D.seq1.eit.inj ), 2 );
        case '10-2'
            glbl_ref = mean( D.seq1.eit.fdata( :, D.seq1.eit.apn: D.seq1.eit.inj ), 2 );
        case '11-2'
            glbl_ref = mean( D.seq2.eit.fdata( :, D.seq2.eit.apn: D.seq2.eit.inj ), 2 );
        case '12-2'
            glbl_ref = mean( D.seq1.eit.fdata( :, D.seq2.eit.apn: D.seq1.eit.inj ), 2 );
    end % end switch
        
end % end function

% ======================================================================= %

function varargout = peakfinder(x0, sel, thresh, extrema, includeEndpoints, interpolate)

    % https://www.mathworks.com/matlabcentral/fileexchange/25500-peakfinder-x0-sel-thresh-extrema-includeendpoints-interpolate
    narginchk(1, 6);
    nargoutchk(0, 2);
    s = size(x0);
    flipData =  s(1) < s(2);
    len0 = numel(x0);
    if len0 ~= s(1) && len0 ~= s(2)
        error('PEAKFINDER:Input','The input data must be a vector')
    elseif isempty(x0)
        varargout = {[],[]};
        return;
    end % end if
    if ~isreal(x0)
        warning('PEAKFINDER:NotReal','Absolute value of data will be used')
        x0 = abs(x0);
    end % end if
    if nargin < 2 || isempty(sel)
        sel = (max(x0)-min(x0))/4;
    elseif ~isnumeric(sel) || ~isreal(sel)
        sel = (max(x0)-min(x0))/4;
        warning('PEAKFINDER:InvalidSel',...
            'The selectivity must be a real scalar.  A selectivity of %.4g will be used',sel)
    elseif numel(sel) > 1
        warning('PEAKFINDER:InvalidSel',...
            'The selectivity must be a scalar.  The first selectivity value in the vector will be used.')
        sel = sel(1);
    end % end if
    if nargin < 3 || isempty(thresh)
        thresh = [];
    elseif ~isnumeric(thresh) || ~isreal(thresh)
        thresh = [];
        warning('PEAKFINDER:InvalidThreshold',...
            'The threshold must be a real scalar. No threshold will be used.')
    elseif numel(thresh) > 1
        thresh = thresh(1);
        warning('PEAKFINDER:InvalidThreshold',...
            'The threshold must be a scalar.  The first threshold value in the vector will be used.')
    end % end if
    if nargin < 4 || isempty(extrema)
        extrema = 1;
    else
        extrema = sign(extrema(1)); % Should only be 1 or -1 but make sure
        if extrema == 0
            error('PEAKFINDER:ZeroMaxima','Either 1 (for maxima) or -1 (for minima) must be input for extrema');
        end % end if
    end % end if
    if nargin < 5 || isempty(includeEndpoints)
        includeEndpoints = true;
    end
    if nargin < 6 || isempty(interpolate)
        interpolate = false;
    end % end if
    x0 = extrema*x0(:); % Make it so we are finding maxima regardless
    thresh = thresh*extrema; % Adjust threshold according to extrema.
    dx0 = diff(x0); % Find derivative
    dx0(dx0 == 0) = -eps; % This is so we find the first of repeated values
    ind = find(dx0(1:end-1).*dx0(2:end) < 0)+1; % Find where the derivative changes sign
    % Include endpoints in potential peaks and valleys as desired
    if includeEndpoints
        x = [x0(1);x0(ind);x0(end)];
        ind = [1;ind;len0];
        minMag = min(x);
        leftMin = minMag;
    else
        x = x0(ind);
        minMag = min(x);
        leftMin = min(x(1), x0(1));
    end % end if
    % x only has the peaks, valleys, and possibly endpoints
    len = numel(x);
    if len > 2 % Function with peaks and valleys
        % Set initial parameters for loop
        tempMag = minMag;
        foundPeak = false;
        if includeEndpoints
            % Deal with first point a little differently since tacked it on
            % Calculate the sign of the derivative since we tacked the first
            %  point on it does not neccessarily alternate like the rest.
            signDx = sign(diff(x(1:3)));
            if signDx(1) <= 0 % The first point is larger or equal to the second
                if signDx(1) == signDx(2) % Want alternating signs
                    x(2) = [];
                    ind(2) = [];
                    len = len-1;
                end
            else % First point is smaller than the second
                if signDx(1) == signDx(2) % Want alternating signs
                    x(1) = [];
                    ind(1) = [];
                    len = len-1;
                end % end if
            end % end if
        end % end if
        % Skip the first point if it is smaller so we always start on a
        %   maxima
        if x(1) >= x(2)
            ii = 0;
        else
            ii = 1;
        end % end if
        % Preallocate max number of maxima
        maxPeaks = ceil(len/2);
        peakLoc = zeros(maxPeaks,1);
        peakMag = zeros(maxPeaks,1);
        cInd = 1;
        % Loop through extrema which should be peaks and then valleys
        while ii < len
            ii = ii+1; % This is a peak
            % Reset peak finding if we had a peak and the next peak is bigger
            %   than the last or the left min was small enough to reset.
            if foundPeak
                tempMag = minMag;
                foundPeak = false;
            end % end if
            % Found new peak that was lager than temp mag and selectivity larger
            %   than the minimum to its left.
            if x(ii) > tempMag && x(ii) > leftMin + sel
                tempLoc = ii;
                tempMag = x(ii);
            end % end if
            % Make sure we don't iterate past the length of our vector
            if ii == len
                break; % We assign the last point differently out of the loop
            end % end if
            ii = ii+1; % Move onto the valley
            % Come down at least sel from peak
            if ~foundPeak && tempMag > sel + x(ii)
                foundPeak = true; % We have found a peak
                leftMin = x(ii);
                peakLoc(cInd) = tempLoc; % Add peak to index
                peakMag(cInd) = tempMag;
                cInd = cInd+1;
            elseif x(ii) < leftMin % New left minima
                leftMin = x(ii);
            end % end if
        end % end while
        % Check end point
        if includeEndpoints
            if x(end) > tempMag && x(end) > leftMin + sel
                peakLoc(cInd) = len;
                peakMag(cInd) = x(end);
                cInd = cInd + 1;
            elseif ~foundPeak && tempMag > minMag % Check if we still need to add the last point
                peakLoc(cInd) = tempLoc;
                peakMag(cInd) = tempMag;
                cInd = cInd + 1;
            end % end if
        elseif ~foundPeak
            if x(end) > tempMag && x(end) > leftMin + sel
                peakLoc(cInd) = len;
                peakMag(cInd) = x(end);
                cInd = cInd + 1;
            elseif tempMag > min(x0(end), x(end)) + sel
                peakLoc(cInd) = tempLoc;
                peakMag(cInd) = tempMag;
                cInd = cInd + 1;
            end % end if
        end % end if
        % Create output
        if cInd > 1
            peakInds = ind(peakLoc(1:cInd-1));
            peakMags = peakMag(1:cInd-1);
        else
            peakInds = [];
            peakMags = [];
        end % end if
    else % This is a monotone function where an endpoint is the only peak
        [peakMags,xInd] = max(x);
        if includeEndpoints && peakMags > minMag + sel
            peakInds = ind(xInd);
        else
            peakMags = [];
            peakInds = [];
        end % end if
    end % end if
    % Apply threshold value.  Since always finding maxima it will always be
    %   larger than the thresh.
    if ~isempty(thresh)
        m = peakMags>thresh;
        peakInds = peakInds(m);
        peakMags = peakMags(m);
    end % end if
    if interpolate && ~isempty(peakMags)
        middleMask = (peakInds > 1) & (peakInds < len0);
        noEnds = peakInds(middleMask);
        magDiff = x0(noEnds + 1) - x0(noEnds - 1);
        magSum = x0(noEnds - 1) + x0(noEnds + 1)  - 2 * x0(noEnds);
        magRatio = magDiff ./ magSum;
        peakInds(middleMask) = peakInds(middleMask) - magRatio/2;
        peakMags(middleMask) = peakMags(middleMask) - magRatio .* magDiff/8;
    end % end if
    % Rotate data if needed
    if flipData
        peakMags = peakMags.';
        peakInds = peakInds.';
    end % end if
    % Change sign of data if was finding minima
    if extrema < 0
        peakMags = -peakMags;
        x0 = -x0;
    end % end if
    % Plot if no output desired
    if nargout == 0
        if isempty(peakInds)
            disp('No significant peaks found')
        else
            figure;
            plot(1:len0,x0,'.-',peakInds,peakMags,'ro','linewidth',2);
        end % end if
    else
        varargout = {peakInds,peakMags};
    end % end if
    
end % end function
