function D= load_HamburgBrain_data(pig)
% -------------------------------------------------------------------------
% DESCRIPTION:
%
%   D= load_HamburgBrain_data(pig)
%
% Each sequence in the data struct is an experimental trial. Each sequence
% has the following fields:
% -------------------------------------------------------------------------
% name: Formatted file name used as title of figures for this sequence
% pig: Pig name
% imgr: Reconstructed images using eit.fdata and imdl for this pig
% eit: EIT part of data struct
%     data: EIT data output from eidors_readdata()
%     fs: EIT data framerate
%     sync1: First syncronization spike in EIT data
%     sync2: Second syncronization spike in EIT data
%     apn: Using the LabChart file, the time relative to the beginning of the first EIT synchronization spike that the apnoea flush occured in the arterial pressure data
%     inj: Same as above, but for the injection flush
%     vnt: Same as above, but for the end of apnoea flush
%     fdata: 9 Hz lowpass-filtered EIT data
% perf: Arterial pressure part of data struct
%     data: Perfusion data
%     tickrate: Arterial pressure framerate
%     apn: Apnoea flush spike in arterial pressure data
%     inj: Injection flush spike in arterial pressure data
%     vnt: End of apnoea flush spike in arterial pressure data
% fdata now calculated in align sequence to trim filter edge effect.
% 
% To change the lowpass filter stopband, go to lowpass_iir function.
% 
% Author:
%   Mark Campbell
%   Carleton University
%   markacampbell@cmail.carleton.ca
% eitfs= 47.68;
maxsz= 0.2; maxh= 2; imgsize= [64 64];
[fmdl, imdl]= mk_pighead_fmdl(maxsz, maxh, imgsize);
% lpass= 9;

% 8.2
if pig== "8.2"
    
    eit_files= {    'EIT_8.2_nativ_1.eit',                  'EIT_8.2_nativ_2.eit',... 
                    'EIT_8.2_30_Min_nach_Embolisation.eit', 'EIT_8.2_rechts_embolisiert.eit'};    
    perf_files= {   'EIT_8.2_nativ_1.mat',                  'EIT_8.2_nativ_2.mat',... 
                    'EIT_8.2_30_Min_nach_Embolisation.mat', 'EIT_8.2_rechts_embolisiert.mat'};
    sync1= [681,    98,     941,    1027];
    sync2= [6370,   4990,   4633,   5264];
    e_apn= [42.5,   14.5,   -9.6,   15.7];
    e_inj= [50.8,   31.1,   12.7,   26.1];
    e_vnt= [98.2,   83.9,   58,     78.8];
    p_apn= [53517,  14282,  7622,   30161];
    p_inj= [61910,  30896,  30011,  40589];
    p_vnt= [109097, 83697,  75311,  93234];

elseif pig== "9.2"
    
    eit_files= {    'EIT_9.2_Nativ_1.eit',                  'EIT_9.2_Nativ_2.eit',...
                    'EIT_nach_Embolisation_1_9.2.eit',      'EIT_nach_Perfusionsminderun_2_9.2.eit'};    
    perf_files= {   'EIT_9.2_Nativ_1.mat',                  'EIT_9.2_Nativ_2.mat',... 
                    'EIT_nach_Embolisation_1_9.2.mat',      'EIT_nach_Perfusionsminderun_2_9.2.mat'};
    sync1= [246,    291,    979,    239];
    sync2= [5239,   5948,   5956,   6382];
    e_apn= [12.3,   12.1,   14.7,   33.7];
    e_inj= [27,     19.9,   20.7,   41];
    e_vnt= [82.1,   98.7,   83.9,   112.1];
    p_apn= [132012, 140143, 173558, 236714];
    p_inj= [146787, 148520, 179538, 243968];
    p_vnt= [201855, 226761, 242739, 315069];
    
elseif pig== "10.2"
    
    eit_files= {    'EIT_nativ_10.2._schleuse.eit',                             'EIT_nativ_2_zvk_10.2.eit',... 
                    'EIT_nach_Perfusionsminderung_Sequenz_3_Schleuse_10.2.eit', 'EIT_nach_Perfusionsminderung_Sequenz_4_ZVK.eit',...
                    'EIT_10.2._nach_6h_Perfusionsmin_Schleuse Sequ5.eit',       'EIT_nach_6h_Perfusionsminderung_10.02 ZVK Sequ6.eit'};    
    perf_files= {   'EIT_nativ_10.2._schleuse.mat',                             'EIT_nativ_2_zvk_10.2.mat',... 
                    'EIT_nach_Perfusionsminderung_Sequenz_3_Schleuse_10.2.mat', 'EIT_nach_Perfusionsminderung_Sequenz_4_ZVK.mat',...
                    'EIT_10.2._nach_6h_Perfusionsmin_Schleuse Sequ5.mat',       'EIT_nach_6h_Perfusionsminderung_10.02 ZVK Sequ6.mat'};    
    sync1= [399,    262,    124,    149,    177,    191];
    sync2= [5232,   6729,   4866,   5390,   5558,   4995];
    e_apn= [7.4,    9.3,    12.8,   11.9,   14,     14.5];
    e_inj= [13.4,   15.6,   18.6,   21.9,   20.1,   20.4];
    e_vnt= [79.1,   85.5,   73,     88.5,   87.7,   81.2];
    p_apn= [13720,  12336,  17460,  13854,  30897,  15381];
    p_inj= [19708,  18832,  23357,  23835,  37773,  21370];
    p_vnt= [85388,  88633,  77865,  90452,  104515, 82125];
    
elseif pig== "11.2"
    
    eit_files= {    'EIT_11.2_Nativ_1_Schleuse.eit',                'EIT_11.2_Nativ_2_ZVK.eit',... 
                    'EIT_Sequenz_3_Schleuse.eit',                   'EIT_Seqeunz4_zvk.eit',...
                    'EIT_11.2_4h_nach_Stroke_Schleuse_Sequ5.eit',   'EIT_4h_nach_Stroke_ZVk_Sequ6.eit'};
    perf_files= {   'EIT_11.2_Nativ_1_Schleuse.mat',                'EIT_11.2_Nativ_2_ZVK.mat',...
                    'EIT_Sequenz_3_Schleuse.mat',                   'EIT_Seqeunz4_zvk.mat',...
                    'EIT_11.2_4h_nach_Stroke_Schleuse_Sequ5.mat',   'EIT_4h_nach_Stroke_ZVk_Sequ6.mat'};  
    sync1= [1225,   211,    261,    55,     495,    169];
    sync2= [6369,   5694,   7008,   5814,   5322,   6666];
    e_apn= [-10.2,  10.7,   25.1,   22.3,   8.1,    12.1];
    e_inj= [11.5,   15.3,   32.2,   66,     18.2,   16.2];
    e_vnt= [89.1,   91.9,   113.7,  104.4,  80.4,   85.8];
    p_apn= [11690,  10714,  36930,  31092,  16000,  19419];
    p_inj= [33376,  15269,  43994,  74955,  26149,  23516];
    p_vnt= [111004, 91873,  125489, 113239, 88339,  93187];
    
elseif pig== "12.2"
    
    eit_files= {    'EIT_12.2. Seqeunz 1 nativ Schleuse.eit',           'EIT_12.2_NATIV_ZVK_Sequ_2.eit',... 
                    'EIT_direkt_nach_Stroke_12.2._Sequ_3_Schleuse.eit', 'EIT_12.2_Sequ_4_nach_Stroke_ZVK.eit',...
                    'EIT_12.2._3,5_nach_Stroke_Sequ_5_Schleuse.eit',    'EIT_12.2._Sequ_6_3,5h_nach_Stroke_ZVK.eit'};
    perf_files= {   'EIT 12.2 sequence 3.mat',                          'EIT 12.2 sequence 4.mat',...
                    'EIT 12.2 sequence 5.mat',                          'EIT 12.2 sequence 6.mat',...
                    'EIT 12.2 sequence 7.mat',                          'EIT 12.2 sequence 8.mat'};     
    sync1= [251,    156,    126,    518,    132,    245];
    sync2= [5754,   7179,   7600,   7420,   6586,   7715];
    e_apn= [8.5,    14.5,   12.7,   13,     17.4,   12.6];
    e_inj= [15,     23.4,   50,     36.2,   37.5,   19.9];
    e_vnt= [94.7,   117,    128,    126.6,  107.6,  117];
    p_apn= [10415,  15490,  15807,  22384,  24336,  12550];
    p_inj= [16873,  24463,  53138,  45609,  44446,  19840];
    p_vnt= [96616,  116350, 131095, 136047, 114530, 116881];

end % end if

D= load_data(eit_files, perf_files, sync1, sync2, e_apn, e_inj, e_vnt, p_apn, p_inj, p_vnt, imdl, pig);

end % end function


function D= load_data(eit_files, perf_files, sync1, sync2, e_apn, e_inj, e_vnt, p_apn, p_inj, p_vnt, imdl, pig)

% Field names
fn= {'seq1', 'seq2', 'seq3', 'seq4', 'seq5', 'seq6'};

for i= 1:length(eit_files)
    
    % load data
    [vv,aux]= eidors_readdata(eit_files{i});
    eitfs= median(1e6/median(diff(aux.t_rel)));
    
    % assign struct fields
    D.(fn{i}).eit.data= vv;
    D.(fn{i}).eit.elec_impedance= aux.elec_impedance/1e3;
    D.(fn{i}).name= char(remove_underscores(eit_files{i}));
    D.(fn{i}).pig= pig;
    D.(fn{i}).eit.fs= eitfs;
    D.(fn{i}).eit.sync1= sync1(i);
    D.(fn{i}).eit.sync2= sync2(i);
    D.(fn{i}).eit.apn= D.(fn{i}).eit.sync1 + round(e_apn(i)* eitfs);
    D.(fn{i}).eit.inj= D.(fn{i}).eit.sync1 + round(e_inj(i)* eitfs);
    D.(fn{i}).eit.vnt= D.(fn{i}).eit.sync1 + round(e_vnt(i)* eitfs);
    D.(fn{i}).perf= load( perf_files{i} );
    D.(fn{i}).perf.apn= p_apn(i);
    D.(fn{i}).perf.inj= p_inj(i);
    D.(fn{i}).perf.vnt= p_vnt(i);
    
    % allign eit and perfusion sequences, trim syncronization spikes
    D.(fn{i})= allign_eit_and_perf(D.(fn{i}));
    
    % adjust recontruction matrix for noisy electrodes
    imdl_comp= compensate_bad_elec(D.(fn{i}).eit.data, imdl);
    
    % use only good measurements
    msel= imdl.fwd_model.meas_select;
    mm = find(msel);
    use_data= D.(fn{i}).eit.data(mm, :);
    
    % filter eit data
    lowpass_cutoff= 8;
    D.(fn{i}).eit.fdata= lowpass_iir(use_data, lowpass_cutoff, D.(fn{i}).eit.fs);
    
    % compensate for noisy electrodes and calculate reconstructed images
    D.(fn{i}).imgr= get_imgr( D.(fn{i}).eit.fdata, D.(fn{i}).eit.inj, imdl_comp );
    
end % end for

end % end function


function seq= allign_eit_and_perf(seq)
    % Find optimal alignment between eit and perfusion data
    % Annotate peaks and troughs of perfusion signal and transfer annotations
    % to EIT data
    valDiffThresh= 400;

    estart= seq.eit.sync1+ round(seq.eit.fs); % desired starting frame is 1 s after first sync
    eend= seq.eit.sync2- round(seq.eit.fs); % desired end frame is 1 s before last sync
    pstart= round((estart/ seq.eit.fs)* seq.perf.tickrate);
    pend= round((eend/ seq.eit.fs)* seq.perf.tickrate);

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
        ecut= estart+ round( shift*  seq.eit.fs);
        pcut= pstart;   
    elseif shift< 0 % Perfusion sequence is ahead of EIT sequence
        ecut= estart;
        pcut= pstart+ abs(round( shift*  seq.perf.tickrate));
    end % end if

    seq.eit.data=   seq.eit.data(:, ecut+1: eend);

    seq.eit.sync1=  seq.eit.sync1- ecut;
    seq.eit.sync2=  seq.eit.sync2- ecut;
    seq.eit.apn=    seq.eit.apn- ecut;
    seq.eit.inj=    seq.eit.inj- ecut;
    seq.eit.vnt=    seq.eit.vnt- ecut;

    seq.perf.data=  movmean(seq.perf.data(:, pcut+1: pend+ pcut), 10);
    seq.perf.apn=   seq.perf.apn- pcut;
    seq.perf.inj=   seq.perf.inj- pcut;
    seq.perf.vnt=   seq.perf.vnt- pcut;

    % Find peaks and valleys in perfusion Signal
    seq.perf.vals=  peakfinder(seq.perf.data, 8, [], -1);

    % Remove false valleys that are too close to true valleys
    for i= 1:length(seq.perf.vals)- 1
        if seq.perf.vals(i+1)- seq.perf.vals(i) < valDiffThresh 
            seq.perf.vals(i+1)= seq.perf.vals(i);
        end % end if
    end % end for

    % Remove valleys within the synchronization flush
    invalids= [seq.perf.apn- 100, seq.perf.inj- 100, seq.perf.vnt- 100];
    for invalid= [seq.perf.apn, seq.perf.inj, seq.perf.vnt]
        rm= find( (seq.perf.vals> invalid) + (seq.perf.vals< invalid+ 1100)==2 );
        if ~isempty(rm)
            seq.perf.vals(rm)= 0;
        end % end if
    end % end if

    seq.perf.vals=  seq.perf.vals(seq.perf.vals~= 0);
    seq.perf.vals=  unique(seq.perf.vals);
    seq.perf.peaks= zeros(1, length(seq.perf.vals)- 1);

    for i= 1:length(seq.perf.peaks)
        start= seq.perf.vals(i);
        stop= seq.perf.vals(i+1);
        loc= find( seq.perf.data== max(seq.perf.data(start: stop)));
        seq.perf.peaks(i)= loc(find((loc> start) + (loc< stop)==2, 1));
    end % end for

    % Remove peaks within the synchronization flush
    for invalid= invalids
        rm= find( (seq.perf.peaks> invalid) + (seq.perf.peaks< invalid+ 1100)==2 );
        if ~isempty(rm)
            seq.perf.peaks(rm)= 0;
        end % end if
    end % end if

    seq.perf.peaks= seq.perf.peaks(seq.perf.peaks~= 0);

    % Transfer these annotations to the EIT data
    seq.eit.peaks=  round( (seq.perf.peaks./ seq.perf.tickrate).*seq.eit.fs );
    seq.eit.vals=   round( (seq.perf.vals./ seq.perf.tickrate).*seq.eit.fs );

end % end function

function imgr= get_imgr(fdata, inj, imdl)
    
    vv= real(fdata);
    vh= vv(:, inj); % reference frame is injection frame
    imgr= inv_solve(imdl, vh, vv);

    % Set colour limits
    clim= max(imgr.elem_data(:));
    imgr.calc_colours.ref_level= 0;
    imgr.calc_colours.lim= clim;

end % end function

function varargout = peakfinder(x0, sel, thresh, extrema, includeEndpoints, interpolate)
% https://www.mathworks.com/matlabcentral/fileexchange/25500-peakfinder-x0-sel-thresh-extrema-includeendpoints-interpolate
%PEAKFINDER Noise tolerant fast peak finding algorithm
%   INPUTS:
%       x0 - A real vector from the maxima will be found (required)
%       sel - The amount above surrounding data for a peak to be,
%           identified (default = (max(x0)-min(x0))/4). Larger values mean
%           the algorithm is more selective in finding peaks.
%       thresh - A threshold value which peaks must be larger than to be
%           maxima or smaller than to be minima.
%       extrema - 1 if maxima are desired, -1 if minima are desired
%           (default = maxima, 1)
%       includeEndpoints - If true the endpoints will be included as
%           possible extrema otherwise they will not be included
%           (default = true)
%       interpolate - If true quadratic interpolation will be performed
%           around each extrema to estimate the magnitude and the
%           position of the peak in terms of fractional indicies. Note that
%           unlike the rest of this function interpolation assumes the
%           input is equally spaced. To recover the x_values of the input
%           rather than the fractional indicies you can do:
%           peakX = x0 + (peakLoc - 1) * dx
%           where x0 is the first x value and dx is the spacing of the
%           vector. Output peakMag to recover interpolated magnitudes.
%           See example 2 for more information.
%           (default = false)
%
%   OUTPUTS:
%       peakLoc - The indicies of the identified peaks in x0
%       peakMag - The magnitude of the identified peaks
%
%   [peakLoc] = peakfinder(x0) returns the indicies of local maxima that
%       are at least 1/4 the range of the data above surrounding data.
%
%   [peakLoc] = peakfinder(x0,sel) returns the indicies of local maxima
%       that are at least sel above surrounding data.
%
%   [peakLoc] = peakfinder(x0,sel,thresh) returns the indicies of local
%       maxima that are at least sel above surrounding data and larger
%       (smaller) than thresh if you are finding maxima (minima).
%
%   [peakLoc] = peakfinder(x0,sel,thresh,extrema) returns the maxima of the
%       data if extrema > 0 and the minima of the data if extrema < 0
%
%   [peakLoc] = peakfinder(x0,sel,thresh,extrema, includeEndpoints)
%       returns the endpoints as possible extrema if includeEndpoints is
%       considered true in a boolean sense
%
%   [peakLoc, peakMag] = peakfinder(x0,sel,thresh,extrema,interpolate)
%       returns the results of results of quadratic interpolate around each
%       extrema if interpolate is considered to be true in a boolean sense
%
%   [peakLoc, peakMag] = peakfinder(x0,...) returns the indicies of the
%       local maxima as well as the magnitudes of those maxima
%
%   If called with no output the identified maxima will be plotted along
%       with the input data.
%
%   Note: If repeated values are found the first is identified as the peak

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
