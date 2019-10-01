function ensemble= get_ensembles(seq, ref, opt)
% -------------------------------------------------------------------------
% DESCRIPTION:
%
%   ensemble= get_ensembles(seq, ref, opt)
%
%   Returns an stuct with fields eit and perfusion, imgr, imgr_time (or any
%   combination of the four). Each struct is an n x m matrix using either
%   the EIT or arterial pressure signals, where n is the number of cardiac
%   cycles and m is the number of frames per cardiac cycle. The user can
%   specify which landmark (peaks (1) or valleys (2)) of the perfusion
%   signal is used as the reference for segmenting cardiac cycles. Use
%   opt.pixels to specifcy pixel index in reconstructed image use
%   opt.channels to specify measurement pair in EIT data.
% -------------------------------------------------------------------------
% PARAMETERS:
%
%   seq:        
%       A swine sequence struct from load_HamburgBrain_data.
%   ref:        
%       which perfusion reference landmark to use. ref= 1(peaks) or
%       ref = 2 (valleys) opt: options
%   opt:        
%       An options struct. Possible fields include:
%       start:
%       stop:
%       pixels:
%       channels:
%       sel:
%           A 1x3 cell containing 0's or 1 's specifying which ensembles to
%           return. Can specify any comination of the ensembles listed
%           here:
%           eit_ensemble:       [1 0 0]
%           perf_ensemble:      [0 1 0]
%           imgr_ensemble:      [0 0 1]
% -------------------------------------------------------------------------
% RETURNS:
%
%   ensemble:   
%       A struct containing unprocessed windowed data
%       with window start referenced to cardiac signal. Depending on
%       opt.sel, this function will return:
%
%       eit_ensemble:
%
%
%       perf_ensemble:
%
%
%       imgr_ensemble:
%
%           imgr_ensemble= cell(n_pixels, n_windows,
%           ensemble_length);
%
%           Can average this ensemble such that each frame is the average
%           of a signle cardiac cycle, and length of time series is
%           n_windows, by taking the mean across the second dimension. Each
%           subsequent frame will be the average of a subsequent cardiac
%           cycle.
%
%           Can average this ensemble such that each frame is the
%           average of one of the frames in every cardiac cycle, and the
%           length of the time series is ensemble_length, by taking the
%           mean across the third dimension.
%
% -------------------------------------------------------------------------
% AUTHOR:
%   Mark Campbell
%   Carleton University
%   markacampbell@cmail.carleton.ca
% -------------------------------------------------------------------------

% Check options
if nargin== 2 && ~exist('opt', 'var')
    opt.start= 1;
    opt.stop= size(seq.eit.fdata, 2);
    opt.pixels= 1: size(seq.eit.fdata, 1);
    opt.sel= [1 1 1];
end % end if

if isfield(opt, 'start')
    start= opt.start;
else
    start= 1;
end % end if

if isfield(opt, 'stop')
    stop= opt.stop;
else
    stop= size(seq.eit.fdata, 2);
end % end if

if isfield(opt, 'pixels')
    use_pixels= opt.pixels;
else
    use_pixels= 1: size(seq.imgr.elem_data, 1);
end % end if

if isfield(opt, 'channels')
    use_channels= opt.pixels;
else
    use_channels= 1: size(seq.eit.fdata, 1);
end % end if

if isfield(opt, 'sel')
    get_eit=        opt.sel(1);
    get_perf=       opt.sel(2);
    get_imgs=       opt.sel(3);
else
    get_eit=        1;
    get_perf=       1;
    get_imgs=       1;
end % end if

% check ensemble reference
if ~exist('ref', 'var')
    ref= 2; 'use valleys';
elseif ref~= 1 && ref~= 2
    disp("Ref must be either 1 or 2");
end % end if

if ref== 1
    eit_intervals= diff(seq.eit.peaks);
    perfusion_intervals= diff(seq.perf.peaks);
    eit_idx= seq.eit.peaks;
    perf_idx= seq.perf.peaks;
else
    eit_intervals= diff(seq.eit.vals); 
    eit_idx= seq.eit.vals;
    perfusion_intervals= diff(seq.perf.vals);
    perf_idx= seq.perf.vals;
end % end if

% exclude large intervals that may be faulty annotations
eit_ensemble_length= round( 1.0* mean(eit_intervals(eit_intervals< seq.eit.fs* 1.25)) ); 
perf_ensemble_length= round( 1.0* mean(perfusion_intervals(perfusion_intervals< seq.perf.tickrate* 1.25)) );

% calculate EIT indices
% Put intervals in frame with eit data and trim as neccessary. Want
% perfusion index and data to be in frame and start with index of 1.
% Let's assume first index of perf ref is now 3. It is offset from the
% desired start of the eit sequence by 3- 1= 2
n_samples= stop- start+ 1;
last_idx= stop- (eit_ensemble_length- 1);

% get rid of indices that are too large
eit_idx= eit_idx( (eit_idx <= last_idx) );
perf_idx= perf_idx(1:length(eit_idx));

% get rid of indices that are too small
eit_idx= eit_idx( (eit_idx- (start-1)) > 0);
perf_start= length(perf_idx)- length(eit_idx) + 1;
perf_idx= perf_idx(perf_start:end);

% set first index to 1.
offset= eit_idx(1)- start;
eit_idx= eit_idx- (eit_idx(1)- 1);

% EIT signal ensembles
if get_eit== 1
    eit_data= seq.eit.fdata(use_channels, start+ offset: stop);
    ensemble.eit_ensemble= mk_ensemble_mtx(eit_data, eit_idx, eit_ensemble_length);
end % end if

% Perfusion signal ensembles
if get_perf== 1
    ensemble.perf_ensemble= mk_ensemble_mtx(seq.perf.data, perf_idx, perf_ensemble_length);
end % end if

% Reconstructed image ensembles.
if get_imgs== 1
    imgr_data= seq.imgr.elem_data(use_pixels, start+ offset: stop);
    ensemble.imgr_ensemble= mk_ensemble_mtx(imgr_data, eit_idx, eit_ensemble_length);
end % end if

end % end function

function ensemble_mtx= mk_ensemble_mtx(data, wnd_idx, ensemble_length)

% idx holds values for the start of each frame. Idx starts at 1, but check
% that the last index + ensemble length is not greater than the length of
% our data

% returns ensemble_mtx: n_pixel x n_window x ensemble_length

    n_samples= size(data, 2);
    n_pixels= size(data, 1);
    check_idx= 1;
    
    while check_idx== 1
        if wnd_idx(end)+ (ensemble_length-1) > n_samples
            wnd_idx= wnd_idx(1:end- 1);
        else
            check_idx= 0;
        end % end if
    end % end while

    n_windows= length(wnd_idx);
    ensemble_mtx= zeros(n_pixels, n_windows, ensemble_length);
    
    for i= 1:n_windows
        idx= wnd_idx(i);
        ensemble_mtx(:, i, :)= data(:, idx: idx+ ensemble_length-1);
    end % end for
     
end % end function
