function show_fft(data, fs, opt)
% This function was written by :
%                               Mark Campbell
%                               Carleton University
% opt.start: The first frame to be included in the time series 
% opt.stop:  The last frame to be included in the time series.
% fs:    The data framerate
if ~exist('opt', 'var')
    opt.start= 1;
    opt.stop= size(data, 2);
end % end if
    figure;
    dr = real(data(:,opt.start:opt.stop));
    tax = linspace(0,fs,size(dr,2)); 
    semilogy(tax,abs(fft(sum(dr, 1))));
    xlim([0,fs/2]);
end % end function