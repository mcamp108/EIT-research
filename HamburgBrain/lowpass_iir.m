function data_out= lowpass_iir(data_in, loPassCutoff, fs)

% Demean and do lowpass iir filter on data_in
% This assumes that EIT data has been inputted with rows= time axis,
% columns= measurement frames (as output from read_eidors)
% 
% Author:
%   Mark Campbell
%   Carleton University
%   markacampbell@cmail.carleton.ca

% data_in= data_in- mean(data_in, 2); % remove dc
data_in= (real(data_in))';

% prefix and append data to prevent filter edge effect
% https://www.mathworks.com/matlabcentral/answers/161223-how-to-remove-transient-effect-in-the-beginning-of-the-filtered-signal
R= 0.1;
Nr= 50;
N= size(data_in, 1);
NR= min(round(N* R), Nr); % At most 50 points

for i=1:size(data_in, 2)
    x1(:, i)= 2* data_in(1, i)- flipud(data_in(2: NR+ 1, i));  % maintain continuity in level and slope
    x2(:, i)= 2* data_in(end, i)- flipud(data_in(end- NR:end- 1, i));
end % end for

data_in= [x1; data_in; x2];

% % from eeglab toolbox
% minfac = 3;    % this many (lo)cutoff-freq cycles in filter
% lopass_filt_order= max(15, minfac* fix(fs/loPassCutoff));
% d= designfilt('lowpassiir', 'SampleRate', fs, 'PassbandRipple', 1, 'FilterOrder',... 
%                lopass_filt_order, 'PassbandFrequency', loPassCutoff, 'DesignMethod', 'cheby1');
% d = designfilt('lowpassfir', 'PassbandFrequency', 8, ...
%                'StopbandFrequency', 10, 'PassbandRipple', 1, ...
%                'StopbandAttenuation', 60, 'SampleRate', 47, ...
%                'DesignMethod', 'equiripple');

d= designfilt('lowpassfir', 'SampleRate', 47.6826, 'PassbandRipple', 1, ...
              'PassbandFrequency', 5, 'StopbandFrequency', 6, 'StopbandAttenuation', 60,...
              'DesignMethod', 'equiripple');
              
filt_data= filter(d, data_in);

data_out= (filt_data(NR+ 1: end- NR, :))';

end % end function
