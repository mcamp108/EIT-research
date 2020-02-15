function data_out= lowpass_iir(data_in, loPassCutoff, fs)

% Do lowpass iir filter on EIT data
% 
% Author:
%   Mark Campbell
%   Carleton University
%   markacampbell@cmail.carleton.ca

data_in = (real(data_in));
n_meas = size(data_in, 1);

pad_sz = round(fs*2);
pad1 = repmat(data_in(1,:), pad_sz, 1);
pad2 = repmat(data_in(end,:), pad_sz, 1);

data_in = [pad1; data_in; pad2];
data_out = lowpass(data_in, loPassCutoff, fs);
data_out = data_out(pad_sz+1: pad_sz+n_meas, :);

end % end function
