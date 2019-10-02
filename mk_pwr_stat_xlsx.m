function mk_pwr_stat_xlsx(file)

% -------------------------------------------------------------------------
% DESCRIPTION:
%
%   RES= mk_pwr_stat_xlsx(file)
%
%   Write time and frequency domain power statistics to excel file, as well
%   as CLI.
%
%   TODO:
%       parameterize function for different window and CLI calc sizes?
%
% -------------------------------------------------------------------------
% PARAMETERS:
% 
%   file:
%       raw eMotiv csv recording file.
%
% -------------------------------------------------------------------------   
% RETURNS:
% 
%   Does not return. Writes the files to a new directory file name.
%
% -------------------------------------------------------------------------   
% AUTHOR:
%   Mark Campbell
%   Neurovine
%   markacampbell@cmail.carleton.ca
%   01.Oct.2019
% -------------------------------------------------------------------------

starting_dir= pwd;
D= emotiv_load_data(file);
fn= fieldnames(D);
header = {'Window number','Time', 'Delta', 'Theta', 'Alpha', 'Beta Lo', 'Beta Hi', 'Gamma', 'Total Power: Time Domain', 'Total Power: Frequency Domain', 'CLI', 'CLI used', '10s Windowed CLI'};
suffix= '_pwr_stats_CLI';
new_dir_name= horzcat(D.name, suffix);
mkdir(new_dir_name);
cd(new_dir_name);

for i= 3:numel(fn)
    fs= 128;
    window_length= 128;
    offset= 64;
    opt.plot= 0;
    nfft= 128;
    dd= D.(fn{i}).data;
    % Total power
    pxx= welch_psd(dd, fs, window_length, offset, nfft, opt); % for entire time series
%     total_power_pxx= sum(pxx);
    % Total Power: Time Domain
    wnd_pwr_time= get_window_power(dd, window_length, offset, fs, 'time'); % var power per window
    % Total Power: Frequency Domain
    wnd_pwr_freq= get_window_power(dd, window_length, offset, fs, 'freq'); % var power per window
    % Window number
    window_idx= (1:length(wnd_pwr_time))';
    % Time
    time= window_idx* 0.5;
    % Delta
    delta= wnd_pwr_freq(:, 1);
    % Theta
    theta= wnd_pwr_freq(:, 2);
    % Alpha
    alpha= wnd_pwr_freq(:, 3);
    % Beta Low
    betaL= wnd_pwr_freq(:, 4);
    % Beta High
    betaH= wnd_pwr_freq(:, 5);
    % Gamma
    gamma= wnd_pwr_freq(:, 6);
    
    % 10s Windowed CLI and CLI used
    n= 20;
    offset= 1;
    thresh= 1000;
    out= calc_windowed_cli(wnd_pwr_freq, wnd_pwr_time, n, offset, thresh);
    cli_wnd= out.cli_wnd;
    cli_use= out.cli_use;
    
    % CLI
    cli= theta./ alpha;
    
    % Make output cell
    RES = [window_idx, time, delta, theta, alpha, betaL, betaH, gamma, wnd_pwr_time, sum(wnd_pwr_freq,2), cli, cli_use, cli_wnd];
    
    % Save
    elec= fn{i};
    save_name= horzcat(D.name, '_', elec, suffix, '.xlsx');
    xlswrite(save_name, RES, D.(fn{i}).name, 'A2' );
    xlswrite(save_name, header, D.(fn{i}).name, 'A1' );

end % end for

cd(starting_dir);

end % end function