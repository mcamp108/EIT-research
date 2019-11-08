function inspect_eit_elec_and_data(eit_file, imdl)
% -------------------------------------------------------------------------
% Description:
%   inspect_eit_elec_and_data(eit_file, imdl)
% -------------------------------------------------------------------------
% Parameters:
%   eit_file:
%       Target data file ending in .eit file extension.
%   imdl:
%       EIDORS inverse model structure
% ------------------------------------------------------------------------- 
% Returns:
%   IQ plot
%   U-shape plot
%   Contact impedance plot
%   Frequency spectrum plot
% ------------------------------------------------------------------------- 
% Author:
%   Mark Campbell
%   Carleton University
%   markacampbell@cmail.carleton.ca
%   27.Sep.2019
% -------------------------------------------------------------------------

[data,auxdata]= eidors_readdata(eit_file);

% find bad electrodes
impedanceFactor =  2.048 / (2^12 * 0.173 * 0.003) / 2^15; % = 0.9633 / 2^15;
elec_impedance= abs(auxdata.elec_impedance* impedanceFactor);
ei= median(elec_impedance, 2);
bad_elecs= find(ei> 400);

fs = 1e6./median(diff(auxdata.t_rel)); %framerate is median dif of time points/ 1000000 (convert to s)
msel= imdl.fwd_model.meas_select;

subplot(325);
    vk = 1e3* mean(real(data), 2);
    by_stim_pair= reshape(vk, 32, 32);
    av_meas= (mean(by_stim_pair, 2));
    plot(av_meas, 'b');
    hold on;
    scatter(bad_elecs, av_meas(bad_elecs), 'r', 'filled');
%     eb= errorbar(1:32, av_meas, max(by_stim_pair, [], 2)-av_meas, min(by_stim_pair, [], 2)- av_meas);
%     set(eb,'Color',[0, 0, 0],'LineStyle','none');
    hold off;
    xlim([1 32]);
    title("Average Measurement for Measurement Pair");
    legend('meas', 'bad elecs');

subplot(321);
    % all meas
    plot(1e3* data(:, 1:100), 'k+'); 
    hold on;
    % used meas
    plot(1e3* data(msel, 1:100), 'b+');
    kk=meas_icov_rm_elecs(imdl, bad_elecs);
    ee = find(diag(kk)~=1);
    mm = find(msel);
    ee = mm(ee);
    plot(1e3*data(ee,1:100), 'r+');
    hold off; 
    box off
    title 'IQ plot';
  
subplot(322);
    vk = 1e3* mean(real(data), 2);
    plot(vk,'k');    % all meas
    hold on;
    vk(~msel,:) = NaN;
    plot(vk,'b');   % selected meas
    vn = NaN*vk; vn(ee,:) = vk(ee,:);
    plot(vn,'ro');
    hold off;
    box off; 
    xlim([0,400]);
    title 'U shapes';

subplot(323);
    bar(1:32, ei);
    hold on;
    eb= errorbar(1:32, ei, max(elec_impedance, [], 2)-ei, min(elec_impedance, [], 2)- ei);
    set(eb,'Color',[0, 0, 0],'LineStyle','none');
    hold off;
    box off; 
    xlim([1, 32]);
    title 'Elec Z';

subplot(324);
    ll = size(data, 2);    
    ft = fft(detrend(data.').', [], 2);
    fax = linspace(0, fs, ll+1); fax(end)=[];
    semilogy(fax, mean(abs(ft)));
    ylim([1e-4, 1]);
    box off; 
    xlim([0, fs/2]);
    title 'Spectrum (Hz)';
    
end % end function