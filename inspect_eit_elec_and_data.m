function inspect_eit_elec_and_data(seq, imdl)

% Description:
%   
%
% Parameters:
%   
%   
% Returns:
%   IQ plot
%   U-shape plot
%   Contact impedance plot
%   Frequency spectrum plot
%   
% Author:
%   Mark Campbell
%   Carleton University
%   markacampbell@cmail.carleton.ca
%   27.Sep.2019

vv= seq.eit.data;
FR= seq.eit.fs;
msel= imdl.fwd_model.meas_select;

subplot(325);
    vk = 1e3* mean(real(vv), 2);
    by_stim_pair= reshape(vk, 32, 32);
    av_meas= (mean(by_stim_pair, 2))';
    bad_elecs= find( abs( av_meas- mean(av_meas) )> std(av_meas));
    plot(av_meas, 'b');
    hold on;
    scatter(bad_elecs, av_meas(bad_elecs), 'r', 'filled');
    xlim([0 33]);
    title("Average Measurement for Measurement Pair");
    legend('meas', 'bad elecs');

subplot(321);
    % all meas
    plot(1e3* vv(:, 1:100), 'k+'); 
    hold on;
    % used meas
    plot(1e3* vv(msel, 1:100), 'b+');
    kk=meas_icov_rm_elecs(imdl, bad_elecs);
    ee = find(diag(kk)~=1);
    mm = find(msel);
    ee = mm(ee);
    plot(1e3*vv(ee,1:100), 'r+');
    hold off; 
    box off
    title 'IQ plot';
  
subplot(322);
    vk = 1e3* mean(real(vv), 2);
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
    ei= abs(seq.eit.elec_impedance);
    aei = median(ei, 2);
    bar(1:32, aei);
    hold on;
    eb= errorbar(1:32, aei, max(ei, [], 2)-aei, min(ei, [], 2)- aei);
    set(eb,'Color',[0, 0, 0],'LineStyle','none');
    hold off;
    box off; 
    xlim([1, 32]);
    title 'Elec Z';

subplot(324);
    ll = size(vv, 2);    
    ft = fft(detrend(vv.').', [], 2);
    fax = linspace(0, FR, ll+1); fax(end)=[];
    semilogy(fax, mean(abs(ft)));
    ylim([1e-4, 1]);
    box off; 
    xlim([0, 12]);
    title 'Spectrum (Hz)';
    
end % end function

function bad_elecs= find_bad_elecs(vv)


    
end % end function