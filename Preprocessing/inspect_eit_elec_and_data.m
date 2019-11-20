function inspect_eit_elec_and_data(eit_file, imdl, thresh)
% -------------------------------------------------------------------------
% Description:
%   inspect_eit_elec_and_data(eit_file, imdl)
% -------------------------------------------------------------------------
% Parameters:
%   eit_file:
%       Target data file ending in .eit file extension or a cell containing
%       the outpus of eidors_readdata {data, auxdata}.
%   imdl:
%       EIDORS inverse model structure
%   thresh:
%       contact impedance threshold at which to label an electrode poor.
% ------------------------------------------------------------------------- 
% Returns:
%   Total boundary voltage plot
%   IQ plot
%   U-shape plot
%   Frequency spectrum plot
%   Contact impedance plot
%   Average measurement for measureemnt pair plot
% ------------------------------------------------------------------------- 
% Author:
%   Mark Campbell
%   Carleton University
%   markacampbell@cmail.carleton.ca
%   27.Sep.2019
% -------------------------------------------------------------------------

if (nargin== 2) && (~exist(thresh, 'var'))
   thresh= 400; 
end % end if

if iscell(eit_file) % input is output of eidors_read_data
    if length(eit_file)==2
        if isstruct(eit_file{2})
            auxdata= eit_file{2};
            data= eit_file{1};
        else
            auxdata= eit_file{1};
            data= eit_file{2};
        end % end if
    else
        disp("Whoops!");
    end % end if
elseif strcmp(eit_file(end-3:end), '.eit')
    [data,auxdata]= eidors_readdata(eit_file);
else
    disp("Unrecognized input.");
end % end if

fs = 1e6./median(diff(auxdata.t_rel)); %framerate is median dif of time points/ 1000000 (convert to s)
msel= imdl.fwd_model.meas_select;

% find bad electrodes
impedanceFactor =  2.048 / (2^12 * 0.173 * 0.003) / 2^15; % = 0.9633 / 2^15;
elec_impedance= abs(auxdata.elec_impedance* impedanceFactor);
ei= median(elec_impedance, 2);
bad_elecs= find(ei> thresh);

kk=meas_icov_rm_elecs(imdl, bad_elecs);
ee = find(diag(kk)~=1);

% total boundary voltage
subplot(4,2,1:4);
    vv= detrend(real(data));
    xax= (1:size(vv, 2))/ fs;
    plot(xax, sum(vv(msel,:), 1));
    ylabel('Voltage (V)')
    xlabel('Time (s)');
    title('Total Boundary Voltage');

% U shapes 
subplot(426);
    vk = 1e3* mean(vv, 2);
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
% IQ plot
subplot(425);
    % all meas
    plot(1e3* data(:, 1:100), 'k+'); 
    hold on;
    % used meas
    plot(1e3* data(msel, 1:100), 'b+');
    mm = find(msel);
    ee = mm(ee);
    plot(1e3*data(ee,1:100), 'r+');
    hold off; 
    box off
    title 'IQ plot';

% median contact impedance
subplot(428);
    b= bar(ei,'FaceColor','flat');
    b.CData(:,:)= repmat([0 0 1], 32, 1);
    be= find(ei>=thresh);
    for e= 1:length(be)
        b.CData(be(e), :)= [1 0 0];
    end
    hold on;
    eb= errorbar(1:32, ei, max(elec_impedance, [], 2)-ei, min(elec_impedance, [], 2)- ei);
    set(eb,'Color',[0, 0, 0],'LineStyle','none');
    hold off;
    box off; 
    xlim([0, 33]);
    title 'Median Elec Z';

% Fourier Series
subplot(427);
    ll = size(data, 2);    
    ft = fft(detrend(data.').', [], 2);
    fax = linspace(0, fs, ll+1); fax(end)=[];
    semilogy(fax, mean(abs(ft)));
    ylim([1e-4, 1]);
    box off; 
    xlim([0, fs/2]);
    title 'Spectrum (Hz)';

% Average Measurement for Measurement Pair
% subplot(427);
%     vk = 1e3* mean(vv, 2);
%     by_stim_pair= reshape(vk, 32, 32);
%     av_meas= (mean(by_stim_pair, 2));
%     plot(av_meas, 'b');
%     hold on;
%     scatter(bad_elecs, av_meas(bad_elecs), 'r', 'filled');
%     hold off;
%     xlim([1 32]);
%     title("Average Measurement for Measurement Pair");
%     legend('meas', 'bad elecs');
    
end % end function