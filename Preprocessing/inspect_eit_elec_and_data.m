function inspect_eit_elec_and_data(eit_file, imdl, thresh)
% -------------------------------------------------------------------------
% Description:
%   inspect_eit_elec_and_data(eit_file, imdl, thresh)
% -------------------------------------------------------------------------
% Parameters:
%   eit_file:
%       Target data file ending in .eit file extension or a cell containing
%       the outpus of eidors_readdata {data, auxdata}.
%   imdl:
%       EIDORS inverse model structure
%   thresh:
%       Electrode score at which to label an electrode poor.
% ------------------------------------------------------------------------- 
% Returns:
%   Total boundary voltage plot
%   IQ plot
%   U-shape plot
%   Frequency spectrum plot
%   Contact impedance plot
%   Average measurement for measurement pair plot
%   Highlights low-quality electrodes and measurements
% ------------------------------------------------------------------------- 
% Author:
%   Mark Campbell
%   Carleton University
%   markacampbell@cmail.carleton.ca
%   27.Sep.2019
% -------------------------------------------------------------------------
% VERSION:
%   1.2.0
% -------------------------------------------------------------------------
% (C) 2019-2021 Mark Campbell. License: GPL version 2 or version 3
% -------------------------------------------------------------------------
msel= imdl.fwd_model.meas_select;
mm  = find(msel);

if nargin == 2
   thresh = 0.25; 
end % end if

impedanceFactor =  2.048 / (2^12 * 0.173 * 0.003) / 2^15; % = 0.9633 / 2^15;
if isstruct(eit_file) && strcmp(thresh, 'hamburg')
    data = eit_file.data;
    elec_impedance = eit_file.elec_impedance;
    fs = eit_file.fs;
    thresh = 0.25;
elseif iscell(eit_file) % input is output of eidors_read_data
    if length(eit_file) == 2
        if isstruct(eit_file{2})
            data = eit_file{1};
            auxdata = eit_file{2};
        else
            auxdata = eit_file{1};
            data = eit_file{2};
        end % end if
        fs = 1e6./median(diff(auxdata.t_rel));     % framerate is median dif of time points/ 1000000 (convert to s)
        
        elec_impedance = auxdata.elec_impedance* impedanceFactor;
    else
        disp("Whoops!");
    end % end if
elseif strcmp(eit_file(end-3:end), '.eit')
    [data, auxdata] = eidors_readdata(eit_file);
    fs = 1e6./median(diff(auxdata.t_rel));         % framerate is median dif of time points/ 1000000 (convert to s)
    elec_impedance = auxdata.elec_impedance * impedanceFactor;
else
    data = eit_file;
end % end if

n_samples = size(data, 2);
data_samp = round(1: (n_samples / 100): n_samples);
elec_impedance = abs(elec_impedance);

specs = struct; 
specs.n = 32;
specs.type = 'elec';
specs.thresh = thresh;
[badElecs, scores] = worst_n_elecs({data}, imdl, specs);
badElecs = badElecs(scores > thresh);

specs.n = 544;
specs.type = 'meas';
[badMeas, scores] = worst_n_elecs({data}, imdl, specs);
badMeas = badMeas(scores > thresh);
% add plot bad measurements
if length(badElecs) >= 1
    ee1 = electrodes_used_meas_or_stim(badElecs, imdl);
end % end if
ee = badMeas;
ei = mean(elec_impedance, 2);

% total boundary voltage
% subplot(5,2,1:2);
%     vv = detrend(real(data));
%     xax= (1:size(vv, 2))/ fs;
%     plot(xax, sum(vv(msel,:), 1));
%     ylabel('Voltage (V)')
%     xlabel('Time (s)');
%     title('Total Boundary Voltage (Raw)');

% U shapes 
% subplot(524);
%     vk = 1e3 * mean(vv, 2);
%     plot(vk,'k');    % all meas
%     hold on;
%     vk(~msel,:) = NaN;
%     plot(vk,'b');   % selected meas
%     vn = NaN * vk; vn(ee1,:) = vk(ee1,:); % measurements from rejected electrodes
%     plot(vn,'mo');
%     vn = NaN * vk; vn(ee,:) = vk(ee,:); % independently rejected measurements
%     plot(vn,'ro');
%     hold off;
%     box off; 
%     xlim([0,400]);
%     title 'U shapes';

% IQ plot
subplot(4,2,[1,3]);%subplot(5,2,[3,5]);
    plot(1e3 * data(:, data_samp), 'k+'); % all meas
    hold on;
    good = ~ismember(1:length(mm), ee1);
    
    plot(1e3 * data(mm(ee1),data_samp), 'r+');  % measurements from rejected electrodes (covers rejected meas)
    plot(1e3 * data(mm(good), data_samp), 'b+'); % used meas
    hold off; 
    box off;
    axis equal;
    title 'IQ plot - rejected electrodes';
subplot(4,2,[5,7]);%subplot(5,2,[7,9]);
    plot(1e3 * data(:, data_samp), 'k+'); % all meas
    hold on;
    plot(1e3 * data(msel, data_samp), 'b+'); % used meas
    plot(1e3 * data(mm(ee),data_samp), 'r+');   % independently rejected measurements
    hold off; 
    box off;
    axis equal;
    title 'IQ plot - rejected measurements';

% median contact impedance
subplot(4,2,[2,4]);
    b           = bar(ei,'FaceColor','flat');
    b.CData(:,:)= repmat([0 0 1], 32, 1);
    be          = badElecs;
    for e= 1:length(be)
        b.CData(be(e), :)= [1 0 0];
    end
    hold on;
    eb = errorbar(1:32, ei, max(elec_impedance, [], 2)-ei, min(elec_impedance, [], 2)- ei);
    set(eb,'Color',[0, 0, 0],'LineStyle','none');
    hold off;
    box off; 
    xlim([0, 33]);
    title 'Median Electrode Impedance';

% Fourier Series
subplot(4,2,[6,8]);
    ll  = size(data, 2);    
    ft  = fft(detrend(data.').', [], 2);
    fax = linspace(0, fs, ll+1); fax(end)=[];
    semilogy(fax, mean(abs(ft)));
    ylim([1e-4, 1]);
    box off; 
    xlim([0, fs/2]);
    title 'Fourier Spectrum (Hz)';

end % end function


function wasUsed = electrodes_used_meas_or_stim(rm_elecs, imdl)
    % find which measurements involved the electrodes to remove (from
    % meas_icov_rm_elecs)
%     meas_icov = [];
%     for stim = imdl.fwd_model.stimulation(:)'
%        mp = stim.meas_pattern;
%        sp = stim.stim_pattern;
%        icovi = ones(size(mp, 1), 1);
%        if any(sp(rm_elecs) ~= 0)
%           icovi = 0 * icovi;
%        else
%           icovi = ~any( mp(:, rm_elecs) ~= 0, 2);
%        end
%        meas_icov = [meas_icov; icovi];
%     end
%     wasUsed = find(diag(meas_icov) ~= 1);
    
    nElecs = length(imdl.fwd_model.electrode);
    nMeas = sum(imdl.fwd_model.meas_select);    
    wasUsed = zeros(nMeas, 1);
    measPerPair = nMeas / nElecs;
    
    count = 1;
    for i = 1: nElecs
        [r, c] = find(imdl.fwd_model.stimulation(i).meas_pattern);
        for j = 1:max(r)
            temp = ismember(rm_elecs, c(r == j)');
            wasUsed(count) = sum(temp) > 0;
            count = count + 1;
        end % end for
    end % end for
    for i = 1: length(rm_elecs)
        thisElec = rm_elecs(i);
        stop  = thisElec * measPerPair;
        start = stop - measPerPair + 1;
        wasUsed(start : stop) = 1;
    end
    wasUsed = find(wasUsed);
end