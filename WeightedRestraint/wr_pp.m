function [D, imdl_comp] = wr_pp(files, imdl, timeSeriesLength)

if nargin == 2
    timeSeriesLength = inf;
end
% Weighted restraint pre-processing
clip = 50;
D = struct;
RMTHRESHOLD = 0.25;
msel = imdl.fwd_model.meas_select;
mm = find(msel);
sequences = cell(length(files), 1);

for i = 1:length(files)
    
    file = files{i};
    
    if isempty(file)
        continue
    end
    field = sprintf('seq_%0.f', i);
    [ dd, D.(field).aux ] = eidors_readdata(file);
    
    % visualize
%     inspect_eit_elec_and_data({dd, D.(field).aux} , imdl, RMTHRESHOLD);
%     sgtitle(remove_underscores(file));
%     savefigdir = 'C:\Users\Mark\Documents\GraduateStudies\LAB\EIT-restraint\zzMC\figures\data quality\electrode quality\';
%     printPDF( sprintf('%s%s.pdf', savefigdir, remove_underscores(file)) );
%     clf();

    % calculate framerate
    D.(field).aux.elec_impedance = D.(field).aux.elec_impedance(:, 1: size(dd, 2));
    t_rel = D.(field).aux.t_rel(:, 1: size(dd, 2));
    D.(field).fs = 1e6 ./ median( diff(t_rel) ); %framerate is median dif of time points/ 1000000 (convert to s)
    
    % filter
    dd = dd( :, 1: min(size(dd, 2), round(D.(field).fs * timeSeriesLength)) ); % ensure all time series are timeSeriesLength or less
    D.(field).data = dd;
    sequences{i} = dd;
    useData = real(dd(mm, :));
    useData = lowpass(useData', 1, D.(field).fs)'; % 1 Hz lowpass
    useData = useData(:, clip: (size(useData, 2) - clip) ); % trim filter edge artifacts
    D.(field).useData = useData;
end % end for

% noise compensation
specs = struct; 
specs.n = 6; 
specs.type = 'elec';
specs.thresh = RMTHRESHOLD;
[imdl_comp, ~, ~] = eqadr(D, imdl, specs);

end % end function