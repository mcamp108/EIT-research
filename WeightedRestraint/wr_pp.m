function [D, imdl_comp]= wr_pp(files, imdl, timeSeriesLength)

% Weighted restraint pre-processing
clip= 50;
D = struct;
RMTHRESHOLD = 0.25;
msel= imdl.fwd_model.meas_select;
mm = find(msel);
sequences = {};
for file = 1:length(files)
    f = files{file};
    field = horzcat('seq', num2str(file));
    [ dd, D.(field).aux ] = eidors_readdata(f);
    
    inspect_eit_elec_and_data({dd, D.(field).aux} , imdl, RMTHRESHOLD);
    sgtitle(remove_underscores(f));
%     savefigdir = 'C:\Users\Mark\Documents\GraduateStudies\LAB\EIT-restraint\zzMC\figures\data quality\electrode quality\';
%     printPDF( sprintf('%s%s data quality Resolve', savefigdir, remove_underscores(f)) );
%     keyboard;
    clf();
    
    D.(field).aux.elec_impedance = D.(field).aux.elec_impedance(:, 1: size(dd, 2));
    t_rel = D.(field).aux.t_rel(:, 1: size(dd, 2));
    D.(field).fs = 1e6 ./ median( diff(t_rel) ); %framerate is median dif of time points/ 1000000 (convert to s)
    dd = dd( :, 1: min(size(dd, 2), round(D.(field).fs * timeSeriesLength)) ); % ensure all time series are timeSeriesLength or less
    D.(field).data = dd;
    sequences{file} = dd;
    useData = real(dd(mm, :));
    useData = lowpass(useData', 1, D.(field).fs)'; % 1 Hz lowpass
    useData = useData(:, clip:(size(useData,2)-clip) ); % trim filter edge artifacts
    D.(field).useData = useData;
end % end for

% find worst N electrodes
[rmElecs, scores] = worst_n_elecs(sequences, imdl, 6);
rmElecs = rmElecs(scores >= RMTHRESHOLD); % change this
imdl_comp= comp_RM_bad_elec(imdl, rmElecs); % remove measurements from noisy electrodes by adjusting imdl
fprintf('electrodes removed with score threshold %s: %s \n', num2str(RMTHRESHOLD), num2str(rmElecs'));

end % end function