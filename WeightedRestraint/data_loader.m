function [ vv , fs, auxdata, stim] = data_loader( fname ) 
% This function loads .eit files from a local directory 
%   The function outputs the data as vv and the 
%   samples/second of the file as fs

    fname = char(fname); % Required due to storing files as string in array
    % Recursively search target directories for names files
    server_data = dir(['svn/lab_data/2019/weighted_restraint_trial/**/',fname]);
    local_data  = dir(['../**/',fname]);
    i=0; 
    while 1; i=i+1; 
        switch i
            case 1; fpath = [local_data.folder,'/',fname]; % First check locally
            case 2; fpath = [server_data.folder,'/',fname]; % Then check server 
            case 3; error(['Cant find file ',fname]); 
        end
       try eidors_readdata(fpath, 'LQ4'); break; end % break if successful load
    end
    [vv, auxdata, stim]= eidors_readdata(fpath,'LQ4'); %  Swisstom EIT 2015+
    vv = real(vv); % Keep only the real component
    
    % Get the sampling rate of the file to return a time array
    tx = (abs(auxdata.t_rel(1)-auxdata.t_rel(2)));
    fs = 1000000/tx; % 1000000 is 1s in swisstom time
end

