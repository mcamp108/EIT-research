function files = load_joaquin_data(rootDir, horse)
% -------------------------------------------------------------------------
% DESCRIPTION:
%   data = load_joaquin_data()
% -------------------------------------------------------------------------
% PARAMETERS:
%   rootDir (str):
%   horse (int):
% -------------------------------------------------------------------------   
% RETURNS:
%   files (struct):
%       EITfile objects
% -------------------------------------------------------------------------   
% AUTHOR:
%   Mark Campbell
%   Carleton University
%   markacampbell@sce.carleton.ca
% -------------------------------------------------------------------------
% VERSION:
%   1.0.0
% -------------------------------------------------------------------------
% load all files for a horse and store in struct
files = files_for_horse(rootDir, horse);
end % end function

function filesOut = files_for_horse(rootDir, horse)
    filesOut = struct;
    filesOut.info = struct;
    switch horse
        
        case 296439
            
            filesOut.info.elecShift = 'lr';
            files       = {'296439_baseline','296439_IE11_T15_min','296439_IE11_T25_min',...
                        '296439_IE11_Perfusion_apnea','296439_IE11_Perfusion_breathing','296439_IE13_T15_min',...
                        '296439_IE13_T25_min','296439_IE13_Perfusion_apnea','296439_IE13_Perfusion_breathing'};
            Pexp        = [1,  1,  1,  0, 0, 1,  1,  0, 0];
            Pin         = [17, 16, 17, 0, 0, 16, 16, 0, 0];
            startIdx    = [1,1,1,...
                        1, 1, 1,...
                        1, 1, 1];
            endIdx      = [-1,-1,-1,...
                        -1,-1,-1,...
                        -1, -1, -1];
            perfStart   = [0, 0, 0,...
                        1200, 951, 0,...
                        0, 1300, 2110];
            perfEnd     = [0, 0, 0,...
                        2200, 2242, 0,...
                        0, 2100, 3845];
                    
        case 296695
            
            filesOut.info.elecShift = 'lr';
            files       = {'296695_Baseline','296695_IE_11_T15_min','296695_IE11_Perfusion_apnea',...
                        '296695_IE11_Perfusion_breathing','296695_IE11_T25_min','296695_IE13_Perfusion_Apnea',...
                        '296695_IE13_Perfusion_breathing','296695_IE13_T15_min','296695_IE13_T25_min'};
            Pexp        = [1,  1,  0, 0, 1,  0, 0, 1,  1];
            Pin         = [20, 20, 0, 0, 20, 0, 0, 22, 22];
            startIdx    = [1,1,1,...
                        1,1,1,...
                        1,1,1];
            endIdx      = [-1,-1,-1,...
                        -1,-1,-1,...
                        -1,-1,-1];
            perfStart   = [0, 0, 1000,...
                        1200, 0, 680,...
                        1200, 0, 0];
            perfEnd     = [0, 0, 2500,...
                        2300, 0, 2550,...
                        2400, 0, 0];
            
        case 296954
            
            filesOut.info.elecShift = -4;
            files       = {'296954_baseline', '296954_IE13_T25_min', '296954_IE13_T15_min',...
                        '296954_IE13_Perfusion_breathing', '296954_IE11_T25_min', '296954_IE11_T15_min',...
                        '296954_IE11_Perfusion_breathing', '296954_IE11_Perfusion_apnea', '296954_IE_13_Perfusion_apnea'};
            Pexp        = [1,  1,  1,  0, 0,  2,  0, 0, 0];
            Pin         = [22, 23, 23, 0, 25, 22, 0, 0, 0];
            startIdx    = [1,1,1,...
                        1,250,1,...
                        1,1,1]; 
            endIdx      = [7424,-1,-1,...
                        -1,-1,-1,...
                        -1,4100,-1];
            perfStart   = [0, 0, 0,...
                        1200, 0, 0,...
                        530, 900, 1400];
            perfEnd     = [0, 0, 0,...
                        250+2400, 0, 0,...
                        2500, 3200, 2900];
            
        case 296955
            
            filesOut.info.elecShift = 0;
            files       = {'296955_baseline','296955_IE_11_perfusion_apnea','296955_IE_11_Perfusion_breathing',...
                        '296955_IE_13_Perfusion_apnea','296955_IE13_T25_min','296955_IE11_T5_min',...
                        '296955_IE11_T10_min','296955_IE11_T25_min','296955_IE13_Perfusion_breathing',...
                        '296955_IE13_T15_min'};
            Pexp        = [1,  0, 0, 0,  1, 0,  1,  1, 0, 28];
            Pin         = [27, 0, 0, 0, 26, 0, 26, 26, 0,  0];
            startIdx    = [1,1,1,...
                        1,1,1,...
                        1,1,1,...
                        1];
            endIdx      = [7313,-1,-1,...
                        -1,-1,-1,...
                        -1,-1,-1,...
                        -1];
            perfStart   = [0, 520, 1600,...
                        600, 0, 0,...
                        0, 0, 1600,...
                        0];
            perfEnd     = [0, 1920, 2700,...
                        2000, 0, 0,...
                        0, 0, 2700,...
                        0];
            
    end
    for f = 1:length(files)
        file = sprintf('%s%.0f%s%s.eit',rootDir,horse,'\',files{f});
        name = strsplit(files{f}, sprintf('%.0f_',horse));
        filesOut.(name{end})            = load_eit(file,startIdx(f),endIdx(f));
        filesOut.(name{end}).perfStart  = perfStart(f);
        filesOut.(name{end}).perfEnd    = perfEnd(f);
        filesOut.(name{end}).Pexp       = Pexp(f);
        filesOut.(name{end}).Pdelta     = Pin(f);
        filesOut.(name{end}).name       = remove_underscores( files{f} );
        filesOut.(name{end}).horseID    = horse;
    end % end for
end

function EIT = load_eit(file, startIdx, endIdx)
    [vv,aux]= eidors_readdata(file);
    if endIdx == -1
        endIdx = size(vv,2);
    end
    t_rel   = aux.t_rel(:, startIdx:endIdx);
    fs      = 1e6 ./ median( diff(t_rel) ); %framerate is median dif of time points/ 1000000 (convert to s)
    vv      = vv(:, startIdx:endIdx);
    clip    = 50;
    rvv     = real(vv);
    if contains(file,'296954')
       rvv = movmedian(rvv, 23, 2); 
    end
    fdata   = lowpass(rvv',1, fs,'StopbandAttenuation',40,'ImpulseResponse','iir')';
    fdata   = fdata(:, clip: (size(fdata,2) - clip) ); % trim filter edge artifacts
    EIT     = struct;
    EIT.fs  = fs;
    EIT.data    = vv;
    EIT.fdata   = fdata;
end % end function


