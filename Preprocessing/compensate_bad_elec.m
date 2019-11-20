function varargout= compensate_bad_elec(eit_file, imdl, thresh)
% -------------------------------------------------------------------------
% Description:
%   imdl_comp= compensate_bad_elec(eit_file, imdl)
%
%   Modify reconstruction matrix of imdl to compensate for noisy or
%   disconnected electrodes using the method from Mamatjan 2017. This
%   function will remove electrodes whose average contact impedance was
%   over 400 ohms. If more than 6 electrodes are identified, only the worst
%   6 are removed.
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
%   imdl_comp:
%       inverse model whose reconstruction matrix has been modified with a
%       series of rank one updates that compensate for the removal of the
%       removed electrodes.
%   vv_prime: 
%       real and complex component (if applicable) of EIT data with
%       measurements involving bad elecs set to 0.
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

imdl_comp= imdl;

% find bad electrodes
impedanceFactor =  2.048 / (2^12 * 0.173 * 0.003) / 2^15; % = 0.9633 / 2^15;
elec_impedance= abs(auxdata.elec_impedance* impedanceFactor);
ei= median(elec_impedance, 2);
bad_elecs= find(ei> thresh);

if length(bad_elecs)> 6
    be_ci= sort(ei(bad_elecs), 'descend');
    bad_elecs= zeros(6, 1);
    for i= 1:6
        bad_elecs(i)= find(ei== be_ci(i));
    end % end for
end % end if

% Find measurements from bad electrodes
kk= meas_icov_rm_elecs(imdl_comp, bad_elecs);
ee = find(diag(kk)~=1); % bad channels of meas_sel


% % remove noisy channels from data
if nargout==2
    ge= find(diag(kk)==1); % good channels of meas_sel
    msel= imdl.fwd_model.meas_select;
    vv_prime= real(data(find(msel), :));
    varargout{2}= vv_prime(ge, :);
end % end if

% do work on reconstruction matrix
imdl_comp.solve_use_matrix.RM_orig= imdl_comp.solve_use_matrix.RM;
PJt= imdl_comp.solve_use_matrix.PJt;
X= imdl_comp.solve_use_matrix.X;
X_star= X;

for i= 1:length(ee)
    channel= ee(i);
    Xi= X_star(:, channel);
    X_star= X_star- Xi* Xi' * inv(X_star(channel, channel)); % update on X_star each time. go back and explicitly zero bad meas
    X_star(channel, :)= 0;
    X_star(:, channel)= 0;
end % end for

imdl_comp.solve_use_matrix.X_star= X_star;
imdl_comp.solve_use_matrix.RM= PJt* X_star;
varargout{1}= imdl_comp;
end % end function