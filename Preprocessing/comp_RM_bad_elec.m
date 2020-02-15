function imdl_comp= comp_RM_bad_elec(imdl, rm_elecs)
% -------------------------------------------------------------------------
% Description:
%   imdl_comp= comp_RM_bad_elec(imdl, rm_elecs)
%
%   Modify reconstruction matrix of imdl to compensate for noisy or
%   disconnected electrodes using the method from Mamatjan 2017. This
%   function will remove electrodes whose average contact impedance was
%   over 400 ohms. If more than 6 electrodes are identified, only the worst
%   6 are removed.
% -------------------------------------------------------------------------
% Parameters:
%   imdl:
%       EIDORS inverse model structure with keep intermediate results.
%   rm_elecs:
%       Electrodes to zero out from the reconstruction matrix.
% -------------------------------------------------------------------------   
% Returns:
%   imdl_comp:
%       inverse model whose reconstruction matrix has been modified with a
%       series of rank one updates that compensate for the removal of the
%       removed electrodes.
% -------------------------------------------------------------------------   
% Author:
%   Mark Campbell
%   Carleton University
%   markacampbell@cmail.carleton.ca
%   27.Sep.2019
% -------------------------------------------------------------------------

imdl_comp = imdl;

% Find measurements from bad electrodes
kk = meas_icov_rm_elecs(imdl_comp, rm_elecs);
ee = find(diag(kk)~=1); % bad channels of meas_sel

% do work on reconstruction matrix
imdl_comp.solve_use_matrix.RM_orig = imdl_comp.solve_use_matrix.RM;
PJt = imdl_comp.solve_use_matrix.PJt;
X = imdl_comp.solve_use_matrix.X;
X_star = X;

for i= 1:length(ee)
    channel= ee(i);
    Xi= X_star(:, channel);
    X_star= X_star- Xi* Xi' * inv(X_star(channel, channel)); % update on X_star each time. go back and explicitly zero bad meas
    X_star(channel, :)= 0;
    X_star(:, channel)= 0;
end % end for

imdl_comp.solve_use_matrix.X_star= X_star;
imdl_comp.solve_use_matrix.RM= PJt* X_star;

end % end function