function [imdl_comp, vv_prime]= compensate_bad_elec(vv, elec_impedance, imdl)
% -------------------------------------------------------------------------
% Description:
%   imdl_comp= compensate_bad_elec(vv, imdl) 
%
%   Modify reconstruction matrix of imdl to compensate for noisy or
%   disconnected electrodes using the method from Mamatjan 2017.
% -------------------------------------------------------------------------
% Parameters:
%   vv:         EIT data with or without complex component.
% -------------------------------------------------------------------------   
% Returns:
%   vv_cleaned: real and complex component (if applicable) of EIT data with
%               measurements involving bad
%               elecs set to 0.
% -------------------------------------------------------------------------   
% Author:
%   Mark Campbell
%   Carleton University
%   markacampbell@cmail.carleton.ca
%   27.Sep.2019
% -------------------------------------------------------------------------

imdl_comp= imdl;

% Find bad electrodes
% vk = 1e3* mean(real(vv), 2);
% by_stim_pair= reshape(vk, 32, 32);
ei= mean(real(abs(elec_impedance)), 2);
bad_elecs= find(ei> 400);

% Identify bad electrodes based on jumps in contact impedance
% d= abs(elec_impedance');
% dm= movmean(d, 5);
% ddt= detrend(dm);
% drng= range(ddt);
% bad_elecs= [bad_elecs; find(drng> 15)'];
% 
% bad_elecs= unique(bad_elecs);

% av_meas= (mean(by_stim_pair, 2))';
% bad_elecs= find( abs( av_meas- mean(av_meas) )> std(av_meas));

% Find measurements from bad electrodes
kk= meas_icov_rm_elecs(imdl_comp, bad_elecs);
ee = find(diag(kk)~=1); % bad channels of meas_sel
ge= find(diag(kk)==1); % good channels of meas_sel

% remove noisy channels from data
msel= imdl.fwd_model.meas_select;
vv_prime= real(vv(find(msel), :));
vv_prime= vv_prime(ge, :);

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

end % end function