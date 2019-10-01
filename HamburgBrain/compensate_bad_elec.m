function imdl_comp= compensate_bad_elec(vv, imdl)

% Description:
%   imdl_comp= compensate_bad_elec(vv, imdl) 
%
%   Modify reconstruction matrix of imdl to compensate for noisy or
%   disconnected electrodes using the method from Mamatjan 2017.
%
% Parameters:
%   vv:         EIT data with or without complex component.
%   
% Returns:
%   vv_cleaned: real and complex component (if applicable) of EIT data with
%               measurements involving bad
%               elecs set to 0.
%   
% Author:
%   Mark Campbell
%   Carleton University
%   markacampbell@cmail.carleton.ca
%   27.Sep.2019

imdl_comp= imdl;

% Find bad electrodes
vk = 1e3* mean(real(vv), 2);
by_stim_pair= reshape(vk, 32, 32);
av_meas= (mean(by_stim_pair, 2))';
bad_elecs= find( abs( av_meas- mean(av_meas) )> std(av_meas));

% Find measurements from bad electrodes
kk=meas_icov_rm_elecs(imdl_comp, bad_elecs);
ee = find(diag(kk)~=1);

% do work on reconstruction matrix
imdl_comp.solve_use_matrix.RM_orig= imdl_comp.solve_use_matrix.RM;
PJt= imdl_comp.solve_use_matrix.PJt;
X= imdl_comp.solve_use_matrix.X;
X_star= X;

for i= 1:length(ee)
    channel= ee(i);
    Xi= X(:, channel);
    X_star= X_star- Xi'.* Xi* inv(X(channel, channel));
end % end for

imdl_comp.solve_use_matrix.X_star= X_star;
imdl_comp.solve_use_matrix.RM= PJt* X_star;

end % end function