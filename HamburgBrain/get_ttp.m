function out_img= get_ttp(seq, ref, opt)
% -------------------------------------------------------------------------
% DESCRIPTION:
%
%   out_img= get_ttp(seq, ref, opt)
%
%   Get an image showing the time taken for each pixel in the image to
%   reach its peak value in the ensemble-averaged cardiac cycle.
% 
%   This function calls get_ensembles with opt.sel= [0 0 1]
% -------------------------------------------------------------------------
% PARAMETERS:
%   seq:     a swine sequence struct 
%   ref:     which perfusion reference landmark to use. ref= 1(peaks) or
%            ref = 2 (valleys) opt: options
%   opt:
%       ensemble:
%           each
%           one
% -------------------------------------------------------------------------
% RETURNS:
%   out_img: the image whose pixel values represent the time to peak, in
%            number of frames.
% -------------------------------------------------------------------------
% AUTHOR:
%   Mark Campbell
%   Carleton University
%   markacampbell@cmail.carleton.ca
%   23.Sep.2019
% -------------------------------------------------------------------------

opt.sel= [0 0 1];
ensemble= get_ensembles(seq, ref, opt);
imgr_e= ensemble.imgr_ensemble;

% which dimension of ensemble to average

if ~isfield(opt, 'ensemble')
    take_mean_of_dim= 3; % set to 'each'
elseif strcmp(opt.ensemble, 'one')
    take_mean_of_dim= 2;
else
    disp("Unrecognized ensemble average dimension");
end % end if

imgr_e= squeeze(mean(imgr_e, take_mean_of_dim));

ttp= zeros(1, size(imgr_e, 1));
for row= 1: size(imgr_e, 1)
    this_row= imgr_e(row, :);
    ttp(row)= var(this_row);
%     ttp(row)= max(this_row)- mean(this_row);
%     ttp(row)= mean(this_row)- min(this_row);
%     ttp(row)= find(this_row== max(this_row));
%     ttp(row)= mean(this_row);
%     ttp(row)= range(this_row);
%     ttp(row)= std(this_row);
end % end for

% xax= 1:size(ensemble_t, 2);
% hold on
% for i= 1:20
%     plot(xax, ensemble_t(i, :));
% end % end if
% hold off
    
ttp= ttp./ max(ttp);
imgr= seq.imgr;
imgr.elem_data= ttp;
clim= size(imgr_e, 2);
% cmin= NaN;
imgr.calc_colours.ref_level= 0;
imgr.calc_colours.clim= clim;
imgr.calc_colours.lim= clim;
out_img= calc_slices(imgr);


end % end function