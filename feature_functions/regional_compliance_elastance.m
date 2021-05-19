function [rCompliance, rElastance] = regional_compliance_elastance(imgs, triads, Pdelta)
% -------------------------------------------------------------------------
% DESCRIPTION:
%   [rCompliance, rElastance] = regional_compliance_elastance(imgs, triads, Pdelta) 
% 
%   Assumes animal in imgs is in dorsal recumbancy.
% 
%   It is worth noting that fEIT images of aeration change also provide an
%   approximation of Crs, when normalized by the difference in airway
%   pressure between peak-inspiration and peak-expiration.
% 
%   Estimates compliance by calculating impedance change from end ex to end
%   in within the lungs relative to the pressure change between end ex and
%   end in for that breath (Arbitrary Units / Pressure)
% -------------------------------------------------------------------------
% PARAMETERS:
%   imgs:
%       EIT images for full time series
%   triads:
%       n x 3 array of indices in images for triad (end_ex1, end_in,
%       end_ex2)
%   pDelta:
%       n x 2 array of pressure differential between (end_ex1 and end_in)
%       and (end_in and end_ex2) for each breath in triads
% -------------------------------------------------------------------------   
% RETURNS:
%   rCompliance (struct):
%       
% -------------------------------------------------------------------------   
% AUTHOR:
%   Mark Campbell
%   Carleton University
%   markacampbell@sce.carleton.ca
% -------------------------------------------------------------------------
% VERSION:
%   1.0.0
% -------------------------------------------------------------------------
% (C) 2019-2021 Mark Campbell. License: GPL version 2 or version 3
% -------------------------------------------------------------------------

roi     = horse_roi();
Cimgs   = -imgs;  % impedance is high
Ctriads = triads;

% move end in to end_ex1 and end_ex2 to end_in to get difference image of
% end_ex2 - end_in
Eimgs       = imgs;   % impedance is low
Etriads     = triads;
Etriads(:,1)= triads(:,2);
Etriads(:,2)= triads(:,3);

rCompliance = do_compliance_elastance_calc(Cimgs, Ctriads, roi, Pdelta);
rElastance  = do_compliance_elastance_calc(Eimgs, Etriads, roi, Pdelta);

end % end function

function rCE = do_compliance_elastance_calc(imgs, triads, roi, Pdelta)
    
    nBreaths    = size(triads,1);
    bothLungs   = (roi.LeftLung + roi.RightLung) > 0;
    lungMask    = find(bothLungs);
    ltVent      = zeros(nBreaths,1);
    ltDors      = zeros(nBreaths,1);
    rtVent      = zeros(nBreaths,1);
    rtDors      = zeros(nBreaths,1);

    ltVentMask  = find(roi.LeftVentral);
    ltDorsMask  = find(roi.LeftDorsal);
    rtVentMask  = find(roi.RightVentral);
    rtDorsMask  = find(roi.RightDorsal);

    for i = 1:nBreaths
        triad       = triads(i, :);
        breathZ     = calc_breath_delta_z(imgs, triad);
        breathZbyP  = breathZ / Pdelta;
        
        lungTotal   = sum(breathZbyP(lungMask));
        ltVent(i)   = sum(breathZbyP(ltVentMask)) / lungTotal;
        ltDors(i)   = sum(breathZbyP(ltDorsMask)) / lungTotal;
        rtVent(i)   = sum(breathZbyP(rtVentMask)) / lungTotal;
        rtDors(i)   = sum(breathZbyP(rtDorsMask)) / lungTotal;
    end % end for
    
    rCE         = struct;
    rCE.ltVent  = ltVent;
    rCE.ltDors  = ltDors;
    rCE.rtVent  = rtVent;
    rCE.rtDors  = rtDors;

end % end function





