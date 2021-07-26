function roi = horse_roi(spec)
% -------------------------------------------------------------------------
% DESCRIPTION:
%   roi = horse_roi()
% -------------------------------------------------------------------------   
% RETURNS:
%   roi (struct):
%       inside, left lung, right lung, dorsal, ventral
% -------------------------------------------------------------------------   
% AUTHOR:
%   Mark Campbell
%   Carleton University
%   markacampbell@cmail.carleton.ca
% -------------------------------------------------------------------------
% VERSION:
%   2.0.0
% -------------------------------------------------------------------------
roi = struct;
if strcmp(spec, 'vq_full')
   roi.BothLungs = vq_full_roi();
elseif strcmp(spec, 'v2')
    roi.BothLungs = both_lungs();
    roi.RightLung = right_lung();
    roi.LeftLung = left_lung();
else
    % old version
    data = load('C:\Users\Mark\Documents\EIT\demo\sedation_2_21_9_STEM.mat');
    temp = data.data.patient.ROI;
    temp.LeftLung(:, 17) = 0;
    roi.BothLungs = temp.RightLung + temp.LeftLung > 0;
    roi.RightLung = temp.RightLung > 0;
    roi.LeftLung = temp.LeftLung > 0;

    lungRows    = find(sum(roi.BothLungs, 2));
    lungHeight  = length(lungRows);
    endVentral  = round(lungHeight / 2);

    roi.Ventral = zeros(size(roi.RightLung));
    roi.Dorsal  = zeros(size(roi.RightLung));
    roi.Ventral(lungRows(1:endVentral), :) = roi.BothLungs(lungRows(1:endVentral), :);
    roi.Dorsal(lungRows(endVentral+1:end), :) = roi.BothLungs(lungRows(endVentral + 1: end), :);

    roi.LeftVentral     = (roi.LeftLung  + roi.Ventral) == 2;
    roi.LeftDorsal      = (roi.LeftLung  + roi.Dorsal)  == 2;
    roi.RightVentral    = (roi.RightLung + roi.Ventral) == 2;
    roi.RightDorsal     = (roi.RightLung + roi.Dorsal)  == 2;    
end

end % end function

function roi_test()
    roi = horse_roi();
    figure();
    combined = roi.LeftVentral + roi.LeftDorsal*2 + roi.RightVentral*3 + roi.RightDorsal*4;
    imagesc(combined);
    colorbar
    axis equal;
    
    figure();
    imagesc(roi.LeftLung+roi.RightLung);
    colorbar
    axis equal;
end

function BL = both_lungs()
   BL = [   
   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
   0   0   0   0   0   0   0   0   0   0   0   1   1   1   0   0   0   0   0   0   1   1   0   0   0   0   0   0   0   0   0   0
   0   0   0   0   0   0   0   0   0   0   1   1   1   1   1   0   0   1   1   1   1   1   1   0   0   0   0   0   0   0   0   0
   0   0   0   0   0   0   0   0   0   1   1   1   1   1   1   1   0   1   1   1   1   1   1   1   0   0   0   0   0   0   0   0
   0   0   0   0   0   0   0   0   0   1   1   1   1   1   1   1   0   1   1   1   1   1   1   1   0   0   0   0   0   0   0   0
   0   0   0   0   0   0   0   0   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   0   0   0   0   0   0   0
   0   0   0   0   0   0   0   0   1   1   1   1   1   1   1   1   0   1   1   1   1   1   1   1   1   0   0   0   0   0   0   0
   0   0   0   0   0   0   0   0   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   0   0   0   0   0   0   0
   0   0   0   0   0   0   0   0   1   1   1   1   1   1   1   1   0   1   1   1   1   1   1   1   1   0   0   0   0   0   0   0
   0   0   0   0   0   0   0   0   1   1   1   1   1   1   1   1   0   1   1   1   1   1   1   1   1   0   0   0   0   0   0   0
   0   0   0   0   0   0   0   1   1   1   1   1   1   1   1   1   0   1   1   1   1   1   1   1   1   0   0   0   0   0   0   0
   0   0   0   0   0   0   0   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   0   0   0   0   0   0   0
   0   0   0   0   0   0   0   0   1   1   1   1   1   1   1   1   0   1   1   1   1   1   1   1   1   0   0   0   0   0   0   0
   0   0   0   0   0   0   0   0   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   0   0   0   0   0   0   0
   0   0   0   0   0   0   0   0   1   1   1   1   1   1   1   1   0   1   1   1   1   1   1   1   1   0   0   0   0   0   0   0
   0   0   0   0   0   0   0   0   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   0   0   0   0   0   0   0
   0   0   0   0   0   0   0   0   0   1   1   1   1   1   1   1   0   1   1   1   1   1   1   1   0   0   0   0   0   0   0   0
   0   0   0   0   0   0   0   0   0   1   1   1   1   1   1   1   0   1   1   1   1   1   1   1   0   0   0   0   0   0   0   0
   0   0   0   0   0   0   0   0   0   0   0   1   1   1   1   0   0   0   1   1   1   1   1   0   0   0   0   0   0   0   0   0
   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0];
    BL = logical(BL);

end

function RL = right_lung()
   RL = [
   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
   0   0   0   0   0   0   0   0   0   0   0   1   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
   0   0   0   0   0   0   0   0   0   0   1   1   1   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
   0   0   0   0   0   0   0   0   0   1   1   1   1   1   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
   0   0   0   0   0   0   0   0   0   1   1   1   1   1   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
   0   0   0   0   0   0   0   0   1   1   1   1   1   1   1   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
   0   0   0   0   0   0   0   0   1   1   1   1   1   1   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
   0   0   0   0   0   0   0   0   1   1   1   1   1   1   1   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
   0   0   0   0   0   0   0   0   1   1   1   1   1   1   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
   0   0   0   0   0   0   0   0   1   1   1   1   1   1   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
   0   0   0   0   0   0   0   1   1   1   1   1   1   1   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
   0   0   0   0   0   0   0   1   1   1   1   1   1   1   1   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
   0   0   0   0   0   0   0   0   1   1   1   1   1   1   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
   0   0   0   0   0   0   0   0   1   1   1   1   1   1   1   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
   0   0   0   0   0   0   0   0   1   1   1   1   1   1   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
   0   0   0   0   0   0   0   0   1   1   1   1   1   1   1   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
   0   0   0   0   0   0   0   0   0   1   1   1   1   1   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
   0   0   0   0   0   0   0   0   0   1   1   1   1   1   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
   0   0   0   0   0   0   0   0   0   0   0   1   1   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0];
    RL = logical(RL);

end

function LL = left_lung()

BL = both_lungs();
RL = right_lung();
temp = BL + RL;
LL = temp == 1;

end

function roi_vq = vq_full_roi()
    roi_vq = [
    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0
    0   0   0   0   0   0   0   0   0   0   1   1   1   1   0   0   0   0   0   1   1   1   1   0   0   0   0   0   0   0   0   0
    0   0   0   0   0   0   0   0   0   1   1   1   1   1   0   0   0   0   0   1   1   1   1   1   0   0   0   0   0   0   0   0
    0   0   0   0   0   0   0   0   0   1   1   1   1   1   1   0   0   0   0   1   1   1   1   1   0   0   0   0   0   0   0   0
    0   0   0   0   0   0   0   0   1   1   1   1   1   1   0   0   0   0   0   1   1   1   1   1   1   0   0   0   0   0   0   0
    0   0   0   0   0   0   0   0   1   1   1   1   1   1   0   0   0   0   0   1   1   1   1   1   1   0   0   0   0   0   0   0
    0   0   0   0   0   0   0   0   1   1   1   1   1   1   0   0   0   0   0   1   1   1   1   1   1   0   0   0   0   0   0   0
    0   0   0   0   0   0   0   1   1   1   1   1   1   1   0   0   0   0   0   1   1   1   1   1   1   0   0   0   0   0   0   0
    0   0   0   0   0   0   0   1   1   1   1   1   1   1   0   0   0   0   0   1   1   1   1   1   1   0   0   0   0   0   0   0
    0   0   0   0   0   0   0   1   1   1   1   1   1   1   0   0   0   0   0   1   1   1   1   1   1   0   0   0   0   0   0   0
    0   0   0   0   0   0   0   1   1   1   1   1   1   1   0   0   0   0   0   1   1   1   1   1   1   0   0   0   0   0   0   0
    0   0   0   0   0   0   0   1   1   1   1   1   1   1   0   0   0   0   0   1   1   1   1   1   1   0   0   0   0   0   0   0
    0   0   0   0   0   0   0   1   1   1   1   1   1   1   0   0   0   0   0   1   1   1   1   1   1   0   0   0   0   0   0   0
    0   0   0   0   0   0   0   1   1   1   1   1   1   1   0   0   0   0   0   1   1   1   1   1   1   0   0   0   0   0   0   0
    0   0   0   0   0   0   0   0   1   1   1   1   1   1   1   0   0   0   1   1   1   1   1   1   0   0   0   0   0   0   0   0
    0   0   0   0   0   0   0   0   1   1   1   1   1   1   1   1   0   0   1   1   1   1   1   1   0   0   0   0   0   0   0   0
    0   0   0   0   0   0   0   0   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   0   0   0   0   0   0   0   0
    0   0   0   0   0   0   0   0   0   1   1   1   1   1   1   1   0   1   1   1   1   1   1   1   0   0   0   0   0   0   0   0
    0   0   0   0   0   0   0   0   0   0   1   1   1   1   1   1   0   1   1   1   1   1   1   1   0   0   0   0   0   0   0   0
    0   0   0   0   0   0   0   0   0   0   0   1   1   0   0   0   0   0   0   1   1   1   0   0   0   0   0   0   0   0   0   0
    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0];
end