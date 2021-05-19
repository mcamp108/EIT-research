function fmdl = shift_electrodes(fmdl, elecShift)
% -------------------------------------------------------------------------
% DESCRIPTION:
%   rotate electrode belt by shifting electrodes. + (clockwise) -
%   (anti-clockwise)
%   'lr': 
%       flip belt about horizontal axis, if electrode numbering ascended
%       clockwise, it will now ascend anti-clockwise.
%   'ud':
%       flip electrodes so electrode 1 move from front to back
% -------------------------------------------------------------------------
% PARAMETERS:
% 
% -------------------------------------------------------------------------   
% RETURNS:
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

oldElectrode = fmdl.electrode;
newElectrode = fmdl.electrode;
nElec = length(oldElectrode);
idx = (1:nElec);
if ischar(elecShift)
    if strcmp(elecShift, 'lr')
        oldIdx = (nElec+1) - idx;
    elseif strcmp(elecShift, 'ud')
        oldIdx = [fliplr(1:nElec/2), fliplr(nElec/2+1:nElec)];
    end
else
    if elecShift <= 0
        oldIdx = [idx(nElec+elecShift+1 : end), idx(1: nElec+elecShift)];
    elseif elecShift == 0
        oldIdx = idx;
    else
        oldIdx = [idx(elecShift+1:end), idx(1:elecShift)];
    end
end
for i = 1:nElec
    newElectrode(i) = oldElectrode( oldIdx(i) );
end
fmdl.electrode = newElectrode;

end % end function