function breaths= findBreaths(vv)
% -------------------------------------------------------------------------
% Description:
%   breaths= findBreaths(vv)
% -------------------------------------------------------------------------
% Parameters:
%   vv:
%       EIT data from output of eidors_readdata
% ------------------------------------------------------------------------- 
% Returns:
%   breaths: A struct with fields
%       insIdx:
%           Indices of inspiration (peaks of total boundary voltage)
%       expIdx:
%           Indices of expiration (valleys of total boundary voltage)
% ------------------------------------------------------------------------- 
% Author:
%   Mark Campbell
%   Carleton University
%   markacampbell@cmail.carleton.ca
%   27.Sep.2019
% -------------------------------------------------------------------------
vv_size= size(vv);

if min(vv_size)> 1
   vv= sum(vv, 1); 
end % end if

vvMean= mean(vv); % zero data
inspir= vv >vvMean;
expir= ~inspir;
insIdx= [];
expIdx= [];
segMax= 0;
segMin= max(vv);
segMinRes= segMin;
for i= 1: length(inspir)
    if inspir(i)== 1
        if vv(i)> segMax
            segMax= vv(i);
            maxIdx= i;
        end % end if
    elseif inspir(i)== 0
        if segMax> 0
            insIdx= [insIdx, maxIdx];
            segMax= 0;
        end % end if
    end % end if
end % end for

for i= 1: length(expir)
    if expir(i)== 1
        if vv(i)< segMin
            segMin= vv(i);
            minIdx= i;
        end % end if
    elseif expir(i)== 0
        if segMin< segMinRes
            expIdx= [expIdx, minIdx];
            segMin= segMinRes;
        end % end if
    end % end if
end % end for

% ensure there are a whole number of breaths and vectors are same length
if length(insIdx) > length(expIdx)
    insIdx= insIdx(1: length(expIdx));
elseif length(insIdx) < length(expIdx)
    expIdx= expIdx(1: length(insIdx));
end % end if

breaths.insIdx= insIdx;
breaths.expIdx= expIdx;

end % end function