function breaths= find_breaths(vv)
% -------------------------------------------------------------------------
% Description:
%   breaths= find_breaths(vv)
% -------------------------------------------------------------------------
% Parameters:
%   vv:
%       EIT data from output of eidors_readdata or the global EIT signal
% ------------------------------------------------------------------------- 
% Returns:
%   breaths: 
%       A breath is defined as a wave that begins at an expiration
%       (valley), passes an inpiration (peak), and ends at an expiration
%       (valley). breaths is a struct with fields:
%       ins_idx:
%           1 x n_breaths array containing indices of inspiration (peaks of
%           total boundary voltage).
%       exp_idx:
%           2 x n_breaths array containing indices of expiration (valleys
%           of total boundary voltage) bounding the breath.
%       
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
vv= detrend(vv);
vv= vv- mean(vv); % zero data
ins_idx= peakfinder(vv, 0.001, [], 1);
exp_idx= peakfinder(vv, 0.001, [], -1);
% for i= 1: length(vv)
%     if inspir(i)== 1
%         if vv(i)> segMax
%             segMax= vv(i);
%             maxIdx= i;
%         end % end if
%     elseif inspir(i)== 0
%         if segMax> 0
%             ins_idx= [ins_idx, maxIdx];
%             segMax= 0;
%         end % end if
%     end % end if
% end % end for

% for i= 1: length(expir)
%     if expir(i)== 1
%         if vv(i)< segMin
%             segMin= vv(i);
%             minIdx= i;
%         end % end if
%     elseif expir(i)== 0
%         if segMin< segMinRes
%             exp_idx= [exp_idx, minIdx];
%             segMin= segMinRes;
%         end % end if
%     end % end if
% end % end for

% if landmark is first point may may be false positive
if ins_idx(1)== 1 || exp_idx(1)== 1
    ins_idx= ins_idx(2:end);
    exp_idx= exp_idx(2:end);
end % end if

% ensure that there is an expiration before the first inspiration and after
% the last inspiration
while exp_idx(1) > ins_idx(1)
    ins_idx= ins_idx(2:end);
end % end while

while ~(ins_idx(end) < exp_idx(end) && ins_idx(end) > exp_idx(end-1))
    ins_idx= ins_idx(1:end-1);
end % end while

if length(ins_idx) ~= length(exp_idx)-1
    disp("Something's up with Jack's breath detection.");
    keyboard;
end % end if

% prepare ouput
breaths.ins_idx= ins_idx;
breaths.exp_idx= zeros(2,length(ins_idx));
for i= 1:length(ins_idx)
    if exp_idx(i)< ins_idx(i) && exp_idx(i+1) > ins_idx(i)
        breaths.exp_idx(:,i)= [exp_idx(i), exp_idx(i+1)];
    end
end % end for i

end % end function