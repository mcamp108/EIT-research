function triads = horse_breath_finder(vv)
% -------------------------------------------------------------------------
% DESCRIPTION:
%
%
%
% -------------------------------------------------------------------------
% PARAMETERS:
% 
%
%
% -------------------------------------------------------------------------   
% RETURNS:
% 
%
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
vv_size = size(vv);
if min(vv_size)> 1
   vv= sum(vv, 1); 
end % end if

vv= movmean(vv, 20);

% vv= vv- mean(vv); % zero data
ins_idx= peakfinder(vv, std(vv)/3, [], 1);
exp_idx= peakfinder(vv, std(vv)/3, [], -1);

% if landmark is first point may be false positive
if ins_idx(1)== 1
    ins_idx= ins_idx(2:end);
elseif exp_idx(1)== 1
    exp_idx= exp_idx(2:end);
end % end if

if ins_idx(end) == length(vv)
    ins_idx= ins_idx(1:end-1);
elseif exp_idx(end) == length(vv)
    exp_idx= exp_idx(1:end-1);
end

% ensure that there is an expiration before the first inspiration and after
% the last inspiration
if isempty(ins_idx) || length(exp_idx) == 1
    triads = [];
    return
end

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

for i = 1:length(ins_idx) 
    if exp_idx(i)< ins_idx(i) && exp_idx(i+1) > ins_idx(i)
        if (i < length(ins_idx)) && ( exp_idx(i+2) < ins_idx(i+1) ) && ( vv(exp_idx(i+2)) < vv(exp_idx(i+1)) )
            breaths.exp_idx(:,i)= [exp_idx(i), exp_idx(i+2)];
        else
            breaths.exp_idx(:,i)= [exp_idx(i), exp_idx(i+1)];
        end
    end
end % end for i

% prepare ouput
breaths.ins_idx= ins_idx;
breaths.exp_idx= zeros(2,length(ins_idx));
for i= 1:length(ins_idx)
    if exp_idx(i)< ins_idx(i) && exp_idx(i+1) > ins_idx(i)
        breaths.exp_idx(:,i)= [exp_idx(i), exp_idx(i+1)];
    end
end % end for i
endEx1      = breaths.exp_idx(1,:)';
endEx2      = breaths.exp_idx(2,:)';
endIn       = breaths.ins_idx';
triadsOld   = [endEx1,endIn,endEx2];
triads = adjust_breath_landmarks(vv, triadsOld);
end % end function


function newTriads = adjust_breath_landmarks(tbv, triads)
    newTriads   = triads;
    nBreaths    = size(triads, 1);
    for i=1: nBreaths
        triad       = triads(i,:);
        breathLen   = triad(3) - triad(1) + 1;
        
        rising      = tbv(triad(1):triad(2));
        riseNorm    = normalize(rising, 'range');
        fitRiseIdx  = find((riseNorm > 0.1) + (riseNorm < 0.9) == 2);
        riseBase    = mean(riseNorm(riseNorm <= 0.05));
        L1          = polyfit(fitRiseIdx, riseNorm(fitRiseIdx), 1);
        yL1         = polyval(L1, 1:breathLen);
        newEndEx1   = find( abs(riseBase - yL1) == min(abs(riseBase - yL1)) );
        
        if triad(1)+newEndEx1-1 >= triad(2)
            newEndEx1 = length(rising) - 1;
        end
        if newEndEx1 > 1
            % find local minimum
            temp    = diff(rising) < 0;
            temp2   = find(temp(1:newEndEx1-1),1,'last') + 1;
            if isempty(temp2)
                newEndEx1 = 1;
            else
                newEndEx1 = temp2;
            end
            newTriads(i,1) = triads(i,1) + newEndEx1 - 1;
        end
        if i > 1
            newTriads(i-1, 3) = newTriads(i,1);
        end        
    end
end
