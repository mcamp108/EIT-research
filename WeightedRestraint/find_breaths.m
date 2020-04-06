function breaths = find_breaths(vv)
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
% pp.min_insp_length = 0.1; % seconds
% pp.min_expi_length = 0.1; % seconds
% pp.FRC_relative_match = 0.5; % match start/end FRC
% pp.FRC_search_window = 0.3;  % search 100ms for FRC
% data.FR = FR;
% 
% lseq = length(seq);
% Fseq= fft(seq);
% Fcutoff = 0.25;
% L = round( Fcutoff * 2 * lseq / data.FR );
% Fseq([1+L+1:end-L])=0; %HPF
% Fseq([1,2,end])=0;     %LPF
% seq1= ifft(Fseq);
% seq1= real(seq1);
% % Flow calc
% flow = diff(seq1);
% thresh = median( abs(flow));
% inout= zeros(lseq,1);
% i = 1;
% inout(i) = sign( flow(i)+eps ); % eps to force not zero
% 
% for i = 2:lseq-1
%  cval= inout(i-1);
%  if cval*flow(i) > -thresh
%     inout(i) = cval;
%  else
%     inout(i) = -cval;
%  end
% end
% inout(lseq) = inout(lseq-1);
% 
% % Shorter forward window because detection is lagging
% breath_window_fwd= round(pp.FRC_search_window/4*data.FR);
% breath_window_bak= round(pp.FRC_search_window*data.FR);
% 
% dinout= diff( inout );
% fdiff = find( diff(inout) );
% fdiff(fdiff<=breath_window_bak      )= []; % too close
% fdiff(fdiff>lseq - breath_window_fwd)= []; % too close
% 
% eexpi= fdiff( (dinout(fdiff)>0) );
% einsp= fdiff( (dinout(fdiff)<0) );
% 
% % Find the best point
% ww= -breath_window_bak:breath_window_fwd;
% 
% for i=1:length(eexpi)
%  wind= seq( eexpi(i)+ ww );
%  ff = find( wind== min(wind) );
%  ff= ff(1)+min(ww)-1;
%  eexpi(i)= eexpi(i) + ff;
% end
% 
% for i=1:length(einsp)
%  wind= seq( einsp(i)+ ww );
%  ff = find( wind== max(wind) );
%  ff= ff(1)+min(ww)-1;
%  einsp(i)= einsp(i) + ff;
% end
% 
% min_insp_length = round(pp.min_insp_length*data.FR);
% min_expi_length = round(pp.min_expi_length*data.FR);
% breaths = [];
% i=1;e=1; 
% 
% while true
%   if i>=length(einsp) && e>=length(eexpi)-1; break; end
%   if einsp(i) < eexpi(e)
%      i=i+1;
%   else
%      rstr =  sprintf('rejecting breath (%d) [%i-%i-%i]: ', ...
%            i, eexpi(e), einsp(i), eexpi(e+1));
%      FRCs= seq1(eexpi(e+[0:1]));
%      TV  = seq1(einsp(i)) - mean(FRCs);
%      if einsp(i) - eexpi(e) < min_insp_length 
%         disp([rstr,'min_insp_length']);
%      elseif eexpi(e+1)-einsp(i) < min_expi_length
%         disp([rstr,'min_expi_length']);
%      elseif abs(diff(FRCs))/TV > pp.FRC_relative_match
%         disp([rstr,'FRC_relative_match']);
%      else %accept breath
%         breaths(end+1,:) = [eexpi(e), einsp(i), eexpi(e+1)]; 
%      end
%      i=i+1; e=e+1;
%   end
% end % end while
%    
% outbreaths.ins_idx= breaths(:,2);
% outbreaths.exp_idx= [breaths(:,1), breaths(:,3)];
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

end % end function