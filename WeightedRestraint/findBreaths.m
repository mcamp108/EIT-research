function breaths= findBreaths(vv)

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

breaths.insIdx= insIdx;
breaths.expIdx= expIdx;
end % end function