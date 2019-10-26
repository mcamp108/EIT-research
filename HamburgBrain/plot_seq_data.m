function plot_seq_data(seq, opt)
% This function was written by :
%                               Mark Campbell
%                               Carleton University
% opt can have fields:
%   pv
%     pv= 1 - plot peaks
%     pv= 2 - plot valleys
%     pv= 3 - plot peaks and valleys
%   usefdata
%     usefdata= 0 - use unfiltered EIT data
%     usefdata= 1 - use filtered EIT data
usefdata= 0;  
figure('units','normalized','outerposition',[0 0 1 1]);

if nargin> 1
    if isstruct(opt)
        hasOpt= true;
    else
        disp("options must be a struct. Options ignored");
        hasOpt= false;
    end % end if
end % end if
if hasOpt
    if isfield(opt, 'pv')
       plot_pv(seq, opt.pv);
    end % end if
    if isfield(opt, 'usefdata')
       usefdata= opt.usefdata;
    end % end if
end

plot_perf_data(seq);
firstplot= get(gca);
plot_sum_data(seq, usefdata);
xlim(firstplot.XLim);
end % end function

function plot_perf_data(seq)
    fs= seq.perf.tickrate;
    x= seq.perf.data;
    xax = (1: length(x))/fs;
    subplot(2,1,1);
    hold on
    plot(xax, x);
    xlabel("Time (s)");
    ylabel("Arterial Pressure (mmHg)");
    title("Arterial Blood Pressure: "+ seq.name);
    hold off
end % end function

function plot_sum_data(seq, usefdata)
    fs= seq.eit.fs;
    if usefdata== 0
        dat= real(sum(seq.eit.data, 1));
    else
        dat= real(sum(seq.eit.fdata, 1));
    end % end if
    xax = (1: size(dat, 2))/fs;
    tsMin= min(dat);
    tsMax= max(dat);
    subplot(2,1,2);
    plot(xax, dat);
    if isfield(seq.eit, 'apn')
        apn= seq.eit.apn/ fs;
        inj= seq.eit.inj/ fs;
        vnt= seq.eit.vnt/ fs;
        hold on 
        plot([apn, apn], [tsMin, tsMax]); 
        plot([inj, inj], [tsMin, tsMax]); 
        plot([vnt, vnt], [tsMin, tsMax]); 
        legend('Total Boundary Voltage', 'Start of Apnoea' , 'Bolus Injection', 'End of Apnoea');
    else
        legend('Total Boundary Voltage');
    end % end if
    
    xlabel("Time (s)");
    ylabel("Voltage (mV)");
    title("Total Boundary Voltage: "+ seq.name);
    
    hold off
end % end function

function plot_pv(seq, pOpt)
subplot(2, 1, 1)
hold on
if pOpt== 1 || 3
    for i= 1:length(seq.perf.peaks)
        loc= seq.perf.peaks(i);
        y= seq.perf.data(loc);
        plot([loc/ seq.perf.tickrate, loc/ seq.perf.tickrate], [y, y+50], 'r');
    end % end for
end % end if
if pOpt== 2 || 3
    for i= 1:length(seq.perf.vals)
        loc= seq.perf.vals(i);
        y= seq.perf.data(loc);
        plot([loc/ seq.perf.tickrate, loc/ seq.perf.tickrate], [y, y-50], 'g');
    end % end for
end % end if
if pOpt< 1 || pOpt> 3
    disp("Uncrecognized pv option");
end % end if
hold off

end % end function
