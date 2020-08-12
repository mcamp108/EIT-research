function plot_seq_data(seq, opt)
% This function was written by :
%                               Mark Campbell
%                               Carleton University
% opt can have fields:
%   pv
%     pv= 1 - plot peaks
%     pv= 2 - plot valleys
%     pv= 3 - plot peaks and valleys

figure('units','normalized','outerposition',[0 0 1 1]);

if nargin == 1 || ~isfield(opt, 'pv')
    opt.pv = 3;
end % end if

subplot(3, 1, 1);
    plot_perf_data(seq);
    plot_pv(seq, opt.pv);

firstplot= get(gca);

subplot(3, 1, 2);
%     hold on
    plot_sum_Z(seq);
    xlim(firstplot.XLim);
subplot(3, 1, 3);
%     hold on
    plot_sum_V(seq);
    xlim(firstplot.XLim);
    
end % end function

function plot_perf_data(seq)

    fs= seq.perf.tickrate;
    x= seq.perf.data;
    xax = (1: length(x))/fs;
    plot(xax, x);
    xlabel("Time (s)");
    ylabel("Arterial Pressure (mmHg)");
    title("Arterial Blood Pressure: "+ seq.name);
    
end % end function

function plot_sum_Z(seq)

    fs= seq.eit.fs;
    dat= sum(seq.imgr.elem_data, 1);
    xax = (1: size(dat, 2))/fs;
    tsMin= min(dat);
    tsMax= max(dat);
    plot(xax, dat);
    if isfield(seq.eit, 'apn')
        apn= seq.eit.apn/ fs;
        inj= seq.eit.inj/ fs;
        vnt= seq.eit.vnt/ fs;
        hold on;
        plot([apn, apn], [tsMin, tsMax]); 
        plot([inj, inj], [tsMin, tsMax]); 
        plot([vnt, vnt], [tsMin, tsMax]); 
        legend('Global Signal (Delta Z)', 'Start of Apnoea' , 'Bolus Injection', 'End of Apnoea');
    else
        legend('Global Signal (Delta Z)');
    end % end if
    
    xlabel("Time (s)");
    ylabel("Conductivity (Delta Z)");
    title("Global Signal: "+ seq.name);
    hold off;
    
end % end function

function plot_pv(seq, pOpt)

hold on;
if pOpt== 1 || 3
    for i= 1:length(seq.perf.peaks)
        loc= seq.perf.peaks(i);
        y= seq.perf.data(loc);
        plot([loc/ seq.perf.tickrate, loc/ seq.perf.tickrate], [y, y+50], 'r');
    end % end for
end % end if
if pOpt== 2 || 3
    vals = unique( seq.perf.vals(:) );
    for i= 1:length(vals)
        loc= vals(i);
        y= seq.perf.data(loc);
        plot([loc/ seq.perf.tickrate, loc/ seq.perf.tickrate], [y, y-50], 'g');
    end % end for
end % end if
if pOpt< 1 || pOpt> 3
    disp("Uncrecognized pv option");
end % end if
hold off;

end % end function

function plot_sum_V(seq)
   
    fs= seq.eit.fs;
    dat= sum(real(seq.eit.fdata), 1);
    xax = (1: size(dat, 2))/fs;
    tsMin= min(dat);
    tsMax= max(dat);
    plot(xax, dat);
    if isfield(seq.eit, 'apn')
        apn= seq.eit.apn/ fs;
        inj= seq.eit.inj/ fs;
        vnt= seq.eit.vnt/ fs;
        hold on;
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
    hold off;

end % end function
