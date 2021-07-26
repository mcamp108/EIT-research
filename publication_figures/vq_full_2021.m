% open images

% how many < 0.1, 0.1 to 0.8, 0.8 to 1.2, > 1.2
fname = 'C:\Users\Mark\Documents\EIT\proj\vq_apnea_breathing\VQ_images_ap_br.xlsx';
saveDir = 'C:\Users\Mark\Documents\EIT\proj\vq_apnea_breathing\figures\results\';
Excel = actxserver('Excel.Application'); 
set(Excel, 'Visible', 0);
set(Excel,'DisplayAlerts',0);
Workbooks = Excel.Workbooks;
[type, sheet_names] = xlsfinfo(fname);
Workbook = Workbooks.Open(fname);
Sheets = Excel.ActiveWorkBook.Sheets;

opts = detectImportOptions(fname);
opts.DataRange = 'A1';
results = struct;
for i = 1:length(sheet_names)
    sn = sheet_names{i};
    A = readmatrix(fname, 'Sheet', sn);
    v1 = sum(A < 0.1, 'all');
    v2 = sum(A >= 0.1 & A < 0.8, 'all');
    v3 = sum(A >= 0.8 & A < 1.2, 'all');
    v4 = sum(A > 1.2, 'all');
    fnL = strsplit(sn, ' ');
    fn = sprintf('%s_%s_%s', fnL{2}, fnL{1}, fnL{3});
    results.(fn) = [v1, v2, v3, v4];
end

Workbook.Save();
Workbook.Close();
Workbooks.Close;
Excel.Quit();
delete(Excel);

fn = fieldnames(results);
horses = {'296439', '296695', '296954', '296955'};

final = struct;
for i=1:length(horses)
    horse = horses{i};
    hmat = zeros(4,4);
    for j = 1:length(fn)
        thisfn = fn{j};
        if contains(thisfn, horse)
            res = results.(thisfn);
            if contains(thisfn, 'apnea') && contains(thisfn, 'T25')
                row = 1;
            elseif contains(thisfn, 'breathing') && contains(thisfn, 'T25')
                row = 2;
            elseif contains(thisfn, 'apnea') && contains(thisfn, 'T50')
                row = 3;
            elseif contains(thisfn, 'breathing') && contains(thisfn, 'T50')
                row = 4;
            else
                disp('error');
            end
            hmat(row, :) = (res / sum(res)) * 100;
        end % end if
    end % end for j
    final.(sprintf('h_%s', horse)) = hmat;
    fg = figure();
    sgtitle(horse);
    bar(hmat, 'stacked');
    disp(horse);
    disp(hmat);
    ylabel('Percentage of pixels');
    ylim([0,100]);
    xlabel('Pixel values');
    ax = gca;
    ax.XTickLabel = {'T25 apnea', 'T25 breathing', 'T50 apnea', 'T50 breathing'};
    legend({'$\dot{V}/Q < 0.1$'; '0.1 <= $\dot{V}/Q < 0.8$'; '$0.8 <= \dot{V}/Q < 1.2$'; '$\dot{V}/Q > 1.2$'}, 'Interpreter', 'latex', 'location', 'bestoutside');
%     ax.XTickLabel = {'VQ < 0.1'; '0.1 <= VQ < 0.8'; '0.8 <= VQ < 1.2'; 'VQ > 1.2'};
%     legend({'IE11 apnea', 'IE11 breathing', 'IE13 apnea', 'IE13 breathing'}, 'location', 'northwest');
    ax.XTickLabelRotation = 15;

    saveas(fg, sprintf('%s%s hist.png', saveDir, horse));
end

    