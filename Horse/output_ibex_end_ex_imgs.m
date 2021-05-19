path = 'C:\Users\Mark\Documents\EIT\Will_Ditcham\22-04-2021_#2\';

files = ls(path);
files = files(3: end, :);

lungs = human_roi();

for i = 1: size(files, 1)
    f = strtrim(files(i, :));
    if contains(f, '.xls') && ~strcmp(f(1:2), 'TI') && ~strcmp(f(1:4), 'EELI')
        ff = horzcat(path, f);
        
        [genNUM, genTXT, genRAW] = xlsread(ff, 'general info');
        fsIdx = find_column(genRAW(:, 1), 'Sampling rate [Hz]');
        fs = genRAW(fsIdx, 2);
        fs = fs{1};
        
        [NUM, TXT, RAW] = xlsread(ff, 'filtered results');
        TXT = TXT(2, :);
        idx = find_column(TXT, 'Start of Inspiration [s]');
        insStartIdx = round(NUM(2: end, idx) .* fs);
        matFile = find_mat_file(ff);
        save_end_ex_imgs(matFile, insStartIdx, lungs);
    end
end
%%
% data = load(matFile);
% data = data.data;
% img = data.measurement.ZeroRef;
% %%
% fResults = xlsread(ff, 'filtered results');
% fs = 50.86;
% trel = 0;
% insStartIdx = round(fResults(2:end, 36) .* fs) + trel;
% insEndIdx = round(fResults(2:end, 37) .* fs) + trel;
% nBreaths = length(insStartIdx);
% 
% delta = insEndIdx(1) - insStartIdx(1);
% 
% % testVal = 31.8462184531766;
% res = nan(200, 1);
% for i=1:2000
%     tmp = img(:, :, i + delta) - img(:, :, i);
%     res(i) = sum(tmp, 'all') - sum(testImg, 'all');
%     if res(i) == 0
%         disp(i);
%     end
% end

% 
% %%
% 
% outImgs = zeros(size(lungs,1), size(lungs,2), nBreaths);
% % insStartIdx = 50:60;
% for i = 1 : nBreaths
%      thisImg = img(:, :, insStartIdx(i));
%      thisImg(thisImg < 0) = 0;
%      outImgs(:, :, i) = thisImg .* lungs;
% end
% 
% thisImg = img(:, :, insStartIdx(1));
% thisImg = thisImg(:);
% notLung = find(lungs == 0);
% 
% % for i = 1:size(img, 3)
% for i = 40: 60
%     thisImg = img(:, :, i);
% %     thisImg (thisImg < 0) = 0;
%     thisLung = thisImg .* lungs;
%     oldSum = sum(thisLung, 'all');
%     
%     temp = thisLung(:);
%     for j = 1:length(notLung)
%         newPixel = thisImg(notLung(j));
%         newSum = round(newPixel + oldSum);
%         if newSum == 160026
%             if notLung(j) > 0
%                 fprintf('i = %.0f j = %.0f\n', i, notLung(j) );
%                 temp2 = lungs;
%                 temp2(notLung(j)) = 2;
%                 figure();
%                 imagesc(temp2);
%             end
%         end
%     end
% end
% %%
% t1 = insStartIdx(1);
% t2 = round(2.47738888100515 * fs);
% im1 = img(:, :, t1);
% im2 = img(:, :, t2);
% tv = img(:, :, t2) - img(:, :, t1);
% tv .* lungs




function idx = find_column(txt, colName)
    for idx = 1:length(txt)
        if strcmp(txt{idx}, colName)
            return
        end
    end
    idx = nan;
end


function mf = find_mat_file(f)
    loc = regexp(f, '.mat');
    stop = loc + 3;
    mf = f(1: stop);
end


function img = load_img_from_mat(matFile)
    data = load(matFile);
    img = data.data.measurement.ZeroRef;
end


function save_end_ex_imgs(matFile, insStartIdx, roi)
    Excel = actxserver('Excel.Application'); 
    set(Excel, 'Visible', 0);
    set(Excel,'DisplayAlerts',0);
    Workbooks = Excel.Workbooks; 
    
    nBreaths = length(insStartIdx);
    img = load_img_from_mat(matFile);
    temp = strsplit(matFile, '\');
    fName = temp(end);
    fName = fName{1};
    eeliName = horzcat('EELI_', fName(1: end-4), '.xlsx');
    outName = strrep(matFile, fName, eeliName);
    for i = 1 : nBreaths
         thisImg = img(:, :, insStartIdx(i)) .* roi;
         thisImg(thisImg < 0) = 0;
         writematrix(thisImg, outName, 'Sheet', num2str(i));
    end
    
    % delete the first 3 empty sheets
    [type, sheet_names] = xlsfinfo(outName);
    Workbook = Workbooks.Open(outName);
    Sheets = Excel.ActiveWorkBook.Sheets;
    index_adjust = 0;
    for i = 1:max(size(sheet_names))
        current_sheet = get(Sheets, 'Item', (i - index_adjust));
        if contains(sheet_names{i}, 'Sheet')
            invoke(current_sheet, 'Delete');
            index_adjust = index_adjust + 1;
        else
            NewRange = current_sheet.Range('A1:AF1');
            NewRange.ColumnWidth = 8;
        end 
    end
    Workbook.Save();
    Workbook.Close();
    Workbooks.Close;
    Excel.Quit();
    delete(Excel);
end