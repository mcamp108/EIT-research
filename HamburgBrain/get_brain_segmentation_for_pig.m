function brain_seg= get_brain_segmentation_for_pig(seq)

% -------------------------------------------------------------------------
% DESCRIPTION:
%
%   brain_seg= get_brain_segmentation_for_pig(seq)
%
% -------------------------------------------------------------------------
% PARAMETERS:
% 
%   seq:        
%       A swine sequence struct from load_HamburgBrain_data.
%
% -------------------------------------------------------------------------   
% RETURNS:
% 
%   brain_seg: 
%       A struct containing fields:
%
%       upperL: 
%           The upper left reference corner for segmentatin of the brain
%           from reconstructed 64 x 64 image of airchill pig.
%           
%       mtx: 
%           The logical matrix representing the pixels in the 64 x 64 image
%           that are included in the brain, referenced from upperL.
%
%       idx:
%           mtx indices converted to pixel indices of imgr.elem_data cell.
%
%       topL:
%           Indices in imgr.elem_data of the top left portion of the brain.
%
%       topR:
%           Indices in imgr.elem_data of the top right portion of the
%           brain.
%
%       botL
%           Indices in imgr.elem_data of the bottom left portion of the
%           brain.
%
%       botR:
%           Indices in imgr.elem_data of the bottom right portion of the
%           brain.
% -------------------------------------------------------------------------   
% AUTHOR:
%   Mark Campbell
%   Carleton University
%   markacampbell@cmail.carleton.ca
% -------------------------------------------------------------------------

brain_seg= struct;
brain_seg.upperL= [9, 21]; % row, column
brain_seg.mtx= [ % segmentatino for 10.2 (used for all pigs for now)
                        0 0 0 0 0 0 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0; % 9
                        0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0; % 10
                        0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0; % 11
                        0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0; % 12
                        0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0; % 13
                        0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0; % 14
                        0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0; % 15
                        1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0; % 16
                        1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1; % 17
                        1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1; % 18
                        1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1; % 19
                        1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1; % 20
                        0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0; % 21
                        0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0; % 22
                        0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0; % 23
                        0 0 0 0 0 0 0 1 1 1 1 1 1 1 0 0 0 0 0 0 0; % 24
                        0 0 0 0 0 0 0 1 1 1 1 1 0 0 0 0 0 0 0 0 0; % 25
                        0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0; % 26
                                                                  ];
brain_seg.idx= convert_matrix_to_elem_data_index(seq.imgr, brain_seg.mtx, brain_seg.upperL);

[hei, wid]= size(brain_seg.mtx);
tbMid= round(hei/2);
lrMid= floor(wid/2);

top= brain_seg.mtx;
top(tbMid+1: end, :)= 0;

bot= brain_seg.mtx;
bot(1:tbMid, :)= 0;

topL= top;
topR= top;
botL= bot;
botR= bot;

topL(:, lrMid+1: end)= 0;
topR(:, 1:lrMid)= 0;
botL(:, lrMid+1:end)= 0;
botR(:, 1:lrMid)= 0;

brain_seg.topL= convert_matrix_to_elem_data_index(seq.imgr, topL, brain_seg.upperL);
brain_seg.topR= convert_matrix_to_elem_data_index(seq.imgr, topR, brain_seg.upperL);
brain_seg.botL= convert_matrix_to_elem_data_index(seq.imgr, botL, brain_seg.upperL);
brain_seg.botR= convert_matrix_to_elem_data_index(seq.imgr, botR, brain_seg.upperL);

end % end function

function idx= convert_matrix_to_elem_data_index(imgr, mtx, upperL)
    imgr.elem_data= imgr.elem_data(:, 1);
    convertIdx= calc_slices(imgr); % convert 64x64 index to elem_data index
    sz_convertIdx= size(convertIdx);
    notNan= isnan(convertIdx(:, :))==0;
    notNan= notNan* 1; % convert to cell
    idx_ref= (upperL(2)-1)* size(convertIdx, 1) + upperL(1)-1;
    for col= 1: size(mtx, 2)
        ispart= find(mtx(:, col))+ sz_convertIdx(1)*(col-1)+ idx_ref;
        notNan(ispart)= 2; % index as part of 64x64 image
    end % end for
    % imgr.elem_data indexing goes from bottom up, left to right. Indexing in
    % notNan goes from top down, left to right. Transform notNan to emel_data
    % by transpose then fliplr
    notNan= fliplr(notNan');
    convertIdx= notNan(notNan~=0);
    idx= find(convertIdx==2);
end % end function