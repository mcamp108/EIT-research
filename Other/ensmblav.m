function [out_img, peaks_used]= ensmblav(imgs, references, options)
%    ENSEMBLE AVERAGING
% imgs is a cell of slices from calc_slices
% references is an array of time points indicating the beginning of an event of interest
% num is the number of repetitions to be shown in the output. Note that
% the number of events averaged for each occurence is length(references)/ num
% If using an img structure, use img.elem_data= ensmblav(img.elem_data, references, options)
if iscell(imgs)
    out_img= {};
    for i= 1:length(imgs)
        out_img{i}= ensmblav(imgs{i}, references, options);
    end
    return
% elseif isstruct(imgs)
%     data= imgs.elem_data;
else
    data= imgs;
end
dims= length(size(data)); % number of dimension 2 or 3
event_lengths= diff(references);
if dims== 2
    ts_length= size(imgs, 2);
else
    ts_length= size(imgs, 3);
end

if nargin >2
    if ~isfield(options, 'num')
        reps = 1;
    else
        reps= options.num;
    end
    
    if isfield(options, 'offset')
        references= references+ options.offset;
        references= references(references> 0);
    end
    
    if ~isfield(options, 'min_frames')
        frame_length= min(event_lengths);
        good_refs= references(1: length(references-1));
    else
        frame_length= options.min_frames;
        good_refs= [];
        for i= 1:length(event_lengths)
           if event_lengths(i)>= options.min_frames && references(i)+ frame_length <= ts_length
              good_refs= [good_refs, references(i)];
           end
        end
    end
    
end

out_img= [];
per_rep= floor(length(good_refs)/reps); % The number of references to include in each average
segment= [];
peaks_used= [];
i= 1;

for r= 1:reps
    new= 1;
    for seg= 1:per_rep
        if new == 1
            if dims== 3, segment= data(:, :, good_refs(i): good_refs(i)+frame_length);
            elseif dims== 2, segment= data(:, good_refs(i): good_refs(i)+frame_length);
            end
            new= 0;
            peaks_used_current_row= [];
        else
            if dims== 3, segment= segment+ data(:, :, good_refs(i): good_refs(i)+frame_length);
            elseif dims== 2, segment= segment+ data(:, good_refs(i): good_refs(i)+frame_length);
            end
        end
        peaks_used_current_row= [peaks_used_current_row, good_refs(i)];
        i= i+ 1;
    end % end for 
    segment= segment/ per_rep;
    peaks_used= [peaks_used; peaks_used_current_row];
    if dims== 3,out_img= cat(3, out_img, segment);
    elseif dims== 2, out_img= cat(2, out_img, segment);
    end

end % end for

end