function out_img= testStimPattern(img, stim, cutPlanes)
% This function was written by :
%                               Mark Campbell
%                               Carleton University
% Show surface voltages of img
out_img= [];
for i= stim
    [img.fwd_model.stimulation, img.fwd_model.meas_select] = mk_stim_patterns(32,1,[0,i],[0,i]);
    img.fwd_model.normalize_measurements = 0;
    testImg= img;
    testImg.fwd_solve.get_all_meas = 1;
    vh=fwd_solve(testImg);
    imgv = rmfield(testImg,'elem_data');
    stimSlices= [];
    for j= 1:32
        imgv.node_data = vh.volt(:,j);
        colours = calc_colours(imgv,[]);
        % Testing depth voltages    
        imgv.calc_colours.backgnd = 0.9*[1,1,1];
        imgv.show_slices.img_cols= 1;
        if j== 1
            stimSlices= show_slices(imgv,cutPlanes.*[1,inf,inf]);
        else
            stimSlices= stimSlices+ show_slices(imgv,cutPlanes.*[1,inf,inf]);
        end % end if
    end % end for
    % Concatenate into matrix showing results from all stim patterns
    if i== 1
        out_img= stimSlices;
    else
        out_img= [out_img, stimSlices];
    end % end if
end % end for
% Process images
foreground= out_img(out_img~=min(out_img));
background= out_img(out_img==min(out_img));
foreground= foreground- min(foreground);
out_img(out_img~=min(out_img))= foreground;
out_img(background)= 1;

end % end function