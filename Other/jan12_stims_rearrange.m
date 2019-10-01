function fmdl = jan12_stims_rearrange(fmdl,L)
    skip4 = {32,1,[0,5],[0,5],{'no_meas_current_next1'},1};
%  [~,fmdl] = elec_rearrange([16,2],'square', fmdl);
    idx = reshape(1:32,L,[])';
    %reshapes the 1x32 matrix into a 4x8 matrix (when L= 4)
    idx(2:2:end,:) = fliplr(idx(2:2:end,:));
    %extract even elements in rows(from 2 -> end, going by 2's) for all
    %rows. flip along vertical axis.
    fmdl.electrode(idx) = fmdl.electrode(:);
    %these electrodes become all electrodes in the model
    [fmdl.stimulation,fmdl.meas_select] = mk_stim_patterns(skip4{:});
