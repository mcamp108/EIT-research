if ~exist('did_startup')
    cd 'C:\Users\Mark\Documents\GraduateStudies\LAB\EIDORS\eidors';
    run 'startup.m';
    did_startup= 1;
end % end if

[fmdl, new_electrode_centers]= mk_humHead_fmdl(numelec, elec_z_plane, 'PigHeadMesh2refine.vol');
img= mk_image(fmdl, 1);