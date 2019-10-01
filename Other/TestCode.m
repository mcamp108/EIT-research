
for L=[2,4];
    imdl = jan12_GREIT_cyl(L);
    % specifies that this is the inverse model
    img = mk_image(imdl); vh = fwd_solve(img);
    %create eidors image object, uses jacobian background for conductivity
    img.elem_data = 1+0.1*elem_select(img.fwd_model,'(x-0.5).^2+y.^2 +(z-2.1).^2<0.1^2'); vi=fwd_solve(img);
    subplot(2,2,L-1);
       show_3d_slices(img); view(10,10);
    subplot(2,2,L-0);
       show_3d_slices(inv_solve(imdl,vh,vi)); view(10,10);
end
