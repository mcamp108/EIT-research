function testAddInclusion(imdl, stim, meas, loc)
% This function was written by :
%                               Mark Campbell
%                               Carleton University
% Parameters:
%   loc: a 1x4 array [x y z r] where x y and z are the coordinates specifying
%   the center of the inclusion, and r is the radius of the inclusion.
if length(loc)<4
    disp("loc must be a 1x4 array in the form [x y z r]");
elseif loc(4)<= 0
    disp("The inclusion radius must be positive.");
end % end if

% Set stim pattern
[imdl.fwd_model.stimulation, imdl.fwd_model.meas_select] = mk_stim_patterns(32,1,[0,stim],[0,meas]);
img1= mk_image(imdl);
% Add inclusion in center of model with diameter = 12 (~10% head radius)
select_fcn = inline('(x- loc(1)).^2 + (y- loc(2)).^2 + (z- loc(3)).^2 < loc(4)^2','x','y','z');
memb_frac = elem_select( img1.fwd_model, select_fcn);
img2= mk_image(img1, 1+ memb_frac); % not sure if this is correct
img2.elem_data = img1.elem_data + memb_frac; % or if this is correct
img2.calc_colours.cb_shrink_move = [0.3,0.6,0.02];
% Visualize true inclusion
show_fem(img2,1); axis tight;
show_3d_slices(img2, [loc(3)- 1, loc(3)+ 1], [loc(1)- 1], [loc(2)+ 1]);
% Simulate voltages
vh= fwd_solve(img1);
vi= fwd_solve(img2);
% Reconstruction and visualize
J = calc_jacobian( calc_jacobian_bkgnd( imdl) );
iRtR = inv(prior_noser( imdl ));
hp = 0.17;
iRN = hp^2 * speye(size(J,1));
RM = iRtR*J'/(J*iRtR*J' + iRN);
imdl.solve = @solve_use_matrix; 
imdl.solve_use_matrix.RM  = RM;
imgr = inv_solve(imdl, vh, vi);
imgr.calc_colours.ref_level = 0; % difference imaging
imgr.calc_colours.greylev = -0.05;
show_fem(imgr);
print_convert('basic_3d_06a.png','-density 60');
show_3d_slices(imgr, [1,1.9], [0.5],[0.5]);
view(-14,13); axis tight; axis equal;
print_convert('basic_3d_06b.png','-density 60');

end % end function