extra={'out','in', [ ...
       'solid cut = plane(0,0,1.1;0,0,1);' ...
       'solid in  = sphere(0,0,1;0.5);' ...
       'solid out = sphere(0,0,1;0.65) and (not in);' ...
       ]};
fmdl= ng_mk_cyl_models(2,[32,1],[0.05],extra);
% fmdl= ng_mk_cyl_models([2,1,0.1],[16,1],[0.05],extra); % fine mesh
stim = mk_stim_patterns(32,1,[0,5],[0,5],{'no_meas_current_next2'},1);
fmdl.stimulation = stim;
img = mk_image(fmdl,0.41); % scalp
img.elem_data(fmdl.mat_idx{2}) = 0.016; % skull
img.elem_data(fmdl.mat_idx{3}) = 0.47; % brain


fmdl2 = fix_internal_mat_idx(fmdl, 2);

% look at points in each mat_idx

m1 = fmdl.nodes( fmdl.elems( fmdl.mat_idx{1},: ), : );
m2 = fmdl.nodes( fmdl.elems( fmdl.mat_idx{2},: ), : );
m3 = fmdl.nodes( fmdl.elems( fmdl.mat_idx{3},: ), : );

% try to assign elems in m3 (brain) to m2 (skull)
wireframe(fmdl, 2);
wireframe(fmdl, 3);
scat3(fmdl, 1);


function fmdl = fix_internal_mat_idx(fmdl, targetMat, replaceMat)
% targetMat: 
%       material for with any internal elements not assigned to the material will be reassigned to this material
% replaceMat: 
%       optional material that will only be reassigned to targetMat
bnd = find_boundary( fmdl.nodes(fmdl.elems( fmdl.mat_idx{targetMat},: ), :) );
for i=1:size(fmdl.elems,1)
    fcs = find_faces( fmdl.elems(i,:) );
    for j=1:4
        triNodes = fcs(j, :);
        has = (fmdl.elems == triNodes(1)) + (fmdl.elems == triNodes(2)) + (fmdl.elems == triNodes(3));
        inElem = find( sum(has,2) == 3 ); % find the two tets this face is a part of
        isMat = sum( fmdl.mat_idx{targetMat} == inElem' );
        if (sum(isMat) == 1) && (length(isMat) ==2)% one is in mat one is not in mat
            assign = inElem(isMat==0);
            % remove assign from starting material
            if nargin == 3
                isInReplace = (fmdl.mat_idx{replaceMat} == assign);
                if sum(isInReplace) == 1
                    fmdl.mat_idx{replaceMat} = fmdl.mat_idx{replaceMat}( ~(fmdl.mat_idx{replaceMat} == assign) );
                    fmdl.mat_idx{targetMat} = [fmdl.mat_idx{targetMat}; assign];
                end
            else
                for m=1:length(fmdl.mat_idx)
                    if m == targetMat
                        continue
                    else
                        fmdl.mat_idx{m} = fmdl.mat_idx{m}( ~(fmdl.mat_idx{m} == assign) );
                    end % end if
                end % end for m
                % assign assign to target material
                fmdl.mat_idx{targetMat} = [fmdl.mat_idx{targetMat}; assign];
            end % end if

        end % end if
    end % end for j
end % end for i

end % end for

function wireframe(fmdl, material)
nodeIdx = fmdl.elems( fmdl.mat_idx{material},: );
hold on;
h = trimesh(nodeIdx, fmdl.nodes(:,1),fmdl.nodes(:,2),fmdl.nodes(:,3));
h.FaceAlpha = 0;
hold off;
end % end function

function scat3(fmdl, mat)
nodes = fmdl.nodes( fmdl.elems( fmdl.mat_idx{mat},: ), : );
hold on;
scatter3(nodes(:,1), nodes(:,2), nodes(:,3));
hold off;
end % end function


function faces = find_faces(tetNodes)
faces = zeros(4,3);
faces(1,:) = tetNodes(1,[1,2,3]);
faces(2,:) = tetNodes(1,[1,4,2]);
faces(3,:) = tetNodes(1,[1,3,4]);
faces(4,:) = tetNodes(1,[2,4,3]);
end % end function

