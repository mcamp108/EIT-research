run 'myStartup.m';
% Load model
maxsz= 0.2; maxh= 2; imgsize= [64 64];
[fmdl, img, imdl]= mk_humHead_fmdl(maxsz, maxh, imgsize);
[sFmdl, sImg]= mkSphereTestMesh(maxsz, maxh);

%%
scalpElems= fmdl.elems( fmdl.mat_idx{1}, :);
scalpElems= scalpElems(:);
scalpCoors= fmdl.nodes(scalpElems, :);
X= fmdl.nodes(:, 1);
Y= fmdl.nodes(:, 2);
Z= fmdl.nodes(:, 3);
idx1= X>0;
idx2= X<10;
idx= (idx1+idx2)==2;
scatter3(X(idx),Y(idx),Z(idx), '.');
axis equal
bounds= [fmdl.nodes];
save('humHeadExternal.txt', 'bounds', '-ascii');
%%
% Test stim patterns on spherical model of head
cutPlanes = [-50;-40;-30;-20;-10;0;10;20;30;40;50];
stim= (0:15);
out_img= testStimPattern(sImg, stim, cutPlanes);
cd 'C:\Users\Mark\Documents\GraduateStudies\LAB\Thesis\Figures';
printopts.resolution= 500;
imagesc(out_img);
axis equal
print_convert(char("sphereTestStimSkip0-14.png"), printopts);
keyboard;

% Testing addition of inclusion
testAddInclusion(imdl, 15, 15, [0 0 0 6]);
