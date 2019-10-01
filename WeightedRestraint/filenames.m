function [ files ] = filenames()
% All the test files 
names= ls;
numFiles= size(names,1)- 2;
for f= 1:numFiles
    idx= f+2;
    files(f,:) = names(idx,:);           
end
