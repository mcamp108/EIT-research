
function adjMat= getAdjMat(image)
% parameters: matrix dimensions m= int, num rows, n= int, num columns.
% made by mnmltype, retrieved from https://stackoverflow.com/questions/3277541/construct-adjacency-matrix-in-matlab
%0-0-0                  
%| | |                    
%0-0-0  4 connections 
%| | |
%0-0-0
% ** TODO:or 8 connectoins if diagonals are to be included*
% returns: an adjacency matrix assuming maximum connections= 4, whose values correspond to edge indices.
[m, n]= size(image);
I_size= m*n;
V = repmat([ones(m-1,1); 0],n, 1); % 1-off diagonal elements
V = V(1:end-1); % remove last zero
U = ones(m*(n-1), 1); % n-off diagonal elements
W = sparse(1:(I_size-1),    2:I_size, V, I_size, I_size)...
  + sparse(1:(I_size-m),(m+1):I_size, U, I_size, I_size);
nans= find(isnan(image));
try nans(1)
    W(nans, :)= 0;
    W(:, nans)= 0;
catch err
end
A = W + W'; % make W symmetric
[s,t,e] = find(triu(A));
m = numel(e);
e(:) = 1:m;
adjMat = sparse([s;t],[t;s],[e;e]);

end % end function