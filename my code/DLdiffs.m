function [ DLvals ] = DLdiffs( X )
% Input: matrix of values
% Output: matrix of values minus their down-left neighbor. Zeros along the
% left column and bottom row

[numrows, numcols] = size(X);

DLvals = zeros(numrows, numcols);
for i = 1:numrows - 1;
    for j = numcols:-1:2;
        DLvals(i,j) = X(i+1,j-1) - X (i,j);
    end
end

end

