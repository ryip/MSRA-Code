function [ DRvals ] = DRdiffs( X )
% Input: matrix of values
% Output: matrix of values minus their down-right neighbor. Zeros along the
% right column and bottom row

[numrows, numcols] = size(X);

DRvals = zeros(numrows, numcols);
for i = 1:numrows - 1;
    for j = 1:numcols - 1;
        DRvals(i,j) = X(i+1,j+1) - X (i,j);
    end
end

end

