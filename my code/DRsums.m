function [ summatrix ] = DRsums( a )
%INPUT: a matrix
%OUTPUT: a matrix with the sum of the elements on each DR diagonal on the
%top and left entries of the matrix and zeros everywhere else.

[numrows, numcols] = size(a);
summatrix = a;

%we split into two parts: the top row, and the left column

for j = 1:numcols
    rownum = 1;
    colnum = j;
    while (rownum < numrows) && (colnum < numcols)
        rownum = rownum + 1;
        colnum = colnum + 1;
        summatrix(1,j) = summatrix(1,j) + summatrix(rownum,colnum);
        summatrix(rownum,colnum) = 0;
    end
end

for i = 2:numrows
    rownum = i;
    colnum = 1;
        while (rownum < numrows) && (colnum < numcols)
        rownum = rownum + 1;
        colnum = colnum + 1;
        summatrix(i,1) = summatrix(i,1) + summatrix(rownum,colnum);
        summatrix(rownum,colnum) = 0;
        end
end


end

