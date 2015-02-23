function [ summatrix ] = DLsums( a )
%INPUT: a matrix
%OUTPUT: a matrix with the sum of the elements on each DR diagonal on the
%top and right entries of the matrix and zeros everywhere else.

[numrows, numcols] = size(a);
summatrix = a;

%we split into two parts: the top row, and the right column

for j = 1:numcols
    rownum = 1;
    colnum = j;
    while (rownum < numrows) && (colnum > 1)
        rownum = rownum + 1;
        colnum = colnum - 1;
        summatrix(1,j) = summatrix(1,j) + summatrix(rownum,colnum);
        summatrix(rownum,colnum) = 0;
    end
end

for i = 2:numrows
    rownum = i;
    colnum = numcols;
        while (rownum < numrows) && (colnum > 1)
        rownum = rownum + 1;
        colnum = colnum - 1;
        summatrix(i,numcols) = summatrix(i,numcols) + summatrix(rownum,colnum);
        summatrix(rownum,colnum) = 0;
        end
end


end

