function [ repeatedmatrix ] = DRrepmat( a )
%INPUT: a matrix with desired entries on the top row and left column
%OUTPUT: a matrix with each diagonal equal to the representative desired
%entry

[numrows, numcols] = size(a);
repeatedmatrix = a;

%we split into two parts: the top row, and the left column

for j = 1:numcols
    rownum = 1;
    colnum = j;
    while (rownum < numrows) && (colnum < numcols)
        rownum = rownum + 1;
        colnum = colnum + 1;
        repeatedmatrix(rownum,colnum) = repeatedmatrix(1,j);
    end
end

for i = 2:numrows
    rownum = i;
    colnum = 1;
        while (rownum < numrows) && (colnum < numcols)
        rownum = rownum + 1;
        colnum = colnum + 1;
        repeatedmatrix(rownum,colnum) = repeatedmatrix(i,1);

        end
end


end

