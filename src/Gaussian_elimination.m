function [res] = Gaussian_elimination(A,b)
%Gaussian elimination using partial pivoting and row IDs to pivot instead
%of actual row swaps.
if nargin ~= 2
    error ("2 inputs required");
elseif min(size(b))~=1
    error ("b must be a vector");
else
    if iscolumn(b)==0 
        b = b';
    end
    [r,c] = size(A);
    if r~=c
        error("Not a square matrix");
    end
    c = 1:r;            %Permutation vector (will have the actual vector IDs)
    res = zeros(r,1);   %Solution vector initialization
    for m=1:1:r-1       %Going along the column
        max_pos = m;         %Stores which element of the col has max value
        for n = m+1:1:r
            if abs(A(c(n),m))>abs(A(c(max_pos),m))
                max_pos = n;
            end
        end                 %Will have pos of largest element of column m
        if A(c(max_pos),m)>A(c(m),m)
            temp = c(max_pos);  %Swapping of the row IDs of the 
            c(max_pos) = c(m);  %Largest element row and current row
            c(m) = temp;
        end
        %Creating the lower triangular matrix by subtracting linear
        %combination of rows
        for n = m+1:1:r
            b(c(n)) = b(c(n)) - b(c(m))*(A(c(n),m)/A(c(m),m));
            A(c(n),:) = A(c(n),:) - A(c(m),:).*(A(c(n),m)/A(c(m),m));
        end
    end
    %Now A is an upper triangular matrix (Using c as index).
    %Solving using back substitution
    for m = r:-1:1
        temp = 0;       %Term to be subtracted
        for n = m+1:1:r
            temp = temp+res(n)*A(c(m),n);
        end
        res(m) = (b(c(m)) - temp)/A(c(m),m);
    end    
end
end

