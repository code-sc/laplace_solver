%Solution to 2D laplace equation by relaxation method
%
%    (1,1)                (1,lx)
%   (x1,y2)--------------(x2,y2)
%      |                    |
%      |                    |
%      |                    |
%      |                    |
%   (x1,y1)--------------(x2,y1)
%   (ly,1)               (lx,ly)
%
step = 1;           %Step size in x and y

x_start = 1;        %Starting point of x    (x1)
x_stop = 100;       %End point of x         (x2)  
x = x_start:step:x_stop;    %Datapoints
lx = length(x);     %Length of x

y_start = 1;        %Starting point of y    (y1)
y_stop = 100;       %End point of y         (y2)
y = y_start:step:y_stop;    %Datapoints
ly = length(y);     %Length of y

A = zeros(ly,lx);   %Datapoints initialised to 0
Max_run = 1000000;  %Maximum number of iterations
counter = 0;        %Loop counter

%Boundary Condition
A(1,:) = zeros(1,lx);       %Boundary 1 [(x1,y2)->(x2,y2)]
A(:,1) = zeros(ly,1);       %Boundary 2 [(x1,y2)->(x1,y1)]
A(lx,:) = zeros(1,lx);      %Boundary 3 [(x1,y1)->(x2,y1)]
A(:,ly) = 10.*ones(ly,1);   %Boundary 4 [(x2,y1)->(x2,y2)]

while counter<Max_run
    temp = A(2:lx-1,2:ly-1);    %Creating a copy to compare
    if rem(counter,2) == 0
        for a = 2:2:lx-1        %The black boxes
            for b = 2:2:ly-1
                A(a,b) = (A(a+1,b)+A(a,b+1)+A(a-1,b)+A(a,b-1))/4;
            end
        end
        for a = 3:2:lx-1
            for b = 3:2:lx-1
                A(a,b) = (A(a+1,b)+A(a,b+1)+A(a-1,b)+A(a,b-1))/4;
            end
        end
    else
        for a = 2:2:lx-1        %The white boxes    
            for b = 3:2:ly-1
                A(a,b) = (A(a+1,b)+A(a,b+1)+A(a-1,b)+A(a,b-1))/4;
            end
        end
        for a = 3:2:lx-1
            for b = 2:2:lx-1
                A(a,b) = (A(a+1,b)+A(a,b+1)+A(a-1,b)+A(a,b-1))/4;
            end
        end
    end
    if max(abs(temp - A(2:lx-1,2:ly-1))) < 1e-6 %Accuracy check
        break
    end
    counter = counter + 1;
end

disp(counter/2);

[X,Y] = meshgrid(x,y);
subplot(2,1,1);
surf(X,Y,A);
subplot(2,1,2);
contour(A);