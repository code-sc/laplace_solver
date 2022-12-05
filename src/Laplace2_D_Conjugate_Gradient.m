step = 1;           %Step size in x and y

x_start = 1;        %Starting point of x    (x1)
x_stop = 100;       %End point of x         (x2)     
x = x_start:step:x_stop;    %Datapoints
lx = length(x);

y_start = 1;        %Starting point of y    (y1)
y_stop = 100;       %End point of y         (y2)
y = y_start:step:y_stop;    %Datapoints
ly = length(y);


b = zeros(ly-2,lx-2);   %Datapoints initialised to 0
Max_run = 1000000;  %Maximum number of iterations
counter = 0;        %Loop counter

%Boundary Condition
b(1,:) = b(1,:)+zeros(1,lx-2);       %Boundary 1 [(x1,y2)->(x2,y2)]
b(:,1) = b(:,1)+zeros(ly-2,1);       %Boundary 2 [(x1,y2)->(x1,y1)]
b(lx-2,:) = b(lx-2,:)+zeros(1,lx-2);      %Boundary 3 [(x1,y1)->(x2,y1)]
b(:,ly-2) = b(:,ly-2)+10.*ones(ly-2,1);   %Boundary 4 [(x2,y1)->(x2,y2)]
b = b(:);            %Vectorising the matrix of boundary condition

lb = length(b);
one_s = -1*ones(lb,1);
temp = spdiags([one_s,-4*one_s,one_s],[-1,0,1],lx-2,lx-2);
A = blkdiag(temp);      %Developing the A matrix in terms of block matrices
for a=1:1:lx-3
    A = blkdiag(A,temp);
end
A = A + spdiags([one_s,one_s],[-(lx-2),(ly-2)],lb,lb);
res = CG(A,b);          %Filling in the rest of the elements
res = reshape(res,[lx-2,ly-2]);

res = [zeros(1,lx-2);res;zeros(1,lx-2)];    %Resustituting the boundary conditions for plotting
res = [zeros(ly,1) res 10.*ones(ly,1)];

[X,Y] = meshgrid(x,y);  %Plotting
subplot(2,1,1);
surf(X,Y,res);
subplot(2,1,2);
contour(res);    