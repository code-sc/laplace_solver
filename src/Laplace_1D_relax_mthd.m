%Solution to 1D Laplace equation by relaxation method

x_start = 0;    %Starting point for x
x_stop = 10;    %Stopping point for x
x_step = 0.1;   %Step size
lx = (x_stop-x_start)/x_step + 1;   %Length of x
x = x_start:x_step:x_stop;          %Datapoints

Max_run = 1000000;                  %Maximum number of iterations
counter = 0;                        %Loop counter
y = zeros(lx,1);                    %Initialize to 0

%Boundary Condition
y(1) = 0;       %Boundary 1
y(lx) = 1;      %Boundary 2

while counter<Max_run
    temp = y(2:(lx-1),1);    %Creating a copy to compare
    for a = 2:1:(lx-1) %Relaxation method
        y(a) = (y(a-1) + y(a+1))/2;
    end
    if max(abs(temp - y(2:(lx-1),1))) < 1e-6
        break
    end
    counter = counter + 1;
end

disp(counter);
plot(x,y,'g.');
xlabel("x axis");
ylabel("f(x)");


