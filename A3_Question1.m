f = @(x,y) ((y/x) - ((y/x)^2)); % f(x) defined where y' = f(x)
x0 = 1; % Initial conditions
y0 = 1;
h = 0.01; %Step size
x1 = 2; %Value for which the functions value is to be estimated
[y,y_out] = RangeKutta(f,x0,y0,h,x_target); % The output y is the estimation for f(x1) and the vector y_out is the estimation across the range of values from x0 to x1 
x_lin = linspace(x0,x1,(x1-x0)/h + 1); %Points on the x-axis for which we will plot the function's approximation
disp("The approximate value of the function at x = " + x1 + " is:");
disp(y);
plot(x_lin,y_out); % We plot our approximation to the function in [x0,x1]

function [y,y_out] = RangeKutta(f,x0,y0,h,x_target)
    n = (x_target - x0)/h;  %number of points in the range
    y_out = [y0];
    y = y0;
    for i = 1:n
        x = x0 + (i-1)*h; %Current x co-ordinate for which we are calcualating the approximate value of the funtion
        %The four terms to be calculated in the range-kutta method
        k1 = h*f(x,y);
        k2 = h*f(x + (h/2),y + (k1/2));
        k3 = h*f(x + (h/2),y + (k2/2));
        k4 = h*f(x + h,y + k3);
        y = y + (1/6)*(k1 + 2*k2 + 2*k3 + k4);
        y_out = [y_out y]; % We append our approximation to functional value to the vector y_out
    end
end