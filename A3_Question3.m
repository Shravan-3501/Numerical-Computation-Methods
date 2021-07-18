%Interpretations: The actual value of the integral is 2*arctan(4) = 2.6516
% We observe that the trapezoidal rule and simpson's rule have a large
% error magnitude, whereas the composite trapezoidal's and composite
% simpson's rule are very close to the actual answer.
% The gauss quadrature method also has a much larger error magnitude as
% compared to the other methods above in part (c). This could be because we
% are using only 2 point gauss quadrature method. If we had used a larger
% value of n for the gauss quadrature method the result would have been
% closer to the actual value.

%Upper and lower limit of the integral are a and b respectively.
a = -4;
b = 4;
% f is the function whose integral is to be approximated.
f = @(x) (1/(1 + x*x));

res = trapezoidal(a,b,1,f);
disp("(a)Using trapezoidal rule:") 
disp(res);

res = simpsons(a,b,2,f);
disp("(b)Using simpson's rule:")
disp(res);

n = 10; 
res = trapezoidal(a,b,n,f);
disp("(c)(i)Using composite trapezoidal rule with n = " + n);
disp(res);

res = simpsons(a,b,10,f);
disp("(c)(ii)Using composite simpson's rule with n = " + n);
disp(res);

res = gauss_quad_2(f,a,b);
disp("(d)Using two point Gauss-Legrande Quadrature:");
disp(res);

disp("Interpretations: The actual value of the integral is 2*arctan(4) = 2.6516");
disp("We observe that the trapezoidal rule and simpson's rule have a large error magnitude, whereas the composite trapezoidal's and composite simpson's rule are very close to the actual answer.");
disp("The gauss quadrature method also has a much larger error magnitude as compared to the other methods above in part (c)");
disp("This could be because we are using only 2 point gauss quadrature method. If we had used a larger value of n for the gauss quadrature method the result would have been closer to the actual value.");


function res = trapezoidal(a,b,n,f)
    h = (b-a)/n; %Step size in terms of end points and n
    res = (f(a)+f(b))/2;
    x = a;
    for i = 1:(n-1)
        x = x + h;
        res = res + f(x); % We keep adding approximate value of the function at the intermediate points
    end
    res = res*h; 
end

function res = simpsons(a,b,n,f)
    h = (b-a)/n; %Step size in terms of end points and n
    res = f(a) + f(b);
    x = a;
    for i = 1:(n-1)
        x = x + h;
        %The formula for simpson's rule with n ordinates
        if(mod(i,2) == 0)
            res = res + 2*f(x); %Alternate terms have co-efficient 2 
        else
            res = res + 4*f(x); %Alternate terms have co-efficient 4
        end
    end
    res = res*(h/3);
end

function res = gauss_quad_2(f,a,b)
    %We directly use the formula derived for 2 point gauss quadrature as
    %derived in the lecture.
    g = @(x) (f(((b-a)*(x/2)) + (a+b)/2));
    res = ((b-a)/2)*(g(1/(3^0.5)) + g(-1/(3^0.5)));
end