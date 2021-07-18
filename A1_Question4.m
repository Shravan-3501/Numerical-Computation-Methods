f = @(x) (tanh(x));  %The function f whose root is to be found is defined here.
f1 = @(x) (sech(x)^2);  % We manually calculate the differention of f and store it as a function here. 
s = 0.1;  
a = -10;  %Initial limits of the interval defined here, where a is the lower limit
b = 15;  % b is the upper limit.

Root = bisection(f,f1,a,b,s);  %The function returns the root of the equation calculated partially using bisection method and partially using Newton-Raphson Method.
disp("Root of the equation is: " + Root); % Printing the root.

function root = bisection(f,f1,a,b,s)
    stop = s*(b-a);   % The stopping criterion for the bisection method.
    while((b-a)>=stop)
        c = (a+b)/2;    
        if(f(a)*f(c) < 0)  %This condition implies that the root is present in the lower half. 
            b = c;     %Hence, we update the upper limit of the interval.
        else    % This implies that the root is present in the upper half of the interval. 
            a = c;  % Hence, we update the lower limit of the interval.
        end
    end
    root = newton(f,f1,c);  %After this, we proceed to find the root using Newton-Raphson method.
end

function guess = newton(f,f1,c)
    epsilon = 10^-10;  %The value of epsilon for the stopping criterion of the newton method.
    guess = c;  
    while(abs(f(guess)) > epsilon) %Stopping criterion for the method.
        difference = f(guess)/f1(guess);
        guess = guess - difference;  % Update step for the newton-raphson method.
    end
end