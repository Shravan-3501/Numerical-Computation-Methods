%p,q,f defined for -u''+p(x)u'+q(x)u = f(x)
p = @(x) (2*x*x - 3*x);
q = @(x) (x);
f = @(x) (5*x + 1);
% The two end points a and b defined below with the initial condition of
% the functions value at the two end points i.e. u(a) = g0 and u(b) = g1
a = 0;
b = 1;
g0 = 0;
g1 = 1;
n = 10; %The number for intervals we wish to divide the interval [a,b] into
disp("Using n = " + n + " the value of the function at the intermediate values is :");
disp(BVP(p,q,f,a,b,g0,g1,n));

function u = BVP(p,q,f,a,b,g0,g1,n)
    h = (b-a)/n; %Step size defined in terms of end points and number of intervals.
    X = []; % X is the vector which contains x values in the interval [a,b] for which we are calculating the functions approximate value
    temp = a;
    while(temp<=b)
        temp = temp+h;
        X = [X temp];
    end
    A = zeros(n-1,n-1);
    B = zeros(n-1,1);
    % We use the formulae to create matrices A and b to solve the (n-1)
    % equations Ax = b
    A(1,1) = 2 + q(X(1))*h*h;
    A(1,2) = (1/2)*p(X(1))*h - 1;
    B(1,1) = (1 + (1/2)*p(X(1))*h)*g0 + (h*h*f(X(1)));
    for i = 2:(n-2)
        A(i,i-1) = -1*(1 + (1/2)*p(X(i))*h);
        A(i,i) = 2 + q(X(i))*h*h;
        A(i,i+1) = (1/2)*p(X(i))*h - 1;
        B(i,1) = h*h*f(X(i));
    end
    A(n-1,n-2) = -1*(1 + (1/2)*p(X(n-1))*h);
    A(n-1,n-1) = 2 + q(X(n-1))*h*h;
    B(n-1,1) = h*h*f(X(n-1)) + (g1*(1 - (1/2)*p(X(n-1))*h));
    guess = zeros(n-1,1); % Initial guess is 0 for all variables
    u = GaussSeidel(A,B,guess); % We use the gauss-seidel method to solve the equations
end

function x = GaussSeidel(A,b,guess)
    epsilon = 10^-10;  % error tolerance
    n = size(A,1);
    x = guess;
    for i = 1:n
        x(i) = b(i);
        for j = 1:n
            if(j ~= i)
                x(i) = x(i) - A(i,j)*x(j); %Iterative formula for gauss-seidel method
            end
        end
        x(i) = x(i)/A(i,i);
    end
    while(error(guess,x) > epsilon)
        guess = x;
        for i = 1:n
            x(i) = b(i);
            for j = 1:n
                if(j ~= i)
                    x(i) = x(i) - A(i,j)*x(j); %Iterative formula for gauss-seidel method
                end
            end
            x(i) = x(i)/A(i,i);
        end
    end
end

%Function to calculate the error between two consecutive guesses using L-infinity norm
function e = error(guess,x)
    e = 0;
    n = size(x,1);
    for i = 1:n
        e = max(e,abs(guess(i) - x(i)));  %L-infinity norm
    end
end