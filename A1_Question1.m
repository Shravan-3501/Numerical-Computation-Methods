% x = (1/a) if and only if (a*x) = 1 which is equivalent to finding the
% root of (a*x) - 1 = 0

% Hence we use newton's method to find the root of the equation (a*x - 1) = 0

a = 5;   % The value of a for which we will calculate the reciprocal.
f = @(x) (a*x - 1); %The function f is defined here.
f1 = @(x) (a); % We calculate and define the differentiation of the function f here.
guess = 0.1;  % The initial guess we give to the Newton Raphson's Algorithm.
Reciprocal = solve(f,f1,guess);  % The value returned by the function is stored in the variable named Reciprocal.
disp("Reciprocal of " + a + " = " + Reciprocal);    % Finally we print the output.

function guess = solve(f,f1,g)
    guess = g;
    epsilon = 10^-10;  % The value of epsilon to be used as the stopping criterion.
    while(abs(f(guess)) > epsilon)  % Condition to stop the algortithm.
       difference = f(guess)/f1(guess);
       guess = guess - difference;      % The update step in the algorithm to update our guess of the root.
    end
end