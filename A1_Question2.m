f = @(x) (x^4 - 16*(x^3) + 89*(x^2) - 194*x + 120);  %Defined the function whose root is to be calculated.
D = @(x) ((f(x + f(x)) - f(x))/f(x)); %The function D to be used for updating guess defined in terms of the function f.
guess = 1.5; % Initial guess as mentioned in the question.
[Root,p] = solve(f,D,guess);  % The function solve returns the root of the equation along with vector of values of p(The order of convergence at each step).
disp("Root of the equation = " + Root); %Printing the root of the equation.
disp("Series of values for p:");
disp(p); %Printing the vector of order of convergence.
disp("Therefore, the value of Order of Convergence is :");
disp(p(size(p,2)));  %Printing the value of the Order of Convergence where the vector p converges.

function [guess_3,p] = solve(f,D,g)
    epsilon = 10^-10;  % Value of epsilon for the stopping criterion.
    guess_0 = g;  % Initial guess.
    p =[];
    %Since to calculate a value of p we need 3 consecutive error terms.
    % To calculate the 3 error terms we need 4 consecutive guesses.
    %Therefore, we run the algorithm 3 times first to obtain first 4 values
    %of guesses and hence the initial approximation to p(The order of convergence); 
    difference = f(guess_0)/D(guess_0);
    guess_1 = guess_0 - difference;
    difference = f(guess_1)/D(guess_1);
    guess_2 = guess_1 - difference;
    difference = f(guess_2)/D(guess_2);
    guess_3 = guess_2 - difference;
    while(abs(f(guess_3)) > epsilon)  %Stopping criterion in terms of epsilon
        en3 = abs(guess_3 - guess_2); % Error term calculated using the formula mentioned in the question.
        en2 = abs(guess_2 - guess_1);
        en1 = abs(guess_1 - guess_0);
        p1 = log(en3/en2);
        p2 = log(en2/en1);
        temp = p1/p2;  % Order of convergence calculated using the formula mentioned in the question.
        p = [p temp]; % The p value added to the vector.
        difference = f(guess_3)/D(guess_3);
        guess_0 = guess_1;                
        guess_1 = guess_2;
        guess_2 = guess_3;
        guess_3 = guess_3 - difference;  % We finally update the guess using Newton-Raphson Method formula.
    end
end