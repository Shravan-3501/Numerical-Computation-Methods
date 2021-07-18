%A and b matrices defined for the system of equations Ax = b
A = [10 1;1 10];
b = [11;11];
guess = [1/2;1/2]; %Initial guess for the solution x
disp("(i)Gauss-Seidel Method");
disp("A = ");
disp(A);
disp("b = ");
disp(b);
x = Gauss_Seidel(A,b,guess); %Solving the equations using gauss-seidel method
disp("x = ");
disp(x); 
disp("A*x = "); %Verifying if the solution is correct
disp(A*x); 

%A and b matrices defined for the system of equations Ax = b
A = [4 1 -1;2 7 1;1 -3 12];
b = [3;19;31];
guess = [0;0;0]; %Initial guess for the solution x
disp("(ii)Gauss-Jacobi Method");
disp("A = ");
disp(A);
disp("b = ");
disp(b);
x = Gauss_Jacobi(A,b,guess); %Solving the equations using Gauss-Jacobi method
disp("x = ");
disp(x); 
disp("A*x = "); %Verifying if the solution is correct
disp(A*x); 

function x = Gauss_Seidel(A,b,guess)
    epsilon = 10^-4; %Error tolerance
    n = size(A,1);
    x = guess;
    for i = 1:n
        x(i) = b(i);
        for j = 1:n
            if(j ~= i)
                x(i) = x(i) - A(i,j)*x(j); %Gauss-Seidel Iterative formula
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
                    x(i) = x(i) - A(i,j)*x(j); %Gauss-Seidel Iterative formula
                end
            end
            x(i) = x(i)/A(i,i);
        end
    end
end

function x = Gauss_Jacobi(A,b,guess)
    epsilon = 10^-4; %Error tolerance
    n = size(A,1);
    x = zeros(n,1);
    for i = 1:n
        x(i) = b(i);
        for j = 1:n
            if(j ~= i)
                x(i) = x(i) - A(i,j)*guess(j); %Gauss-Jacobi Iterative Formula
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
                    x(i) = x(i) - A(i,j)*guess(j); %Gauss-Jacobi Iterative Formula
                end
            end
            x(i) = x(i)/A(i,i);
        end
    end
end

%Function to calculate error between two consecutinve guesses using L-infinity norm
function e = error(guess,x)
    e = 0;
    n = size(x,1);
    for i = 1:n
        e = max(e,abs(guess(i) - x(i))); %L-infinty norm over a vector
    end
end