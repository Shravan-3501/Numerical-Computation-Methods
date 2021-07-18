%Input vectors a,b,c,d
a = [0;7;2;7;1]; 
b = [1;9;2;4;8];
c = [8;1;6;4;0];
d = [9;7;1;9;2];
%Here, it is more logical to take input as these column vector a,b,c,d
%instead of the usual n*n matrix A because since these system of linear
%equations are triadiagonal majority of the entries in the matrix A will be
%0 hence wasting a lot of memory.
disp("Vector a = ");
disp(a);
disp("Vector b = ");
disp(b);
disp("Vector c = ");
disp(c);
disp("Vector d = ");
disp(d);
[X,flag] = tridiagonal(a,b,c,d); %Output from the function after solving the tridiagonal system of linear equations

%The output flag signifies if the specified system can be solved by this
%method or not. If flag = 1, it means we have encountered a case where the
%denominator becomes 0 and hence this method cannot be used.
if(flag == 0) %Therefore, we only display the output X if flag = 0 
    disp("The solution for the system of equations X = ");
    disp(X);
else %If flag = 1 we display an error message instead
    disp("Denominator becomes 0, hence unable to find solution");
end

%Function to solve a tridiagonal system of linear equations.
function [X,flag] = tridiagonal(a,b,c,d)
    flag = 0;
    n = size(a,1);
    X = zeros(n,1);
    c_dash = zeros(n,1);
    d_dash = zeros(n,1);
    
    %Using the formula given in the question we calculate the vectors
    %c_dash and d_dash
    for i = 1:n
        if(i == 1)
            c_dash(i,1) = c(i,1)/b(i,1);
            d_dash(i,1) = d(i,1)/b(i,1);
        else
            den = b(i,1) - a(i,1)*c_dash(i-1,1);
            if(den == 0) %If the denominator becomes 0 at any stage we return from the function call with the value of flag as 1 signifying an error
                flag = 1;
                return;
            end
            c_dash(i,1) = c(i,1)/den;
            num_d = d(i,1) - a(i,1)*d_dash(i-1,1);
            d_dash(i,1) = num_d/den;
        end
    end
    %Using the formula given in the question we obtain the solution vector X
    for i = n:-1:1
        if(i == n)
            X(i,1) = d_dash(i,1);
        else
            X(i,1) = d_dash(i,1) - (c_dash(i,1)*X(i+1,1));
        end
    end
end