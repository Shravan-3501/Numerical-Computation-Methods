%OBSERVATION MADE:Since the matrix defined in the question is symmetric, the row and column
%norm will always be equal and hence the conditional numbers obtained using
%the row and column norms will also be equal as is observed in the output.

%Iterating over the values of n i.e 3,4,5,6
for n = 3:6
    [row_cond,col_cond,eucl_cond] = conditional(n);
    disp("Conditional number using row norm for n = " + n + " is : " + row_cond);
    disp("Conditional number using column norm for n = " + n + " is : " + col_cond);
    disp("Conditional number using eucledian norm for n = " + n + " is : " + eucl_cond);
end
disp("Since the matrix defined in the question is symmetric, the row and column norm are always equal"); 
disp("and hence the conditional numbers obtained using the row and column norms are also equal as is observed in the output above.");

%Function to calculate conditional number of a matrix using all the 3 norms
%defined i.e Row norm, Column norm, Eucledian norm
function [row_cond,col_cond,eucl_cond] = conditional(n)
    %For a given n, row_cond is the conditional number obtained by using
    %row norm. col_cond is the conditional number obtained by using
    %Column norm and eucl_cond is the conditional number obtained by using
    %Eucledian norm. 
    H_n = zeros(n,n);  %H_n is the matrix prepared as asked in the question for a given n 
    for i = 1:n
        for j = 1:n
            H_n(i,j) = 1/(i+j-1);
        end
    end
    
    I_n = inv(H_n); % I_n is the inverse of the matrix H_n
    
    r1 = calc_row_norm(H_n);  % Calculate the row norms for H_n and I_n using a separate function defined below.
    r2 = calc_row_norm(I_n);
    row_cond = r1*r2;  %Calculate the row conditional number using the formula row_norm(H_n)*row_norm(I_n) 
    
    c1 = calc_col_norm(H_n); % Calculate the column norms for H_n and I_n using a separate function defined below.
    c2 = calc_col_norm(I_n);
    col_cond = c1*c2; %Calculate the column conditional number using the formula col_norm(H_n)*col_norm(I_n)
    
    e1 = calc_eucl_norm(H_n); % Calculate the eucledian norms for H_n and I_n using a separate function defined below.
    e2 = calc_eucl_norm(I_n);
    eucl_cond = e1*e2; %Calculate the eucledian conditional number using the formula eucl_norm(H_n)*eucl_norm(I_n)
end

%Function to calculate the column norm of any matrix A
function col_norm = calc_col_norm(A)
    A = abs(A);
    x = sum(A,1); % We sum over all the columns of A(i.e axis = 1)
    col_norm = max(x); %The maximum sum is the col_norm according to the definition of column norm
end

%Function to calculate the row norm of any matrix A
function row_norm = calc_row_norm(A)
    A = abs(A);
    x = sum(A,2); % We sum over all the rows of A(i.e axis = 2)
    row_norm = max(x); %The maximum sum is the row_norm according to the definition of row norm
end

%Function to calculate the eucledian norm of any matrix A
function eucl_norm = calc_eucl_norm(A)
    s = 0;
    n = size(A,1);
    for i = 1:n
        for j = 1:n
            s = s + (A(i,j)*A(i,j)); %Simply add squares of all elements of the matrix A
        end
    end
    eucl_norm = s^0.5; %The eucledian norm is the root of sum of squares of all elements of A 
end