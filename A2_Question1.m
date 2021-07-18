A = [7 2 4;1 9 5;6 3 8];  %Input Matrix A
b = [2;3;5]; %Input Vector b such that the solution to AX=b is to be found
[B,L,U,b,flag] = decompose(A,b);  % Function call to decompose A into lower and upper matrix
A = B; % The matrix B is a modified form of A if A had some diagonal entry 0 and if modification was needed flag is 1 else 0.
%In the function "decompose" we add an additional functionality where in case
% some A(i,i) = 0,we swap i_th row of A and b with some other row where the
% pivot element is not 0. Hence, the function also returns an updated
% matrix A and updated vector b.

% Displaying the results from the decompose function
if(flag == 1)  % Flag = 1 implies that modification was needed
    disp("Modified Matrix A = ");
    disp(A);
else
    disp("Matrix A = ");
    disp(A);
end
disp("Matrix L = ");
disp(L);
disp("Matrix U = ");
disp(U);
disp("The product of L and U is = ");
disp(L*U);  % Verifying if the matrix A has been decomposed correctly

X = solve(L,U,b); %Function call to solve the system of linear equations AX = b
%Displaying the result
if(flag == 1) % Flag = 1 implies b vector was also modified with A
    disp("Modified vector b = ");
    disp(b);
else
    disp("Vector b = ");
    disp(b);
end
disp("Solution X = ");
disp(X);
disp("The product of A and X ="); %Checking if the product A*X is equal to b
disp(A*X);

if(flag == 1) %Explaining the role of the flag variable
    disp("Matrix A and vector B were modified because some diagonal element in matrix A was zero and rows in A and b were needed to be swapped");
end

%Function to decompose A into L and U such that A = L*U where L is a lower
%triangular matrix and U is an upper triangular matrix/
function [B,L,U,b,flag] = decompose(A,b)
    flag = 0;
    n = size(A,1);
    M = zeros(n,n); %Matrix M to store the values m_ij
    B = A;
    for k = 1:n
        for i = k:n %We check the rows below to find a row such that the pivot element is non zero
            if(A(i,k) ~= 0)
                if(i ~= k) %If i = k, it implies that the current row is good to work with and hence no swaps needed
                    flag = 1;  % If even 1 swap was made, it means that A has been modified and hence flag = 1
                    B([i k],:) = B([k i],:); % We maintain the swaps made in a copy of A called B
                    A([i k],:) = A([k i],:);  %We swap the rows accordingly
                    b([i k],:) = b([k i],:);  %Vector b is swapped accordingly to maintain the original state of the system of linear equations.
                end
                break;
            end
        end   
        for i = 1:n
            M(i,k) = A(i,k)/A(k,k); %Values for the k_th column of matrix M
        end
        for i = (k+1):n
            for j = 1:n
                A(i,j) = A(i,j) - M(i,k)*A(k,j); % Calculating (A_ij)^{k+1} for all entries of A from row k+1 to n
            end
        end
    end
    L = zeros(n,n);
    U = zeros(n,n);
    %Using the result from the lecture slides we fill in the values for the
    %upper and lower matrices L and U using the matrices M and A
    for i = 1:n
        for j = 1:i
            L(i,j) = M(i,j);
        end
        for j = i:n
            U(i,j) = A(i,j);
        end
    end
end

%Function to calculate the solution to the system of equations AX = b which is
%equivalent to LUX = b
function X = solve(L,U,b)
    n = size(L,1);
    %We need to find solutions to the LUX = b 
    %Let Z = UX
    %Therefore, we first find solution to LZ = b which can be found easily
    %using the forward substitution method.
    Z = zeros(n,1);
    for i = 1:n
        sum = b(i,1);
        j = i-1;
        while(j>0)
            sum = sum - (L(i,j)*Z(j,1));
            j = j-1;
        end
        Z(i,1) = sum/L(i,i);
    end
    
    %Since we have a Z such that UX = Z, we can find X easily using the
    %backward substitution method.
    X = zeros(n,1);
    for i = n:-1:1
        sum = Z(i,1);
        j = i+1;
        while(j <= n)
            sum = sum - (U(i,j)*X(j,1));
            j = j+1;
        end
        X(i,1) = sum/U(i,i);
    end
end