f = @(x) (tanh(x));  % The function f is defined here.
a = 13;  % First term of the sequence u such that f(a)>0
b = -7;  % First term of the sequence v such that f(b)<0 
[u,v] = solve(f,a,b);  % The function returns two vectors each storing one sequence.
%disp("The function f is:" + f);
disp("Sequence u = ");   % Finally we print both the sequences.
disp(u);
disp("The value of the function f for the sequence u are:"); 
disp(f(u));
disp("Sequence v = ");
disp(v);
disp("The value of the function f for the sequence v are:");
disp(f(v));
function [u,v] = solve(f,a,b)
    u = [];  % Empty vectors are initialized to store the sequences.
    v = [];
    count = 1; 
    while(count<=10)  % We limit the sequences to the first 10 terms as mentionede in the question.
        u = [u a];  % We keep adding the next terms to the corresponding sequences.
        v = [v b];
        num = a*f(b) - b*f(a);  
        den = f(b) - f(a);
        w = num/den;  % The value of w is calculated as mentioned in the question.
        if(f(w)*f(a) > 0)  % We check for the update condition of a and b as mentioned in the question.
            a = w;
        else
            b = w;
        end
        count = count + 1; % The count variable is incremented until it crosses the upper limit.
    end
end