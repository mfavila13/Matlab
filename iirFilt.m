function y_n = iirFilt(x_n)

%Marcus Favila
%ECE 437 Computer Project
%07 Dec 2019

% 'clear iifFilt' should be entered in Command Window before running 
% Project_Part1.m 
% if not cleared, peristent variables will cause calculations to be
% inaccurate resulting in an unrepresentative output signal for the 
% filter

global a b M N wmax
%These values are used to store the previous values for the calculation
%of y(n)
persistent w;
persistent X;
persistent Y;
persistent count;
%addend refers to the term inside the summation at some k or j or l
B_addend = 0;
A_addend = 0;
W_addend = 0;
%term refers to the sum of all terms in the summation for range of k or j
%or l
B_term = 0;
A_term = 0;
W_term = 0;

%count is analogous to n in y(n) or x(n)
%count is set to 1 on first call of function: persistent variables should
%be cleared with every run of a script that calls the iirFilt function
if isempty(count)
    count = 1;
end

if isempty(w)
    w=zeros(1,N+1);
    wmax = 0;
end
if isempty(X)
    X=zeros(1,N+1);
end
if isempty(Y)
    Y=zeros(1,N+1);
end
X(count) = x_n;

%Calculation of the x(n) summation
for k = 1:M+1
    if (count-k) < 1
        B_addend = 0;
    else
        B_addend = b(k)*X(count-k+1);
    end
    B_term = B_term + B_addend;
end
%disp(B_term)
%Calculation of the y(n) summation
for j = 2:N+1
    if (count-j) < 1
        A_addend = 0;
    else
        A_addend = a(j)*Y(count-j+1);
    end
    A_term = A_term + A_addend;
end
%disp(A_term)

for l = 2:N+1
    if(count-l+1) < 1
        W_addend = 0;
    else
        W_addend = (a(l)*w(count-l+1));
    end
    W_term = W_term + W_addend;
end
w(count) = -1*W_term + X(count);
if w(count) > wmax
    wmax = w(count);
end
%calulates the final y(n) value 
y_n = B_term - A_term;

Y(count) = y_n;
disp(count)
count = count + 1;
%pause(2);
end

    

