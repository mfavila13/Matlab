function y_n = iirFilt(x_n)
global a b M N wmax
persistent w;
persistent X;
persistent Y;
persistent count;
%persistent B_addend;
%persistent A_addend;
B_addend = 0;
A_addend = 0;
B_term = 0;
A_term = 0;
if isempty(count)
    count = 1;
end
if isempty(w)
    w=zeros(1,N+1);
end
if isempty(X)
    X=zeros(1,N+1);
end
if isempty(Y)
    Y=zeros(1,N+1);
end
%if isempty(A_addend)
%    A_addend=zeros(1,N+1);
%end
%if isempty(B_addend)
%    B_addend=zeros(1,N+1);
%end

X(count) = x_n;

%Y(count) = y_n;

for k = 1:M+1
    if (count-k) < 1
        B_addend = 0;
    else
        B_addend = b(k)*X(count-k+1);
    end
    B_term = B_term + B_addend;
end
disp(B_term)
for j = 2:N+1
    if (count-j) < 1
        A_addend = 0;
    else
        A_addend = a(j)*Y(count-j+1);
    end
    A_term = A_term + A_addend;
end
disp(A_term)

y_n = B_term - A_term;
Y(count) = y_n;
disp(count)
count = count + 1;
%pause(2);
end

% CONTINUE IMPLEMENTING THE IIR DF-2 FILTER ... 
    

