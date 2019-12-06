function y_n = iirFilt(x_n)
global a b M N wmax
persistent w;

if isempty(w)
    w=zeros(1,N+1);
end

% CONTINUE IMPLEMENTING THE IIR DF-2 FILTER ... 
    

