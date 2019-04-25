function z = fourierdescp(s) 
[np, nc] = size(s); 
if nc ~= 2 
    error('S must be of size np-by-2.'); 
end 
if np/2 ~= round(np/2); 
    s(end + 1, :) = s(end, :); 
    np = np + 1; 
end
x = 0:(np - 1); 
m = ((-1) .^ x)'; 
% center transform
s(:, 1) = m .* s(:, 1); s(:, 2) = m .* s(:, 2); 
% convert to complex value
s = s(:, 1) + i*s(:, 2); 
% compute descriptors
z = fft(s);
end