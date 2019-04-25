function s = ifourierdescp(z, nd) 

np = length(z); 
if nargin == 1 | nd > np 
    nd = np; 
end
x = 0:(np - 1); 
m = ((-1) .^ x)'; 
% number of descriptors defined to be used
d = round((np - nd)/2);  
z(1:d) = 0; 
z(np - d + 1:np) = 0; 
% compute inverse and convert to coordinates. 
zz = ifft(z); 
s(:, 1) = real(zz); 
s(:, 2) = imag(zz); 
% reverse centering
s(:, 1) = m.*s(:, 1); 
s(:, 2) = m.*s(:, 2);
end
