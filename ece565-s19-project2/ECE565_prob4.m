function ECE565_prob4(file)
% show original image
image = imread(file); 
figure(40)
imshow(image)
title('Original Image')

% outer boundary
B=bwboundaries(image); 
d=cellfun('length',B); 
[max_d,k]=max(d); 
b=B{1}; 
[M,N]=size(image); 
g1=bound2im(b,M,N); 
figure(41) 
imshow(g1) 
title('Extracted Boundary of Chromosome')

% fourier descriptor 
z=fourierdescp(b) 

% reconstruction using 50% of the total possible descriptors 
z1344=ifourierdescp(z,1344) 
z1im1344=bound2im(z1344,M,N); 
figure(42) 
imshow(z1im1344) 
title('Reconstruction using 50% of Total Possible Descriptors')

% reconstruction using 1% of the total possible descriptors 
z26=ifourierdescp(z,26) 
z1im26=bound2im(z26,M,N); 
figure(43) 
imshow(z1im26) %%Inverse Fourier descriptor 
title('Reconstruction using 1% of Total Possible Descriptors')
end

