function ECE565_prob4(file)
image = imread(file); 
figure(40)
imshow(image)
%Extract the outer boundary of gB and display the results as a binary image. 
B=bwboundaries(image); 
d=cellfun('length',B); 
[max_d,k]=max(d); 
b=B{1}; 
[M,N]=size(image); 
g1=bound2im(b,M,N); 
figure(41) 
imshow(g1) 
%Compute the fourier descriptor 
z=fourierdescp(b); 
%Reconstruct the boundary using 50% of the total possible descriptors 
z1344=ifourierdescp(z,1344); 
z1im1344=bound2im(z1344,M,N); 
figure(42) 
imshow(z1im1344) 
%Reconstruct the boundary using 1% of the total possible descriptors 
z26=ifourierdescp(z,26); 
z1im26=bound2im(z26,M,N); 
figure(43) 
imshow(z1im26) %%Inverse Fourier descriptor 
end

