f = imread('chromosome.tif'); 
figure
imshow(f)
%Extract the outer boundary of gB and display the results as a binary image. 
B=bwboundaries(f); 
d=cellfun('length',B); 
[max_d,k]=max(d); 
b=B{1}; 
[M,N]=size(f); 
g1=bound2im(b,M,N); 
figure 
imshow(g1) 
%Compute the fourier descriptor 
z=frdescp(b); 
%Reconstruct the boundary using 50% of the total possible descriptors 
z1344=ifrdescp(z,1344); 
z1im1344=bound2im(z1344,M,N); 
figure 
imshow(z1im1344) 
%Reconstruct the boundary using 1% of the total possible descriptors 
z26=ifrdescp(z,26); 
z1im26=bound2im(z26,M,N); 
figure 
imshow(z1im26) %%Inverse Fourier descriptor 


