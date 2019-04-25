function ECE565_prob3(file)
% show original image
image = imread(file); 
figure(30)
imshow(image);
title('Original Image')

w=fspecial('average',9); % built in matlab function: 2D filter averaging filter with size 9 x 9
g = imfilter(image,w,'symmetric'); % built in matlab function: filtering by smoothing image using symmetric option
figure(31) 
imshow(g,[]) 
level = graythresh(image); 
title('Smoothed image (9x9 averaging filter)')

% binary image 
BW = im2bw(g,level); 
figure(32) 
imshow(BW) 
title('Binary image')

% outer boundary
B=bwboundaries(BW); 
d=cellfun('length',B); 
[max_d,k]=max(d); 
b=B{1}; 
[M,N]=size(BW); 
g1=bound2im(b,M,N); 
figure(33)
imshow(g1) 
title('Outer Boundary')

% subsample boundary, display points 
[s,su]=bsubsamp(b,50); 
g2=bound2im(s,M,N); 
figure(34) 
imshow(g2) 
title('Outer Boundary, subsampled (points only)')

% subsample boundary, connect points
cn=connectpoly(s(:,1),s(:,2)); 
g3=bound2im(cn,M,N); 
figure(35) 
imshow(g3) 
title('Outer Boundary, subsampled (points connected)')
fchcode(su,4)
fchcode(su,8)
end
