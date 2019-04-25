function ECE565_prob1(file) 
image = imread(file);
% mean intensity of image calculated 
[counts,N]=imhist(image);
i=1; 
mu_1=cumsum(counts);
T(i)=(sum(N.*counts))/mu_1(end); 
T(i)=round(T(i));
% mean above and below threshold calculated
mu_2=cumsum(counts(1:T(i)));
MBT=sum(N(1:T(i)).*counts(1:T(i)))/mu_2(end);
mu_3=cumsum(counts(T(i):end));
MAT=sum(N(T(i):end).*counts(T(i):end))/mu_3(end);
i=i+1;
T(i)=round((MAT+MBT)/2);

while abs(T(i)-T(i-1))>=1 
    mu_2=cumsum(counts(1:T(i)));
    MBT=sum(N(1:T(i)).*counts(1:T(i)))/mu_2(end);
    mu_3=cumsum(counts(T(i):end));
    MAT=sum(N(T(i):end).*counts(T(i):end))/mu_3(end);
    i=i+1; T(i)=round((MAT+MBT)/2); Threshold=T(i);
end

% normalize threshold
level = (Threshold - 1) / (N(end) - 1);
BW = im2bw(image,level);
figure(10)
imshow(BW)
title('Segemented Image')
end