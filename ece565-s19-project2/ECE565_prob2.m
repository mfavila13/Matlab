function ECE565_prob2(file)
% implementation of Otsu's algorithm
image = imread(file);
mgk = 0;
mt = 0;
[m,n] = size(image);
h = imhist(image);
prob = h/(m.*n);

for i=1:1:256 
    if prob(i)~=0
        low_var=i;
        break 
    end
end

for i=256:-1:1 
    if prob(i)~=0 
        high_var=i;
        break 
    end
end

for k = 1:256
    prob_1(k)=sum(prob(1:k)); 
    prob_2(k)=sum(prob(k+1:256));
    mgk=(k-1)*prob(k)+mgk;
end

for k=1:256
    mean_1(k)=sum((k-1)*prob(1:k))/prob_1(k); 
    mean_2(k)=sum((k-1)*prob(k+1:256))/prob_2(k); 
end
for k =1:256 
    var(k)=prob_1(k)*(mean_1(k)-mgk)^2+prob_2(k)*(mean_2(k)-mgk)^2; 
end
[y,T]=max(var(:));
T=T+low_var;
g=image; 
g1=find(g>=T); 
g(g1)=255; 
g2=find(g<T); 
g(g2)=0;
figure(20)
imshow(g)
title('Otsu thresholding algorithm')
end