function ECE565_prob2(file)
image = imread(file);
mgk = 0;
mt = 0;
[m,n] = size(image);
h = imhist(image); pi = h/(m.*n);
for i=1:1:256 
    if pi(i)~=0
        lv=i;
        break 
    end
end
for i=256:-1:1 
    if pi(i)~=0 hv=i;
        break 
    end
end
lh = hv - lv;
    for k = 1:256
        p1(k)=sum(pi(1:k)); 
        p2(k)=sum(pi(k+1:256));
    end
for k=1:256
    m1(k)=sum((k-1)*pi(1:k))/p1(k); 
    m2(k)=sum((k-1)*pi(k+1:256))/p2(k); 
end
for k=1:256 
    mgk=(k-1)*pi(k)+mgk;
end
for k =1:256 
    var(k)=p1(k)*(m1(k)-mgk)^2+p2(k)*(m2(k)-mgk)^2; 
end
[y,T]=max(var(:));
T=T+lv;
g=image; 
g1=find(g>=T); 
g(g1)=255; 
g2=find(g<T); 
g(g2)=0;
figure(20)
imshow(g)
end