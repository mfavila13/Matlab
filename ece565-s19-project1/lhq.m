A=imread('test1.tif'); 
figure('NumberTitle', 'off', 'Name', 'Original Image')
imshow(A);
Img = A;

%Defining the Window size 
m=221;
n=221;
midval=round((m*n)/2);

%Preparation of Matix 

%Step 1: extend matrix for zeros buffer
in=0;
for i=1:m
    for j=1:n
        in=in+1;
        if (in==midval)
            Pad_M = i-1;
            Pad_N = j-1;
            break
        end
    end
end
%Step 2: Pad image with zeros
B = padarray(A,[ Pad_M, Pad_N]);
for i = 1:size(B,1) - ((Pad_M*2) +1)
    for j = 1:size(B,2) - ((Pad_N*2) +1)
        cdf = zeros(256,1);
        inc = 1;
        for x = 1:m
            for y = 1:n
                if (inc == midval)
                    element = B(i+x-1, j+y-1) + 1;
                end
                position = B(i+x-1, j+y-1)+1;
                cdf(position) = cdf(position) + 1;
                inc = inc+1;
            end
        end
        for l = 2:256
            cdf(l) = cdf(l) + cdf(l-1);
        end
            Img(i,j) = round(cdf(element)/m*n)*255;
    end
end

figure('NumberTitle', 'off', 'Name', 'Manual Local Histogram Equalization')
subplot(2,2,3)
imhist(Img);
subplot(2,2,4)
imhist(A)
subplot(2,2,1)
imshow(Img);



