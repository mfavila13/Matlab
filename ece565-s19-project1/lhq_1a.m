function G = localhisteq(f,m,n)

A = imread(f);
img = A;

% define window size
M = 7;
N = 7;
M = m;
N = n
win_size = M*N;
hori_pad = round(M/2)-1;
vert_pad = round(N/2)-1;
H = size(A,1);
L = size(A,2);

% add padding to image matrix
padA = padarray(A,[vert_pad,hori_pad],0,'pre');
padA = padarray(padA,[vert_pad,hori_pad],0,'post');
comA = A;
values = zeros(M,N);
occurance = zeros(1,256);
prob = zeros(1,256)
disp('Padded Array Created: Beginning cdf')
for i = 1:H
    for j = 1:L
        % calculate the probability density in window
        % calculate the number of occurances of each pixel value in the
        % window
        for k = 1:M
            for l = 1:N
                values(k,l) = padA(i+hori_pad+k-hori_pad-1,j+vert_pad+l-vert_pad-1);
                val = values(:);
                
            end
        end
        [a,b] = hist(val,unique(val));
        for p = 1:256
            for o = 1:size(b,1)
                if isempty(b) == 0
                    if b(o) == p 
                        occurance(p) = a(o);
                    end
                else
                    occurance(p) = 0;
                end
            end
            
        %disp(m)    
        
        %count the total number of elements in the window
        
        %divide the two
            if isempty(b) == 0
                prob(p) = occurance(p)/win_size;
            else 
                prob(p) = 0;
            end
        end
        %calculate the cumulative probability
        cfact = sum(prob)*255;
        %replace the middle value with the cdf
        comA(i,j) = cfact;
        
        
    end
    disp(i)
end
figure('NumberTitle', 'off', 'Name', '1.c.')
imshow(comA)
G = comA;
end
