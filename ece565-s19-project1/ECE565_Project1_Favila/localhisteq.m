% Marcus Favila
% localhisteq() function that applies local histogram eqaulization of
% window size m x n on image f
% returns image G
% default window size is 3x3 if omitted.

function G = localhisteq(f,m,n)

A = imread(f);
img = A;

% if arguments are left void, assign default window size of 3x3
if nargin <2
    m = 3;
    n = 3;
    disp('*** Window size set to default (3x3) ***')
end
if rem(m,2) ==0
    disp('*** ERROR: Dimension ''m'' must be odd... ***')
    m = input('Please enter an odd value for m: ');
end
if rem(n,2) ==0
    disp('*** ERROR: Dimension ''n'' must be odd... ***')
    n =  input('Please enter an odd value for n: ');
end
M = m;
N = n;
%define winow size
win_size = M*N;
hori_pad = round(N/2)-1;
vert_pad = round(M/2)-1;
H = size(A,1);
L = size(A,2);
% add padding to image matrix
padA = padarray(A,[vert_pad,hori_pad],0,'pre');
padA = padarray(padA,[vert_pad,hori_pad],0,'post');
%variables for computations
comA = A; %ouput image reference
values = zeros(M,N); %lists values of image within window
occurance = zeros(1,256); % lists occurance for each intenisty value within the window
prob = zeros(1,256); %probability densities
disp('Padded Array Created: Beginning local histogram equalization')
for i = 1:H
    for j = 1:L
        %obtain values within window
        for k = 1:N
            for l = 1:M
                values(k,l) = padA(i+(l-1),j+(k-1));
                val = values(:);      
            end
        end
        %compute the number of ocurances within the window, list values in
        %ascending order
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
            % calculate the probability density of each value
            if isempty(b) == 0
                prob(p) = occurance(p)/win_size;
            else 
                prob(p) = 0;
            end
        end
        %calculate the cumulative distribution value
        cfact = sum(prob)*255;
        %replace the middle value with the cdv
        comA(i,j) = cfact;      
    end
    if i == H/4
        disp('25% complete...')
    end
    if i == H/2
        disp('50% complete...')
    end
    if i == 3*H/4
        disp('75% complete...')
    end
    if i == H
        disp('100% complete...outputing image')
    %disp(i) %used to keep track of processing progress 
    end
end
%final output
G = comA;
figure('NumberTitle', 'off', 'Name', 'Local Histogram Equalization')
imshow(G);
