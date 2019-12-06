
%PART2
% clear all variables
clear
%read audio file 
[XpF,Fs] = audioread('Signal_plus_feedback.wav');

%plot audio file
figure(1);
plot(XpF);

%play audio 
%sound(XpF,Fs);
L=256; 
%design FIR filter of order L (i.e. length L+1)
b_fir=fir1(L,[0.04 0.08],'stop'); 

%use block size N that is the next power of 2 compared to size of b_fir
% e.g. 256, 512, 1024, 2048 ... 
t = nextpow2(length(b_fir));
N = 2^t; 

% make FIR length equal to block-size, i.e. 
h=zeros(N,1);
%copy coefficients 
h(1:length(b_fir))=b_fir; 
H=fft(h); % H is N-point fft

num_Blocks = length(XpF)/N;
%initialize output vector
y=[]; 
for block = 1:num_Blocks,
    %get block of input audio data
    % (NOTE: in a real-time implementation, you would get only 1 
    % block of data at a time to work with) 
    x = XpF((block-1)*N+1:block*N);
    
    X=fft(x);
    % multiply N-point ffts
    Y=X.*H;
    
    %get the N-point ifft of the output
    y_fft=ifft(Y);
    
    % append block to the output 
    y=[y;y_fft];
    block
end

%plot final output 
figure(2);
plot(y);

%play final output
%sound(y,Fs); 

%write the output to a wave file
audiowrite('FilteredSignal_Part2.wav',y,Fs);