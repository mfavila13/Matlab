%Marcus Favila
%ECE 437 Computer Project
%07 Dec 2019
%PART2

% clear all variables
clear
clear iirFilt;
time = [];
xmin  = 32500;
xmax = 42500;

%read audio file 
[XpF,Fs] = audioread('Signal_plus_feedback.wav');
for i = 1:length(XpF)
    time(i) = i/Fs;
end

%plot audio file
figure('Name','Signal with Feedback','NumberTitle','off');
plot(time(1,xmin:xmax),XpF(xmin:xmax,1));
title('Signal With Feedback')
xlabel('Time (secs) ')
ylabel('Signal')

%play audio 
%sound(XpF,Fs);
%_______________________________________________________________________

%First FIR Filter: DFT/IDFT values provided

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
figure('Name','Signal: FIR Filter Applied: {original block size}','NumberTitle','off');
plot(time(1,xmin:xmax),y(xmin:xmax,1));
title('Signal: FIR Filter Applied: {original block size}')
xlabel('Time (secs) ')
ylabel('Signal')
%play final output
%sound(y,Fs); 

%write the output to a wave file
audiowrite('FilteredSignal_Part2a_DFT_IDFT_originalblock.wav',y,Fs);

%_________________________________________________________________________

%Second Filter: time-domain filtered output
y2=[];
b2 = b_fir;
a2=[1];
y2 = filter(b2,a2,XpF);

%plot the output audio
figure('Name','Signal: FIR Filter Applied: {matlab 1-D digital filter}','NumberTitle','off');
plot(time(1,xmin:xmax),y2(xmin:xmax));
title('Signal: FIR Filter Applied: {matlab 1-D digital filter}')
xlabel('Time (secs) ')
ylabel('Signal')
%play the output audio
%sound(y2,Fs); 

%save the output audio to wave-file
audiowrite('FilteredSignal_Part2b_FIR_filter_timedomain.wav',y2,Fs);


%__________________________________________________________________________

%Third Filter: Modified N: DFT/IDFT based filtering
%Chosen value of N: 262144 or 2^18

%use block size N that is the next power of 2 compared to size of b_fir
% e.g. 256, 512, 1024, 2048 ... 
t = nextpow2(length(b_fir));
N3 = 2^(2*t); 

% make FIR length equal to block-size, i.e. 
h=zeros(N3,1);
%copy coefficients 
h(1:length(b_fir))=b_fir; 
H=fft(h); % H is N-point fft

num_Blocks3 = length(XpF)/N3;
%initialize output vector
y3=[]; 
for block3 = 1:num_Blocks3,
    %get block of input audio data
    % (NOTE: in a real-time implementation, you would get only 1 
    % block of data at a time to work with) 
    x = XpF((block3-1)*N3+1:block3*N3);
    
    X=fft(x);
    % multiply N-point ffts
    Y=X.*H;
    
    %get the N-point ifft of the output
    y_fft=ifft(Y);
    
    % append block to the output 
    y3=[y3;y_fft];
    block3
end

%plot final output 
figure('Name','Signal: FIR Filter Applied: {modified block size; (N = 2 2t) }','NumberTitle','off');
plot(time(1,xmin:xmax),y(xmin:xmax,1));
title('Signal: FIR Filter Applied: {modified block size; (N = 2 2t) }')
xlabel('Time (secs) ')
ylabel('Signal')
%play final output
%sound(y,Fs); 

%write the output to a wave file
audiowrite('FilteredSignal_Part2c_DFT_IDFT_modifiedblock.wav',y,Fs);
%_______________________________________________________________________
%{

% Altering the size of N such that the number of blocks used results in
% memory errors as array becomes too large for program to compute
% This section is omitted to reflect that

%use block size N that is the next power of 2 compared to size of b_fir
% e.g. 256, 512, 1024, 2048 ... 
t = nextpow2(length(b_fir));
N3 = 2^(3*t); 

% make FIR length equal to block-size, i.e. 
h=zeros(N3,1);
%copy coefficients 
h(1:length(b_fir))=b_fir; 
H=fft(h); % H is N-point fft

num_Blocks3 = length(XpF)/N3;
%initialize output vector
y3=[]; 
for block3 = 1:num_Blocks3,
    %get block of input audio data
    % (NOTE: in a real-time implementation, you would get only 1 
    % block of data at a time to work with) 
    x3 = XpF((block3-1)*N3+1:block3*N3);
    
    X=fft(x3);
    % multiply N-point ffts
    Y=X.*H;
    
    %get the N-point ifft of the output
    y_fft=ifft(Y);
    
    % append block to the output 
    y3=[y;y_fft];
    block3
end

%plot final output 
figure(5);
plot(y3(xmin:xmax,1));
title('modified block size overlap-add (N = 2 3t) Method')
%play final output
%sound(y,Fs); 

%write the output to a wave file
audiowrite('FilteredSignal_Part2e.wav',y,Fs);
%}
%___________________________________________________________________

%Same function as Part1
%IIR Filter for comparison with FIR filters
for i=1:length(XpF),
    %read input signal sample-by-sample
    x_n = XpF(i);
    % call the iirFilt function to produce the output sample
    y_iir(i)=iirFilt(x_n);
end

%plot the output audio
figure('Name','Signal: IIR Filter Applied: LTI filter using DFII','NumberTitle','off');
plot(time(1,xmin:xmax),y_iir(1,xmin:xmax));
title('Signal: IIR Filter Applied: LTI filter using DFII');
xlabel('Time (secs) ')
ylabel('Signal')
audiowrite('FilteredSignal_Part2d_IIR_Filter.wav',y,Fs);
