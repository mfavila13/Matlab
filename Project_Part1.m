%Project - part 1 : DO NOT MODIFY
%clear variables
clear

global a b M N 

M = 4; %order of numerator
N = 4; %order of denominator

%filter coefficients 
b= [0.995304813340108, -3.904721241753850, 5.820302744670216, -3.904721241753850, 0.995304813340108];
a =[1.000000000000000, -3.913921547130277, 5.820292921365054, -3.895520936377423, 0.990619449985378];

%read input audio 
[XpF,Fs] = audioread('Signal_plus_feedback.wav');

%plot input audio
figure(1);
plot(XpF);
%play input audio
%sound(XpF,Fs);

%initialize output to all zeros
%y=zeros(1,length(XpF));

for i=1:length(XpF),
    %read input signal sample-by-sample
    x_n = XpF(i);
    % call the iirFilt function to produce the output sample
    y(i)=iirFilt(x_n); 
end

%plot the output audio
figure(2);
plot(y);

%play the output audio
%sound(y,Fs); 

%save the output audio to wave-file
audiowrite('FilteredSignal_Part1.wav',y,Fs);

