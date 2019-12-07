%Marcus Favila
%ECE 437 Computer Project
%07 Dec 2019

%Project Part 1: Questions

global a b M N wmax
[XpF,Fs] = audioread('FilteredSignal_Part1.wav');
for i = 1:length(XpF)
    time(i) = i/Fs;
end

figure('Name','IIR Filter: frequency response','NumberTitle','off');
freqz(b,a)
title('IIR Filter: frequency response')

figure('Name','IIR Filter: frequency response: logarithmic scale','NumberTitle','off');
[h,w] = freqz(b,a,length(XpF),'whole',Fs);
semilogx(w/pi,20*log10(abs(h)));
ax = gca;
ax.YLim = [-180, 20];
ax.XLim = [0,200000];
title('IIR Filter: frequency response')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')

figure('Name','IIR Filter: zero-pole plot','NumberTitle','off');
[vz,vp]=zplane(b,a);
title('IIR Filter: zero-pole plot')
zeros = vz
poles = vp