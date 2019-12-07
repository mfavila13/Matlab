%Marcus Favila
%ECE 437 Computer Project
%07 Dec 2019

%Project Part 2: Questions
L = 256;
b_fir=fir1(L,[0.04 0.08],'stop');
a_fir=[1];

[XpF,Fs] = audioread('FilteredSignal_Part1.wav');
for i = 1:length(XpF)
    time(i) = i/Fs;
end

figure('Name','FIR Filter: frequency response','NumberTitle','off');
freqz(b_fir);
title('FIR Filter: frequency response')

figure('Name','FIR Filter: frequency response: logarithmic scale','NumberTitle','off');
[h,w] = freqz(b_fir,a_fir,length(XpF),'whole',262144);
semilogx(w/pi,20*log10(abs(h)));
ax = gca;
ax.YLim = [-180, 20];
ax.XLim = [0,20000];
title('FIR Filter: frequency response')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')

figure('Name','FIR Filter: zero-pole plot','NumberTitle','off');
[vz,vp]=zplane(b_fir,a_fir);
title('FIR Filter: zero-pole plot')
zeros = vz
poles = vp
