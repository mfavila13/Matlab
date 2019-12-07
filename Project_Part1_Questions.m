function part1q = Project_Part1_Questions()
% Marcus Favila
%Project Part 1: Questions

%clear variables
clear
clear iifFilt
global a b M N wmax
[XpF,Fs] = audioread('FilteredSignal_Part1.wav');
for i = 1:length(XpF)
    time(i) = i/Fs;
end

[h,w] = freqz(b,a,length(XpF),'whole',Fs);


figure(6)
semilogx(w/pi,20*log10(abs(h)));
ax = gca;
ax.YLim = [-180, 20];
ax.XLim = [0,20000];

figure(7)
[vz,vp]=zplane(b,a);
end