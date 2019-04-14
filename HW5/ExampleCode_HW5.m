clc;
clear all;
close all;

%% Q1: Create time-courses to compare the effects of time slice correction on even and odd slices
% Keep in mind that your results are dependent on the choice of reference
% slice. Slice 1 is the preferred reference slice to observe the effect of
% time slice correction

N_timepoints=131;
N_slice=21;
N_x=64;
N_y=64;
rawvolume=zeros(N_x,N_y,N_slice,N_timepoints);
corrvolume=zeros(N_x,N_y,N_slice,N_timepoints);

% Load raw data and corrected data into 4D matrices
for t=1:N_timepoints
    s=sprintf('FUNC4_vol.%03d.img',t);
    fd=fopen(s,'rb');
    
    for k=1:N_slice
        for i=1:N_x
            for  j=1:N_y
                 rawvolume(i,j,k,t)=fread(fd,1,'short');       
            end
        end
    end
    fclose(fd);
end

for t=1:N_timepoints
    s=sprintf('aFUNC4_vol.%03d.img',t);
    fd=fopen(s,'rb');
    
    for k=1:N_slice
        for i=1:N_x
            for  j=1:N_y
                 corrvolume(i,j,k,t)=fread(fd,1,'short');       
            end
        end
    end
    fclose(fd);
end

% Collect data from 4D matrices
even_raw=squeeze(rawvolume(30,36,10,:));
even_corr=squeeze(corrvolume(30,36,10,:));
odd_raw=squeeze(rawvolume(30,36,9,:));
odd_corr=squeeze(corrvolume(30,36,9,:));

% Create plot to see differences
subplot(2,1,1);
plot(1:131,even_raw,'b',1:131,even_corr,'r');
title('Voxel in even slice (36,30,10) before and after time correction');
legend('Original','Corrected');

subplot(2,1,2);
plot(1:131,odd_raw,'b',1:131,odd_corr,'r');
title('Voxel in odd slice (36,30,9) before and after time correction');
legend('Original','Corrected');

%% Q2: Detect head motion

front=squeeze(corrvolume(51,31,13,:)); 
back=squeeze(corrvolume(8,32,11,:));
left=squeeze(corrvolume(28,48,11,:));
right=squeeze(corrvolume(29,14,11,:));
top=squeeze(corrvolume(36,24,19,:));
bottom=squeeze(corrvolume(36,31,5,:));

figure;
subplot(6,1,1);
plot(1:N_timepoints,front);
title('Anterior');
xlabel('Time (s)');
ylabel('BOLD signal');

subplot(6,1,2);
plot(1:N_timepoints,back);
title('Posterior');
xlabel('Time (s)');
ylabel('BOLD signal');

subplot(6,1,3);
plot(1:N_timepoints,left);
title('Left');
xlabel('Time (s)');
ylabel('BOLD signal');

subplot(6,1,4);
plot(1:N_timepoints,right);
title('Right');
xlabel('Time (s)');
ylabel('BOLD signal');

subplot(6,1,5);
plot(1:N_timepoints,top);
title('Superior');
xlabel('Time (s)');
ylabel('BOLD signal');

subplot(6,1,6);
plot(1:N_timepoints,bottom);
title('Inferior');
xlabel('Time (s)');
ylabel('BOLD signal');

% Refer to PDF for explanation

%% Q3: Correcting for bulk head motion with SPM and experimenting with different interpolation algorithms

% Refer to PDF for explanation