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
%selection for plot, in order: x, y, even slice, odd slice
sel = [40,15,14,15];
% Load raw data and corrected data into 4D matrices
for t=1:N_timepoints
    s=sprintf('FUNC4_vol.%03d.img',t);
    fd=fopen(s,'rb');
    a_s=sprintf('aFUNC4_vol.%03d.img',t);
    a_fd=fopen(a_s,'rb');
    for k=1:N_slice
        for i=1:N_x
            for  j=1:N_y
                 rawvolume(i,j,k,t)=fread(fd,1,'short');
                 corrvolume(i,j,k,t)=fread(a_fd,1,'short'); 
            end
        end
    end
    waitbar(t/N_timepoints)
    fclose(fd);
    fclose(a_fd);
end



% Collect data from 4D matrices
even_raw=squeeze(rawvolume(sel(1),sel(2),sel(3),:));
even_corr=squeeze(corrvolume(sel(1),sel(2),sel(3),:));
odd_raw=squeeze(rawvolume(sel(1),sel(2),sel(4),:));
odd_corr=squeeze(corrvolume(sel(1),sel(2),sel(4),:));

% Create plot to see differences
figure(1)
subplot(2,1,1);
plot(1:131,even_raw,'b',1:131,even_corr,'r');
title('Voxel in even slice (40,15,14) before and after time correction');
legend('Original','Corrected');

subplot(2,1,2);
plot(1:131,odd_raw,'b',1:131,odd_corr,'r');
title('Voxel in odd slice (40,15,15) before and after time correction');
legend('Original','Corrected');

%% Q2: Detect head motion

front=squeeze(corrvolume(30,50,12,:)); 
back=squeeze(corrvolume(31,9,5,:));
left=squeeze(corrvolume(48,24,8,:));
right=squeeze(corrvolume(15,24,8,:));
top=squeeze(corrvolume(33,28,21,:));
bottom=squeeze(corrvolume(34,31,1,:));

figure(2)
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
o_fd=fopen('FUNC4_vol.011.img','rb');
a_fd=fopen('aFUNC4_vol.011.img','rb');
t_fd=fopen('r_tri_aFUNC4_vol.011.img','rb');
f_fd=fopen('r_4th_aFUNC4_vol.011.img','rb');
n_fd=fopen('r_nnb_aFUNC4_vol.011.img','rb');

N_slice=21;
N_x=64;
N_y=64;
for k=1:N_slice
    for i=1:N_x
        for  j=1:N_y
             o_volume(i,j,k)=fread(o_fd,1,'short');
             a_volume(i,j,k)=fread(a_fd,1,'short');
             t_volume(i,j,k)=fread(t_fd,1,'short');
             f_volume(i,j,k)=fread(f_fd,1,'short');
             n_volume(i,j,k)=fread(n_fd,1,'short');
        end
    end
end
fclose(fd);

figure('NumberTitle','off','Name','Different Reslice Options Axial View')
slice_no = 3;
subplot(2,3,1)
org_image=zeros(N_x,N_y);
org_image=a_volume(1:64,1:64,slice_no);
imagesc(org_image);
title('Time Corrected Volume')

subplot(2,3,2)
org_image=t_volume(1:64,1:64,slice_no);
imagesc(org_image);
title('Trilinear Interpolation')

subplot(2,3,3)
org_image_f=f_volume(1:64,1:64,slice_no);
imagesc(org_image_f);
title('4th Degree B-Spline Interpolation')

subplot(2,3,4)
org_image_n=n_volume(1:64,1:64,slice_no);
imagesc(org_image_n);
title('Nearest Neighbor Interpolation')

subplot(2,3,5)
org_image_o=o_volume(1:64,1:64,slice_no);
imagesc(org_image);
title('Original Image')

subplot(2,3,6)
imagea = org_image_f - org_image_o ;
imagesc(imagea)
colormap(gray);
title('Pixel Variation from original to 4th Degree interpolation')
%----------------------------------------------------
figure('NumberTitle','off','Name','Different Reslice Options Frontal View')
subplot(2,3,1)
Nx_slice = 50;
org_image=zeros(N_x,N_slice);
org_image_f=zeros(N_x,N_slice);
org_image_o=zeros(N_x,N_slice);
image=zeros(N_x,N_slice);
org_image=flip(rot90(squeeze(a_volume(Nx_slice,1:64,1:21)),3));
imagesc(org_image);
title('Time Corrected Volume')

subplot(2,3,2)
org_image=flip(rot90(squeeze(t_volume(Nx_slice,1:64,1:21)),3));
imagesc(org_image);
title('Trilinear Interpolation')

subplot(2,3,3)
org_image_f=flip(rot90(squeeze(f_volume(Nx_slice,1:64,1:21)),3));
imagesc(org_image_f);
title('4th Degree B-Spline Interpolation')

subplot(2,3,4)
org_image=flip(rot90(squeeze(n_volume(Nx_slice,1:64,1:21)),3));
imagesc(org_image);
title('Nearest Neighbor Interpolation')

subplot(2,3,5)
org_image_o=flip(rot90(squeeze(o_volume(Nx_slice,1:64,1:21)),3));
imagesc(org_image_o);
title('Original Image')

subplot(2,3,6)
imagef = (org_image_f - org_image_o);
imagesc(imagef)
colormap(gray);
title('Pixel Variation from original to 4th Degree interpolation')
%---------------------------------------------------
figure('NumberTitle','off','Name','Different Reslice Options Sagital View')
subplot(2,3,1)
Ny_slice = 30;
org_image=zeros(N_x,N_slice);
org_image_f=zeros(N_x,N_slice);
org_image_o=zeros(N_x,N_slice);
image=zeros(N_x,N_slice);
org_image=rot90(flip((squeeze(a_volume(1:64,Ny_slice,1:21)))));
imagesc(org_image);
title('Time Corrected Volume')

subplot(2,3,2)
org_image=rot90(flip((squeeze(a_volume(1:64,Ny_slice,1:21)))));
imagesc(org_image);
title('Trilinear Interpolation')

subplot(2,3,3)
org_image_f=rot90(flip((squeeze(a_volume(1:64,Ny_slice,1:21)))));
imagesc(org_image_f);
title('4th Degree B-Spline Interpolation')

subplot(2,3,4)
org_image=rot90(flip((squeeze(a_volume(1:64,Ny_slice,1:21)))));
imagesc(org_image);
title('Nearest Neighbor Interpolation')

subplot(2,3,5)
org_image_o=rot90(flip((squeeze(o_volume(1:64,Ny_slice,1:21)))));
imagesc(org_image_o);
title('Original Image')

subplot(2,3,6)
images = (org_image_f - org_image_o );
imagesc(images)
colormap(gray);
title('Pixel Variation from original to 4th Degree interpolation')
%---------------------
figure('NumberTitle','off','Name','Comparison of different interpolation methods')
slice_no = 3;

org_image_o=o_volume(1:64,1:64,slice_no);
org_image_a=zeros(N_x,N_y);
org_image_t=zeros(N_x,N_y);
org_image_f=zeros(N_x,N_y);
org_image_n=zeros(N_x,N_y);

subplot(2,2,1)

org_image_a=a_volume(1:64,1:64,slice_no);
image = org_image_a - org_image_o
imagesc(image);
title('Time Corrected Volume')

subplot(2,2,2)
org_image=t_volume(1:64,1:64,slice_no);
image = org_image_t - org_image_o
imagesc(image);
title('Trilinear Interpolation')

subplot(2,2,3)
org_image_f=f_volume(1:64,1:64,slice_no);
image = org_image_f - org_image_o
imagesc(image);
title('4th Degree B-Spline Interpolation')

subplot(2,2,4)
org_image_n=n_volume(1:64,1:64,slice_no);
image = org_image_n - org_image_o
imagesc(image);
title('Nearest Neighbor Interpolation')
colormap(gray);


% Refer to PDF for explanation