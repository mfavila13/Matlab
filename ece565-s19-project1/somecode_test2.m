clc;
clear all;
close all;
i1=imread('test2.tif');
i1=rgb2gray(i1);

f1=fftn(i1);
mag1=abs(f1);
s=log(1+fftshift(f1));
phase1=angle(f1);

r1=ifftshift(ifftn(mag1));
r2=ifftn(exp(1i*phase1));
figure,imshow(i1);
figure,imshow(s,[]);
figure,imshow(uint8(r1));
figure,imshow(r2,[]);
r2=histeq(r2);
r3=histeq(uint8(r2));     
figure,imshow(r2);
figure,imshow(r3);