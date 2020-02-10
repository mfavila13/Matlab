% Marcus Favila
% problem 1. c.
% use of MATLAB histeq() of test image, 
% displays the grayscale image 
img = imread('test1.tif');
H = histeq(img);
figure('NumberTitle', 'off', 'Name', 'Use of histeq()')
imshow(H)