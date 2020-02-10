img = imread('test2.tif');
FT_img = fft2(img);
FT_imgR = real(FT_img);
FT_imgI = imag(FT_img);
inv_img = uint8(ifft2(FT_img));
M_spec = abs(((FT_imgR)^2 + (FT_imgI^2))^(1/2));
P_angle = atan(FT_imgI/FT_imgR);
P_conj = conj((P_angle));

M_img = uint8(ifft2(M_spec));
P_img = uint8(ifft2(P_angle));
M_P_img = uint8(ifft2(M_spec*exp(1j*P_angle)));
M_Pconj_img = uint8(ifft2(M_spec*exp(1j*P_conj)));

figure('NumberTitle', 'off', 'Name', 'Test 2');
subplot(6,2,1)
imshow(img);title('Original Image');
subplot(6,2,2)
imshow(FT_img);title('FT of Original Image');

subplot(6,2,3)
imshow(M_img);title('Magnitude of FT image');

subplot(6,2,4)
imshow(P_img);title('Phase of FT image');

subplot(6,2,5);
imshow(M_P_img);title('Image obtained from Magnitude and Phase angle');

subplot(6,2,6);
imshow(M_Pconj_img);title('Image obtained from Magnitude and Phase angle conjugate');


