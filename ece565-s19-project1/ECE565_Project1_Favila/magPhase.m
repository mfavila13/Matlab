% Marcus Favila
% Code that calculates:
% - FT of test2.tif
% - Magnitude spectrum of FT of test2.tif
% - Phase angle of FT of test2.tif
% - displays various plots obtained from inverse FT of components of
%   the FT of test2.tif

A = imread('test2.tif');
img = complex(A);
F = fft2(img);

for i = 1:size(A,1)
    for k = 1:size(A,2)
        M_ft(i,k) = abs(F(i,k));
        theta_ft(i,k) = angle(F(i,k));
        conj_theta(i,k) = conj(theta_ft(i,k));
        phase_comp(i,k) = exp(1j*theta_ft(i,k));
        phase_comp_conj(i,k) = conj(exp(1j*theta_ft(i,k))); 
        MP_img(i,k) = M_ft(i,k) * phase_comp(i,k);
        MPconj(i,k) = M_ft(i,k) * phase_comp_conj(i,k);

    end
end
img_M = uint8(ifft2(M_ft));
img_P = uint8(ifft2(phase_comp));
img_MP = uint8(ifft2(MP_img));
img_MPconj = uint8(ifft2(MPconj));


figure('NumberTitle', 'off', 'Name', 'Problem 2.')
subplot(3,2,1)
imshow(A);title('Orginal Image');

subplot(3,2,2)
imshow(F);title('FT of Original Image');

subplot(3,2,3)
imshow(img_M);title('Image: inverse FT of magnitude');

subplot(3,2,4)
imshow(img_P);title('Image: inverse FT of phase term ');

subplot(3,2,5)
imshow(img_MP);title('Image: inverse FT of Mag and Phase Angle');

subplot(3,2,6)
imshow(img_MPconj);title('Image: inverse FT of Mag and Phase Angle conjugate');

