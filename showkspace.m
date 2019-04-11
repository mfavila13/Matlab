%this program displays the k-space and image-space representation of 
%a particular slice of a volume.
%usage: showkspace('FUNC2_vol.006.img',10)

function showkspace(s,slice_no)

fd=fopen(s,'rb');
%for FUNC2 volume
N_slice=14;

%for FUNC3 volume
%N_slice=7;
N_x=64;

N_y=64;

volume=zeros(N_x,N_y,N_slice); % volume is a 3-D array

for k=1:N_slice

    for i=1:N_x
 
        for  j=1:N_y

             volume(i,j,k)=fread(fd,1,'short');
        
        end

    end

end

fclose(fd);

org_image=zeros(N_x,N_y);

org_image=volume(1:end,1:end,slice_no);

figure;

imagesc(org_image);

colormap(gray);

title('Image');

kspace=zeros(N_x,N_y);

kspace=fftshift(fft2(org_image));


figure;

imagesc(20*log10(abs(kspace)));

colormap(gray);

title('K-space representation in log10 scale');