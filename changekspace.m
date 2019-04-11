%this program allows one to add a spike in k-space at a desired location
%usage: changekspace('FUNC3_vol.002.img',2,2,3)
%first number in the parentheses is x-location in k-space (0-32)
%second number in the parentheses is y-location in k-space (0-32)
%third number in the parentheses is the slice number of the volume.
%Center of k-space is at (x,y)=(0,0);

function changekspace(s,x,y,slice_no)

if ( (x>(32))||(y>32)||(x<(-32))||(y<(-32)))
    disp('invalid argument: x and y should be in the range (-32:32)');
    return
else

fd=fopen(s,'rb');


%For FUNC3 volumes
N_slice=7;


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

org_kspace=zeros(N_x,N_y);

org_kspace=(fftshift(fft2(org_image)));

mod_kspace=org_kspace;
%%% vertical line 
mod_kspace(1:16,33)=50e2;
mod_kspace(48:64,33)=50e2;
mod_kspace(33,33) = org_kspace(33,33);

mod_kspace(33-y,x+33)=50e2;

mod_kspace(33+y,-x+33)=50e2;

mod_image=100+(ifft2(mod_kspace));

maximum= max(max(abs(mod_image)));

figure;

imagesc(abs(org_image));

colormap(gray);

caxis([0 maximum]);

title('Original Image');

figure;

imagesc(abs(mod_image));

colormap(gray);

title('Modified Image');

caxis([0 maximum]);

figure;

imagesc(20*log10(abs(org_kspace)));

colormap(gray);

title('Original K-space');

figure;

imagesc(20*log10(abs(mod_kspace)));

colormap(gray);

title('Modified K-space');

end
