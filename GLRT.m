clc
clear all
close all

addpath('D:\YanliCHEN\MATLAB\CODE\encryption_camera_attribution_test\data');
imagefiles = dir('D:\YanliCHEN\MATLAB\CODE\encryption_camera_attribution_test\data\*.jpg');  
image_numfiles = length(imagefiles);
im_num=image_numfiles;

load('coef_D70.mat');
%load('ee_coef_D70.mat');
p0=mean(coef(:,:));

channel=1;
scrambling_encryption= 0;
noise_encryption= 0;
scale_res=1.25;

alpha=0.00001;
tau=norminv(1-alpha);

tic
for i=1:im_num
    image_name = imagefiles((i)).name;
    im = imread(image_name);
    
    z0 = double(im(:,:,channel));
    [r c] = size(z0); 
    z0 = z0(1:8*floor(r/8),1:8*floor(c/8)); 
    clear r c;
     
    sigma=4*function_stdEst2D(z0,4);
    [z,z_approx,z_res] = denoising(z0,sigma);
     
      %preprocessing
    [B,Cl,X,BCX]=preprocessing(z,z_approx,z_res);
    [z z_approx z_res B Cl X BCX]=preprocessing_SelectCenter(z,z_approx,z_res,B,Cl,X,BCX);
     z=z.*BCX;
     z_approx=z_approx.*BCX;
     z_res=z_res.*BCX;
      
      %noise encryption
     [e1_z,e1_z_approx,e1_z_res] = scrambling_encryption_image(scrambling_encryption,z,z_approx,z_res);
     [e_z,e_z_approx,e_z_res] = noise_encryption_image(noise_encryption,e1_z,e1_z_approx,e1_z_res,scale_res);

     %GLRT
     LambN2(i,:)= encryption_GLRT(e_z,e_z_res,e_z_approx,p0,coef);

end
toc
load chirp
sound(y,Fs)
