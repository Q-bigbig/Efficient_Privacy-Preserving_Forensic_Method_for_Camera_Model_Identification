clc
clear all
close all

addpath('D:\YanliCHEN\MATLAB\CODE\encryption_camera_attribution_test\data\fingerprints');
imagefiles = dir('D:\YanliCHEN\MATLAB\CODE\encryption_camera_attribution_test\data\fingerprints\*.jpg');  
image_numfiles = length(imagefiles);
im_num=image_numfiles;

channel=1;
scrambling_encryption=1;
noise_encryption= 1;
scale_res=1.25;

tic
for i=1:im_num
    image_name = imagefiles(i).name;
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
     [z,z_approx,z_res] = noise_encryption_image(noise_encryption,z,z_approx,z_res,scale_res);
     [e_z,e_z_approx,e_z_res] = scrambling_encryption_image(scrambling_encryption,z,z_approx,z_res);

%      %fingerprints extraction
     [Set mu_hat v_hat zk num d e] = extract_leveljpeg(e_z,e_z_approx,e_z_res);  
     index = find( abs(zk-mu_hat) <= 3*sqrt(v_hat./num));     
     if numel(index) >= 10
            [coef_ab p]= est_paramjpeg(mu_hat(index),v_hat(index),e(index));
     else
         coef_ab=[100 100];
     end
     
     coef(i,:)=coef_ab;  

end
toc
load chirp
sound(y,Fs)

 
 