function [e_z,e_z_approx,e_z_res] = revision_pure_scrambling_encryption_image(encryption,z,z_approx,z_res)

if encryption
    [M,N]   = size(z);
    Num=M*N;
    squM=randsample(Num,Num)';
    squMtmp=[1:1:Num;squM];
    
    z0=z(:);
    z_approx0=z_approx(:);
    z_res0=z_res(:);
    
    e_z0 (squMtmp(1,:)) = z0 (squMtmp(2,:));
    e_z_approx0 (squMtmp(1,:)) = z_approx0 (squMtmp(2,:));
    e_z_res0 (squMtmp(1,:)) = z_res0 (squMtmp(2,:));
    
    e_z=reshape(e_z0,M,N);
    e_z_approx=reshape(e_z_approx0,M,N);
    e_z_res=reshape(e_z_res0,M,N);
    

    figure;
    subplot(2,3,1);imshow(z./255); title('inquiry image');hold on;
    subplot(2,3,4);imshow(e_z./255); title('encrypted image');hold on;
    subplot(2,3,2);imshow(z_approx./255);title('denoised image');hold on;
    subplot(2,3,5);imshow(e_z_approx./255);title('encrypted denoised image');hold on;
    subplot(2,3,3);imshow(z_res);title('noise');hold on;
    subplot(2,3,6);imshow(e_z_res);title('encrypted noise');hold on;
else
    e_z=z;
    e_z_approx=z_approx;
    e_z_res=z_res;
    
%     figure;
%     subplot(2,3,1);imshow(z./255); title('z');hold on;
%     subplot(2,3,4);imshow(e_z./255); title('e_z without encryption');hold on;
%     subplot(2,3,2);imshow(z_approx./255);title('z_approx');hold on;
%     subplot(2,3,5);imshow(e_z_approx./255);title('e_z_approx without encryption');hold on;
%     subplot(2,3,3);imshow(z_res);title('z_res');hold on;
%     subplot(2,3,6);imshow(e_z_res);title('e_z_res without encryption');hold on;
end