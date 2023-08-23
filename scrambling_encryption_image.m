function [e_z,e_z_approx,e_z_res] = scrambling_encryption_image(encryption,z,z_approx,z_res)

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
 
else
    e_z=z;
    e_z_approx=z_approx;
    e_z_res=z_res;
    
end