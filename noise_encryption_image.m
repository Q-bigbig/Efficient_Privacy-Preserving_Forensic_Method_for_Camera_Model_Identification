function [e_z,e_z_approx,e_z_res] = noise_encryption_image(encryption,z,z_approx,z_res,scale_res)

if encryption
    e_z_res=scale_res.*z_res;
    e_z_approx=1*z_approx;
else
    e_z_res=z_res;
    e_z_approx=z_approx;
end
    e_z=e_z_approx+e_z_res;
end