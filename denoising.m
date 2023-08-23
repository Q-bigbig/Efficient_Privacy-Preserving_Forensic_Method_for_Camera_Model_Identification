function  [z,z_approx,z_res] = denoising(z,sigma)

z_res = noiseextractfromgrayscale(z,sigma); %z_res = single(z_res);
%z_res = double(WienerInDFT(z_res,std2(z_res)));
z_approx = z - z_res;
end