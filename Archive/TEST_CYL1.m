%%%%%%%%%%%%% SIMULATION %%%%%%%%%%%%%%%%
clear all
% close all
clc
%% Loading image
ImNum = 200;
y = imread('image_t1_dev_00001.jpg');
% figure; 
% imshow(y);
z=rgb2gray(y);
z=im2double(z);
figure; 
imshow(z);

%% Mu,a,b unknown
a0 = -0.0012 ; b0 = 0.11;
% a1 = -0.0025 ; b1 = 0.20;

% a0 = 4.5e-5*(4095/255)^2 ; b0 = 5e-7*(4095/255)^2;
% a1 = 7.5e-5*(4095/255)^2 ; b1 = 1e-6*(4095/255)^2;
m=2;
n=2;
N=256;
K=0.083;%K=deta.^2/12;

gamma=1;
figure;
for i = 1:ImNum
%     z = z + sqrt(a0*y_s2+b0*z)./255.*randn(size(z));
%     figure; 
%     imshow(z);
    
    sigma=function_stdEst2D(z);
    
%     [Set BW mu_hat v_hat zk num d e] = extract_leveljpeg_simu(z,m,n,sigma,N);
    [Set BW mu_hat v_hat zk num d e] = extract_leveljpeg(z,m,n,sigma,N);
    
%     hold on;
    plot(mu_hat(:),v_hat(:),'x')
%     coef = est_paramjpeg(mu_hat,v_hat,e,K,gamma);

    
    
% %     [a_hat0 b_hat0 var_a_hat0 var_b_hat0 cov0] = cam_param_est(im0);
% %     [a_hat1 b_hat1 var_a_hat1 var_b_hat1 cov1] = cam_param_est(im1);
% %     a2(i)=a_hat0;
% %     b2(i)=b_hat0;
% %     GLR0(i,:) = GLR_mu_a_b_unknown(im0,a0,b0,a_hat0, b_hat0, var_a_hat0, var_b_hat0, cov0);
% %     GLR1(i,:) = GLR_mu_a_b_unknown(im1,a0,b0,
% a_hat1, b_hat1, var_a_hat1, var_b_hat1, cov1);
end

% Plot ROC
beta = zeros(N,4);
for j = 1:4
    tmp = sort(GLR0(:,j));
    for i = 1:N
        beta(i,j) = sum(GLR1(:,j) >= tmp(i))/N;
    end
end

vecx = 1-(0:N-1)./N;
figure;
hold on
plot(vecx,beta(:,1),'g','linewidth',2)
plot(vecx,beta(:,2),'r','linewidth',2)
plot(vecx,beta(:,3),'m','linewidth',2)
plot(vecx,beta(:,4),'y','linewidth',2)
legend('50','100','200','500')

figure;
plot(a2(1,:),b2(1,:),'g*');


 %  fun = @(p,mu) max(p(1)./(p(3)^2).*(mu.^(2-p(3))) + p(2)./(p(3)^2).*(mu.^(2-2*p(3))) + K ,eps) ;
 %  t = 10:0.1:220;
 %  figure;
 %  hold on; plot(mu_hat,v_hat,'x')
 %  hold on; plot(t,fun(coef,t),'g')