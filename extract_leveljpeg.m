function [Set mu_hat v_hat zk num d e] = extract_leveljpeg(z,z_approx,z_res)

% Segmentation
z_smt=z_approx;
img = z; img = img(:);
img_approx = z_approx; img_approx = img_approx(:);
img_res = z_res; img_res = img_res(:);
img_smt = z_smt; img_smt = img_smt(:);

u = linspace(0,255,256);
delta = u(2) - u(1);
Set = zeros(numel(img),4);
ind = 1;
% for i = floor(0.1*N):delta:u(end)
for i = 25:length(u)
    tmp = find((img_smt>=u(i)-delta/2) & (img_smt<u(i)+delta/2));
    if size(tmp,1) >= 100 
        Set(ind:ind+size(tmp,1)-1,:) = [ i*ones(size(tmp)) img_approx(tmp) img_res(tmp)  img(tmp) ];
        ind = ind + size(tmp,1);
    end
end
Set = Set(Set(:,1) > 0,:); clear tmp;

mu_hat = zeros(length(u),1);
v_hat = zeros(length(u),1);
num = zeros(length(u),1);
zk = zeros(length(u),1);
d = zeros(length(u),1);
e = zeros(length(u),1);
tau = 3;

for i = 1:length(u)
   ind = find(Set(:,1) == i); 
   tmp = Set(ind,2:4);   
   if size(tmp,1) >= 2 
        m1 = median(tmp(:,1));
        s1 = 1.4826*mad(abs(tmp(:,1)));
        s2 = 1.4826*mad(tmp(:,2),1);   
        index1 = find( abs(tmp(:,1) - m1) > tau*s1 | abs(tmp(:,2)) > tau*s2 );
        Set(ind(index1),:) = 0;
        while numel(index1) > 0 && s2 ~= 0 
            index2 = find(abs(tmp(:,1) - m1) <= tau*s1 & abs(tmp(:,2)) <= tau*s2 );
            tmp = tmp(index2,:);
            ind = ind(index2);     
            m1 = median(tmp(:,1));
            s1 = 1.4826*mad(abs(tmp(:,1)));
            s2 = 1.4826*mad(tmp(:,2),1);       
            index1 = find( abs(tmp(:,1) - m1) > tau*s1 |  abs(tmp(:,2)) >= tau*s2 );
            Set(ind(index1),:) = 0;
        end
        
       if size(tmp,1) >= 100
           num(i) = size(tmp,1);           
           mu_hat(i) = mean(tmp(:,1));
           v_hat(i) = var(tmp(:,2));           
           %v_hat(i) = (1.4826*mad(tmp(:,2),1))^2;
           zk(i) = mean(tmp(:,3));
           d(i) = 1/num(i);
           e(i) = 2/(num(i));
           
       else
           Set(Set(:,1) == i,:) = 0;
       end
       
   end
end
Set = Set(Set(:,1) > 0,:);
index = find(mu_hat > 0 );
mu_hat = mu_hat(index);
v_hat = v_hat(index);
zk = zk(index);
num = num(index);
d = d(index);
e = e(index);
%plot(mu_hat,v_hat,'x')
clear z z_res z_smt z_approx z_lap img img_res img_approx img_smt Cl X;
end