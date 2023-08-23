function LambN2= encryption_GLRT(z,z_res,z_approx,p0,coef)

[Set mu_hat v_hat zk num d e] = extract_leveljpeg(z,z_approx,z_res);        
u = linspace(0,255,256);
% [Set mu_hat v_hat zk num d e] = extract_leveljpeg_CYL(z,z_res,z_approx,X,m,n,N,levelNum,TH,offset);
index0 = find(abs(zk-mu_hat) <= 3*sqrt(v_hat./num) );
% numel(index0) 
if numel(index0) >= 10
[c p_ini]= est_paramjpeg(mu_hat(index0),v_hat(index0),e(index0));
p1_hat = [c(1) c(2)];
% clear index0;

% % fun = @(p,mu) max(p(1).*(mu.^(2)) + p(2).*(mu) + K ,eps);
% %     t = 1:0.1:255;
% %     figure;
% %     hold on; plot(mu_hat(index0),v_hat(index0),'x');
% %     hold on; plot(t,fun(p1_hat,t),'LineWidth',2);  
    
% mean_pa=mean(p_ini(:,1));
% var_pa=var(p_ini(:,1));
% mean_pb=mean(p_ini(:,2));
% var_pb=var(p_ini(:,2));
% indpa=find((abs(p_ini(:,1)-mean_pa) <= 3*sqrt(var_pa))& (abs(p_ini(:,2)-mean_pb) <= 3*sqrt(var_pb)));
% p1a=p_ini(indpa,1);
% p1b=p_ini(indpa,2);

% covmat = covariance_matrixCYL(p0,mu_hat,e,K);
% covmat = cov(p1a(:,1),p1b(:,1));
covmat = cov(p_ini(:,1),p_ini(:,2));
% covmat = cov(coef(:,1),coef(:,2));
va = covmat(1,1);
vb = covmat(2,2);
covab = covmat(1,2);

pix = Set(:,4);
ind = 1;
ind1 = 1;
tmpmu = zeros(size(pix));
tmpv = zeros(size(pix));
tmpnum = zeros(size(pix));
tmpzk = zeros(size(pix));
for i = 1:length(u)
    x = find(Set(:,1) == i);
    if numel(x) >= 2
    tmpmu(ind:ind+numel(x)-1) = mu_hat(ind1)*ones(numel(x),1);
    tmpv(ind:ind+numel(x)-1) = v_hat(ind1)*ones(numel(x),1);
    tmpnum(ind:ind+numel(x)-1) = num(ind1)*ones(numel(x),1);
    tmpzk(ind:ind+numel(x)-1) = zk(ind1)*ones(numel(x),1);
    ind = ind + numel(x);
    ind1 = ind1 + 1;
    end
end
clear mu_hat v_hat num zk;
mu_hat = tmpmu; v_hat = tmpv; num = tmpnum; zk = tmpzk;
clear tmpmu tmpv tmpnum x ind ind1 tmpzk;
 
index = find(abs(zk-mu_hat) <= 3*sqrt(v_hat./num) );
if numel(index) > 50
pix = pix(index);    
mu_hat = mu_hat(index);
v_hat = v_hat(index);
zk = zk(index);
num = num(index);

indexp = find(abs(pix-mu_hat) <= 3*sqrt(v_hat) );
if numel(indexp) > 2
pix = pix(indexp);    
mu_hat = mu_hat(indexp);
v_hat = v_hat(indexp);
zk = zk(indexp);
num = num(indexp);

fun = @(p,mu) max(p(1).*(mu.^2) + p(2).*mu +1/12 ,eps) ;

v0 = fun(p0,mu_hat);
v1 = fun(p1_hat,mu_hat);
h1 = log(v0./v1);
h2 = 1./v0 - 1./v1;

varf = mu_hat.^4 * va + mu_hat.^2 * vb + 2* mu_hat.^3 * covab;


rhok = (zk-mu_hat).^2;
Lambk = 1/2 * h1 +  num/2 .* h2.*rhok  ;
E0k = 1/2 * h1 +  1/2 * h2.*v0 ;
V0k = 1/2 .*(h2.^2).*(v0.^2);
E1k = 1/2 * h1 +  1/2 * h2.*v1 ;
V1k = 1/2 .*(h2.^2).*(v1.^2);
% V1k_b = 1/4 * varf./(v1.^2) + 3/4 * varf./(v1.^4) .* (v0.^2) + 1/2 .*(h2.^2).*(v0.^2);%original
V1k_b = 1/4 * varf./(v1.^2) + 3/4 * varf./(v1.^4) .* (v1.^2) + 1/2 .*(h2.^2).*(v1.^2);%cyl
V0k_b = 1/4 * varf./(v1.^2) + 3/4 * varf./(v1.^4) .* (v0.^2) + 1/2 .*(h2.^2).*(v0.^2);

rho = (pix-mu_hat).^2;
Lamb = 1/2 * h1 +  1/2 .* h2.*rho  ;
E0 = 1/2 * h1 +  1/2 * h2.*v0 ;
V0 = 1/2 .*(h2.^2).*(v0.^2);
E1 = 1/2 * h1 +  1/2 * h2.*v1 ;
V1 = 1/2 .*(h2.^2).*(v1.^2);
V0_b = 1/4 * varf./(v1.^2) + 3/4 * varf./(v1.^4) .* (v0.^2) + 1/2 .*(h2.^2).*(v0.^2);
V1_b = 1/4 * varf./(v1.^2) + 3/4 * varf./(v1.^4) .* (v1.^2) + 1/2 .*(h2.^2).*(v1.^2);


id = find( abs(Lambk - E1k) <= 3*sqrt(V1k));
if(size(id,1)>50)
x = abs(Lamb(id) - E1(id));
s = sort(x);
id1 = find(x <= s(round(0.2*numel(x)))); 
E0 = E0(id(id1));
V0 = V0(id(id1));
V0_b = V0_b(id(id1));
Lamb = Lamb(id(id1));
Np = min(5000,numel(id1));
tmp = randperm(numel(id(id1))); tmp = tmp';
LambN = (sum(Lamb(tmp(1:Np))) - sum(E0(tmp(1:Np))))/sqrt(sum(V0(tmp(1:Np))));
LambN2 = (sum(Lamb(tmp(1:Np))) - sum(E0(tmp(1:Np))))/sqrt(sum(V0_b(tmp(1:Np))));


% E0 = E0(id(id1));
% V0 = V0(id(id1));
% V0_b = V0_b(id(id1));
% Lamb = Lamb(id(id1));
% % Np = min(5000,numel(id1));
% % tmp = randperm(numel(id(id1))); tmp = tmp';
% % LambN = (sum(Lamb(tmp(1:Np))) - sum(E0(tmp(1:Np))))/sqrt(sum(V0(tmp(1:Np))));
% % LambN2 = (sum(Lamb(tmp(1:Np))) - sum(E0(tmp(1:Np))))/sqrt(sum(V0_b(tmp(1:Np))));
% 
% LambN = (sum(Lamb) - sum(E0))/sqrt(sum(V0));
% LambN2 = (sum(Lamb) - sum(E0))/sqrt(sum(V0_b));
% % Np = min(5000,numel(E0));
% E0 = E0(id);
% V0 = V0(id);
% V0_b = V0_b(id);
% Lamb = Lamb(id);
% Np = min(5000,numel(id));
% tmp = randperm(numel(id)); tmp = tmp'; 
   
% LambN = (sum(Lamb(tmp(1:Np))) - sum(E0(tmp(1:Np))))/sqrt(sum(V0(tmp(1:Np))));
% LambN2 = (sum(Lamb(tmp(1:Np))) - sum(E0(tmp(1:Np))))/sqrt(sum(V0_b(tmp(1:Np))));
% % Np = min(5000,numel(E0));
% % tmp = randperm(Np); tmp = tmp'; 
% % LambN = (sum(Lamb(tmp(1:Np))) - sum(E0(tmp(1:Np))))/sqrt(sum(V0(tmp(1:Np))));
% % LambN2 = (sum(Lamb(tmp(1:Np))) - sum(E0(tmp(1:Np))))/sqrt(sum(V0_b(tmp(1:Np))));

% if numel(id1) >= 5000
%     Lamb = Lamb(id(id1));
%     E0 = E0(id(id1));
%     V0 = V0(id(id1));
%     V0_b = V0_b(id(id1));
%     Np = 5000;
%     %Np = [50 100 200 500 1000 2000 numel(id1)] ;
%     tmp = randperm(numel(id(id1))); tmp = tmp';
%     for i = 1:numel(Np)        
%         LambN(i) = (sum(Lamb(tmp(1:Np(i)))) - sum(E0(tmp(1:Np(i)))))/sqrt(sum(V0(tmp(1:Np(i)))));
%         LambN2(i) = (sum(Lamb(tmp(1:Np(i)))) - sum(E0(tmp(1:Np(i)))))/sqrt(sum(V0_b(tmp(1:Np(i)))));
%     end
% else
%    LambN = zeros(numel(Np),1);
%    LambN2 = zeros(numel(Np),1);
% end


else 
    LambN = 0;
    LambN2 = 0;
    p1_hat = [0 0];
end

else 
    LambN = 0;
    LambN2 = 0;
    p1_hat = [0 0];
end

else 
    LambN = 0;
    LambN2 = 0;
    p1_hat = [0 0];
end

else 
    LambN = 0;
    LambN2 = 0;
    p1_hat = [0 0];
end
end