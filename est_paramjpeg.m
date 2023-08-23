function [coef p_ini]= est_paramjpeg(mu_hat,v_hat,e)

K=0.0833;
K_mode=1;
opts = optimset('MaxFunEvals',100000,'MaxIter',100000,'TolX',1e-8,'TolFun',1e-8);
    % Initial solution
    ind = 1;
    p_ini = zeros(length(mu_hat)*(length(mu_hat)+1)/2,2);
    for i = 1:length(mu_hat)
        for j = i+1:length(mu_hat)
            A=[mu_hat(i)^2 mu_hat(i); mu_hat(j)^2 mu_hat(j)];
%             B = [(v_hat(i)) ; (v_hat(j))];
            B = [(v_hat(i)-K ) ; (v_hat(j)-K )];
            p_ini(ind,:) = linsolve(A,B);
            ind = ind + 1;
        end
    end
%     nan_p_ini=find(p_ini(:,1)==NAN);
%     p_ini(nan_p_ini,:)=[];

mean_pa=mean(p_ini(:,1));
var_pa=var(p_ini(:,1));
mean_pb=mean(p_ini(:,2));
var_pb=var(p_ini(:,2));
indpa=find((abs(p_ini(:,1)-mean_pa) <= 3*sqrt(var_pa)) & (abs(p_ini(:,2)-mean_pb) <= 3*sqrt(var_pb)));
p1a=p_ini(indpa,1);
p1b=p_ini(indpa,2);
p_ini=p_ini(indpa,:);
    p_ini1 = mean(p_ini);
%     p_ini1=[p_ini1,K];
    switch K_mode
        case 0
            f = @(p) costfun_nonK(mu_hat,v_hat,p,e);
        case 1
            f = @(p) costfun_K(mu_hat,v_hat,p,e,K); 
        otherwise
            f = @(p) costfun_K(mu_hat,v_hat,p,e,K);      
    end
    coef = fminsearch(f,p_ini1,opts);
    
    
return

function E   = costfun_nonK(mu_hat,v_hat,p,e,K)

eps = 0.00000000000001;
fun = @(mu) max(p(1).*(mu.^2) + p(2).*mu,eps) ;
E = sum(log(e.*(fun(mu_hat).^2)) + ((v_hat-fun(mu_hat)).^2)./(e.*(fun(mu_hat).^2)) );
return

function E   = costfun_K(mu_hat,v_hat,p,e,K)

eps = 0.00000000000001;
fun = @(mu) max(p(1).*(mu.^2) + p(2).*mu+K  ,eps) ;
E = sum(log(e.*(fun(mu_hat).^2)) + ((v_hat-fun(mu_hat)).^2)./(e.*(fun(mu_hat).^2)) );%2*pi*
return



