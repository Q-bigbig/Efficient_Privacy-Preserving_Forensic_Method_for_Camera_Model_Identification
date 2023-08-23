function [zs z_approxs z_ress Bs Cls Xs BCXs]=preprocessing_SelectCenter(z,z_approx,z_res,B,Cl,X,BCX)

m=4;
n=4;
offset=1;
bsize=8;
switch offset
    case 0 
        %output for later used
        zs = z(m:bsize:end,n:bsize:end); 
        z_approxs = z_approx(m:bsize:end,n:bsize:end); 
        z_ress = z_res(m:bsize:end,n:bsize:end); 
        
        Bs = B(m:bsize:end,n:bsize:end); 
        Cls = Cl(m:bsize:end,n:bsize:end); 
        Xs = X(m:bsize:end,n:bsize:end); 
        BCXs = BCX(m:bsize:end,n:bsize:end); 
        
    case 4
         %output for later used
        zs = z(1:end,1:end);
        z_approxs = z_approx(1:end,1:end);
        z_ress = z_res(1:end,1:end);
        
        Bs = B(1:end,1:end);
        Cls = Cl(1:end,1:end);
        Xs = X(1:end,1:end);
        BCXs = BCX(1:end,1:end);
        
    otherwise %only 1 2 3
        zs=zeros(size(z,1)/bsize*(2*offset),size(z,2)/bsize*(2*offset));
        z_approxs=zeros(size(z,1)/bsize*(2*offset),size(z,2)/bsize*(2*offset));
        z_ress=zeros(size(z,1)/bsize*(2*offset),size(z,2)/bsize*(2*offset));
        
        Bs=zeros(size(z,1)/bsize*(2*offset),size(z,2)/bsize*(2*offset));
        Cls=zeros(size(z,1)/bsize*(2*offset),size(z,2)/bsize*(2*offset));
        Xs=zeros(size(z,1)/bsize*(2*offset),size(z,2)/bsize*(2*offset));
        BCXs=zeros(size(z,1)/bsize*(2*offset),size(z,2)/bsize*(2*offset));
%         tic
        for ci=1:(2*offset)
            for ri=1:(2*offset)
                zs(ri:(2*offset):end,ci:(2*offset):end)=z(m-offset+ri:bsize:end,n-offset+ci:bsize:end);
                z_approxs(ri:(2*offset):end,ci:(2*offset):end)=z_approx(m-offset+ri:bsize:end,n-offset+ci:bsize:end);
                z_ress(ri:(2*offset):end,ci:(2*offset):end)=z_res(m-offset+ri:bsize:end,n-offset+ci:bsize:end);
                
                Bs(ri:(2*offset):end,ci:(2*offset):end)=B(m-offset+ri:bsize:end,n-offset+ci:bsize:end);
                Cls(ri:(2*offset):end,ci:(2*offset):end)=Cl(m-offset+ri:bsize:end,n-offset+ci:bsize:end);
                Xs(ri:(2*offset):end,ci:(2*offset):end)=X(m-offset+ri:bsize:end,n-offset+ci:bsize:end);
                BCXs(ri:(2*offset):end,ci:(2*offset):end)=BCX(m-offset+ri:bsize:end,n-offset+ci:bsize:end);
            end
            
        end
%         toc            
end
end


