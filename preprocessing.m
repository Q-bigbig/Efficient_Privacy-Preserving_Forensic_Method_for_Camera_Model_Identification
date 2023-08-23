function [B,Cl,X,BCX]=preprocessing(z,z_approx,z_res)
N=256;
%Clipping
Cl = z <= 0.1*N | z >= 0.9*N;
Cl = im2vec(Cl,8);
Cl = sum(Cl); 
Cl = Cl';
Cl = reshape(Cl,size(z,1)/8,size(z,2)/8);
Cl = kron(Cl,ones(8)); 
% Cl = Cl(:);

% Uniform block
B = im2vec(z,8);
B = mean(B); 
B = repmat(B,64,1);
B = im2vec(z,8) - B;
B = B ~= 0; 
B = sum(B); B = B';

B = reshape(B,size(z,1)/8,size(z,2)/8);
B = kron(B,ones(8)); 
% B = B(:);


% Edge detection
sig = 1.4826*mad(z_res(:),1);
x = bdct(z_approx);
x = im2vec(x,8);
x = x(2:64,:);
x = 1.4826*mad(x,1);
X = x > sig; 
y = x(X==0);
if numel(y) > 1 %monia: avant ça était if numel(y) > 0 (j'ai changé 0 en 1): si y contient un seul element donc numel(y)=1 ->floor(0.8)=0 ->Y(0) ->erreur "??? Attempted to access Y(0); index must be a positive integer or logical."
Y = sort(y); 
seuil = Y(floor(0.9*numel(y)));
X = x > seuil;
end
% BW = reshape(X,size(z,1)/8,size(z,2)/8);
% BW = kron(BW,ones(8)); 
% X = X(:);

    X=X';
    X = reshape(X,size(z_approx,1)/8,size(z_approx,2)/8);
    X = kron(X,ones(8));
% XX = reshape(X,size(z,1)/8,size(z,2)/8);
% XX = kron(XX,ones(8)); XX = XX(:);

BCX=(X == 0 & Cl == 0 & B~=0);
end