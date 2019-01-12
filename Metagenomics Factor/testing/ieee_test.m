function [outputArg1,outputArg2] = ieee_test(A,U,V)
%testing function for algorithm found in ieee paper

[m,n] = size(A);
[~,p] = size(U); %do rank finding later
%{
sigma_u = (p-1)/(p^2*(p+1));
u = sqrt(1/p)*ones(m,p); S = sigma_u*ones(m,p);
v = 1/4*ones(p,n); C = (3/16)*ones(p,n);
%}

u = zeros(m,p) + 1/m; S = ones(m,p);
v = zeros(p,n); C = ones(p,n);

w = sqrt(1/p)*u*v;
x = (1/p)*S*C + (u.^2)*C + S*(v.^2);

invsigma = (1/p)*((x.^(-1))'*(u.^2 + S) - (((A - w)./x).^2)'*S)';
sigma = invsigma.^(-1);
T = sigma.*(sqrt(1/p)*(((A - w)./x)'*u)' + v.*((1/p)*((x.^(-1))'*(u.^2))'));

invZ = (1/p)*(x.^(-1)*((v.^2 + C)') - (((A - w)./x).^2)*C');
Z = invZ.^(-1);
W = Z.*(sqrt(1/p)*(((A - w)./x)*v') + u.*((1/p)*((x.^(-1))*(v.^2'))));

old_v = v;
v = mean_normal_denoise(sigma,T,1,0);
C = var_normal_denoise(sigma,T,1);

old_u = u;
u = mean_normal_denoise(Z,W,1,0);
S = var_normal_denoise(Z,W,1);


old_w = w;
old_x = x;
w = sqrt(1/p)*u*v - ((A - w)./x).*((u.*old_u)*C + S*(v.*old_v));
x = (1/p)*S*C + (u.^2)*C + S*(v.^2);

invsigma = (1/p)*((x.^(-1))'*(u.^2 + S) - (((A - w)./x).^2)'*S)';
sigma = invsigma.^(-1);
T = sigma.*(sqrt(1/p)*(((A - w)./x)'*u)' + v.*((1/p)*((x.^(-1))'*(u.^2))')) - old_v.*((1/p)*S'*(((A - w)./x).*((A - old_w)./old_x)));

invZ = (1/p)*(x.^(-1)*((v.^2 + C)') - (((A - w)./x).^2)*C');
Z = invZ.^(-1);
W = Z.*(sqrt(1/p)*(((A - w)./x)*v') + u.*((1/p)*((x.^(-1))*(v.^2')))) -  old_u.*((1/p)*(((A - w)./x).*((A - old_w)./old_x))*C');

v = ieee_Vdenoise(sigma,T);
C = v.*(1-v);

u = ieee_Udenoise(Z,W);
S = (sigma_u*Z)./(sigma_u + Z);
end

