function [outputArg1,outputArg2] = U_MCMC(U , M , S)
%U_MCMC is a nonlinearity designed after the optimal nonlinearity given by
%state evolution, it uses importance sampling
%UNFINISHED
[m,p] = size(U);
X = randg(1,m,p,1000);
X = X./sum(X,2);

W = sqrt(1/2*pi*det(S))*exp(-(1/2)*(M*X - U)*S^(-1)*(M*X - U));
end

