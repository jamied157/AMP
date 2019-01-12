function [ X ] = gbm_factor( A , q  , ell)
%GBM_FACTOR The amp algorithm to reconstruct the factors of the 
%gaussian block model
n = length(A(1,:));
[V,~] = eigs(A,q-1);
X_old = [sqrt(n)*V,zeros(n,1)];X_older = [sqrt(n)*V,zeros(n,1)];
for i = 1:20
    X = A*gbm_nl(X_old,ell) - gbm_nl(X_older,ell)*mean(gbm_dnl(X_old,ell),3)';
    X_older = X_old;
    X_old = X;
end
end

