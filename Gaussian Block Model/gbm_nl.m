function [ output ] = gbm_nl( X , ell)
%GBM_NL the non linearity function used in the amp factorisation algotithm of the gaussian block model

[n,q] = size(X);
output = zeros(n,q);
for i = 1:n
    s = sum(exp(X(i,:)));
    for j = 1:q
        output(i,j) = ell*((q*exp(X(i,j))/s) - 1);
    end
end

end

