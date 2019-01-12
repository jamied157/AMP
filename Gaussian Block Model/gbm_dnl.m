function [ output ] = gbm_dnl( X , ell )
%GBM_DNL the differential of the non-linearity function used in the amp
%factorisation of the gaussian block model
[n,q] = size(X);
output = zeros(q,q,n);
for i = 1:n
    s = sum(exp(X(i,:)));
    for j = 1:q
        for k = 1:q
            if j == k
                output(j,k,i) = ell*exp(X(i,j))*q*(s - exp(X(i,j)))/(s^2);
            else
                output(j,k,i) = -ell*q*exp(X(i,j) + X(i,k))/(s^2);
            end
        end
    end
end

end

