function [ A , A_0 ,x] = ex_gbm( q , n , ell)
%EX_GBM generates a noisy matrix from the gaussian block model
s = randi(q,1,n);

x = zeros(n,q);
I = eye(q);
P = I - 1/q*ones(q);

for i = 1:n
    x(i,:) = P*I(s(i),:)';
end
A_0 = (q/n)*x*x';
A = ell*A_0 + GOE(n);

end

