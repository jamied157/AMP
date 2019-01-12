function [ output ] = gbm_se( n, q , x , x_0 , ell)
%GBM Provides a state evolution analysis of the factorisation problem in
%section 3 of the paper
output = zeros(1,20);
%initialise the algorithm
M = (x_0)'*x/n;
I = eye(q);
Q = (1/ell)*I;
Q(q,q) = 0;
P = eye(q) - 1/q*ones(q);
%calculate the first state evolution estimate
output(1) = gbm_se_overlap(M,Q,q);

for i = 1:20
    monte_m = zeros(q,q,1000,q);
    monte_q = zeros(q,q,1000,q);
    for j = 1:1000
        G_m = randn(q,1);
        G_q = randn(q,1);
        for k = 1:q
            monte_m(:,:,j,k) = gbm_nl((q*M*I(:,k) + sqrtm(Q)*G_m)',ell)'*I(:,k)'*P;
            monte_q(:,:,j,k) = gbm_nl((q*M*I(:,k) + sqrtm(Q)*G_q)',ell)'*gbm_nl((q*M*I(:,k) + sqrtm(Q)*G_q)',ell);
        end
    end
    M = ell*mean(mean(monte_m,3),4);
    Q = mean(mean(monte_q,3),4);
    output(i+1) = gbm_se_overlap(M,Q,q);
end

