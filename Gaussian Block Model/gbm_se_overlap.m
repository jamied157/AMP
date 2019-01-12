function [ output ] = gbm_se_overlap( M, Q , q)
%GBM_SE_OVERLAP given the state evolution iterates M and Q this calculates
%the estimated overlap between the amp estimate. 
I = eye(q);
overlap = zeros(1,1000);
P = eye(q) - 1/q*ones(q);
for j = 1:1000
    G = randn(q,1);
    a = q*M + sqrtm(Q)*G;
    b = P;
    overlap(j) = sum(sum(a*.b));
end
output = mean(mean(overlap));
end

