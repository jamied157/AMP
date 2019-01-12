function [ output ] = gbm_overlap( X , x )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[~,q] = size(x);

P = perms(1:q);
I = eye(q);
output = zeros(1,factorial(q));

for i = 1:factorial(q)
    a = I(P(i,:),:);
    output(i) = sum(sum(X.*(x*a)))/(norm(X,'fro')*norm(x,'fro'));
end
output = max(output);
