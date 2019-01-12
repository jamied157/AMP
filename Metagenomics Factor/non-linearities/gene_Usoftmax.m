function [ U ] = gene_Usoftmax( U , b)
%GENE_SOFTMAX applies the softmax noninearity to U
[m,r] = size(U);

for i = 1:m
    U(i,:) = exp(b*U(i,:))/sum(exp(b*U(i,:)));
end
U = U*sqrt(r*(r+1)/(2*m));
end

