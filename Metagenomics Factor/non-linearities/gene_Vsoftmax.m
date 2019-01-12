function [ V ] = gene_Vsoftmax( V , b )
%GENE_4SOFTMAX applies the softmax noninearity to V
[n,r] = size(V); n = n/4;

for i = 0:n-1
    for j = 1:r
        v = (1/b)*V(4*i+1:4*i+4,j);
        v = v - max(v);
        v = exp(v);
        v = v/sum(v);
        V(4*i+1:4*i+4,j) = v;
    end
end
V = sqrt(1/n)*V;
end