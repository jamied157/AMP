function [ J ] = gene_Vdoftmax( V , b)
%gene_Vdoftmax Calculates the jacobian for the V matrix
[n,r] = size(V); n = n/4;
J = zeros(r,r);
for j = 0:r-1
    J(4*j+1:4*(j+1),4*j+1:4*(j+1)) = gene_Udoftmax(reshape(V(:,j+1),4,n)',b);
end

end

