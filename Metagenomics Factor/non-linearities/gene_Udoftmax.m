function [ J ] = gene_Udoftmax( U , b)
%GENE_DOFTMAX is the mean jacobian of the softmax function applied to U
[m,r] = size(U);
J = zeros(r,r,m);
P = exp(b*U);
%the aim is to loop over the rows of U calculating the jacobian at each one
%and then averaging them at the end
for k = 1:m
    s = sum(P(k,:));
    for i = 1:r
        J(i,i,k) = b*P(k,i)*(s - P(k,i))/(s^2);
        for j = i+1:r
            J(i,j,k) = -(b*P(k,i)*P(k,j))/(s^2);
            J(j,i,k) = J(i,j,k);
        end
    end
end

J = mean(J,3)*sqrt(r*(r+1)/(2*m));
end

