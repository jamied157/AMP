function [ u , v ] = gene_orthog_factor( A , g , s)
% GENE FACTOR Our AMP factorisation algorithm for the genomics factor model

[m,n] = size(A); n = n/4; 
A = sqrt(s/g)*(1/m)*A;

[~,lambda,v_0] = svd(A);
r = sum(sum(lambda>2));
for i = 1:r
    if sum(v_0(:,i)) < 0
        v_0(:,i) = -v_0(:,i);
    end
end
%lambda = lambda(1:r,1:r);
u = zeros(m,r,21);
v = zeros(4*n,r,21);

v(:,:,1) = v_0(:,1:r);
r = gene_Vsoftmax(v(:,:,1),1);
u(:,:,1) = A*r*(r'*r)^(-1);
%lambda = new_lambda(A,u(:,:,1),v(:,:,1),s,g);

for i = 2:21
    
    v(:,:,i) = A'*gene_Usoftmax(u(:,:,i-1),1) - gene_Vsoftmax(v(:,:,i-1),1)*gene_Udoftmax(u(:,:,i-1),1)'; 
    
    u(:,:,i) = A*gene_Vsoftmax(v(:,:,i),1);
    
    %lambda = new_lambda(A,u(:,:,i),v(:,:,i),lambda,s,g);
end
end