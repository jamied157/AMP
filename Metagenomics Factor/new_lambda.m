function [Lambda] = new_lambda( A , U , V ,s,g)
%NEW LAMBDA finds the estiamted noise matrix (lambda) in the meta genomics
%model
[m,nu] = size(A);
[~,p] = size(U);
A = reshape(A,[m*nu,1]);
M = zeros(p,m*nu);
U = gene_Usoftmax(U,1,s,nu/4,g);
V = gene_Vsoftmax(V,1,nu/4);
for i = 1:p
    M(i,:) = reshape(U(:,i)*V(:,i)',[1,m*nu]);
end

Lambda = diag(M*M'\(M*A));
end
