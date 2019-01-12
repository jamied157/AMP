function [ v , lambda ] = prin_vec( A , k)
%prin_vec - finds k largest eigenvectors of a matrix
[v,lambda] = eigs(A,k);

end

