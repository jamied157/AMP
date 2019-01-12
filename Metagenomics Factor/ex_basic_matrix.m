function [ A , U , V] = ex_basic_matrix( m , p , n , g )
%Generates a more basic matrix one can try to use amp to solve

%We first generate the sample/taxa matrix
U = randn(m,p);
%then the taxa/sites matrix
V = randn(n,p);
%multiply out wih noise
A = U*V' + sqrt(g)*randn(m,n);

end

