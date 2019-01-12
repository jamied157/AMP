function [ Y , U , V ] = ex_gene_matrix( m , p , n , noise , ell)
%EX_GENE_MATRIX generates an example matrix of dimensions nxp with rank r
%(i.e the number of taxa) if noise = 'g' the matrix is perturbed by random
%gaussian matrix variance ell, if noise = 'm' we sample the matrix from a multinomial
%distribution with number of samples ell

%We first generate the sample/taxa matrix
U = randg(1,m,p);
U = U./sum(U,2);
%then the taxa/site matrix
V = randi(4,p,n);
I = eye(4);
V = I(V,:);
V = reshape(V',[4*n,p]);
%then multiply out
Y = U*V';
%adding noise
if noise == 'g'
    Y = Y + sqrt(ell)*randn(m,4*n);
elseif noise == 'm'
    A = mnrnd(ell,reshape(Y',[4,n*m])');
    Y = reshape(A',[4*n,m])'/ell;
    
    
end

