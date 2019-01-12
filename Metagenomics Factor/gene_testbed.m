function [ A,U,V,u,v ] = gene_testbed( m , p , n , noise , ell)
% GENE_TESTBED: a function to test the factorisation algorithm gene_factor
% plots the error in U and V as a function of the iteration number
[A,U,V] = ex_gene_matrix(m,p,n,noise,ell);

[u,v] = gene_factor(A,1/4,ell);

U_loss = squeeze(sum(sum((sqrt(2*m/(p*(p+1)))*gene_Usoftmax(u,1)-U).^2,1),2))/(m*p);
V_loss = squeeze(sum(sum((sqrt(n)*gene_Vsoftmax(v,1)-V).^2,1),2))/(p*4*n);

figure(1);
plot(0:20,U_loss)
title(['MSE of U estimate - s = ',num2str(ell)])

hold on
figure(2);
plot(0:20,V_loss)
title(['MSE of V estimate - s = ',num2str(ell)])

hold on
end

