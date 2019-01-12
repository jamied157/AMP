function [A, x,overlap ] = gbm_plt( q , n , ell )
%GBM_PLT a function that plots the success of a factorisation of a gaussian
%block model with the amp algorithm, q is the rank, n is the size and ell
%is the lambda signal ratio. It plots how the overlap changes with number
%of iterations.

%{
THIS IS PREVIOUS CODE TO PLOT FOR MULTIPLE VALUES OF LAMBDA/ELL
overlap = zeros(1,length(ell));
for i = 1:length(ell)
    [A ,~ ,x] = ex_gbm(q,n,ell(i));
    X = gbm_factor(A,q,ell(i));
    overlap(i) = gbm_overlap(X,x);
end
plot(ell,overlap)
%}
[A,~,x] = ex_gbm(q,n,ell);

[V,~] = eigs(A,q-1);
X_old = [sqrt(n)*V,zeros(n,1)];X_older = [sqrt(n)*V,zeros(n,1)];
overlap = zeros(21,1);
overlap(1) = gbm_overlap(X_old,x);

for i = 1:20
    X = A*gbm_nl(X_old,ell) - gbm_nl(X_older,ell)*mean(gbm_dnl(X_old,ell),3)';
    overlap(i+1) = gbm_overlap(X,x);
    X_older = X_old;
    X_old = X;
end

plot(0:1:20,overlap,'ko','MarkerFaceColor','r')
title(['q = ',num2str(q),' n = ',num2str(n),' \lambda = ' , num2str(ell)])
xlabel('Iteration')
ylabel('Overlap')
ylim([0,1])
