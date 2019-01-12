function [ A,U,V,u,v ] = test_generic( m , p , n , noise , s , g)
% GENE_TESTBED: a function to test the factorisation algorithm gene_factor
% plots the error in U and V as a function of the iteration number

[A,U,V] = ex_gene_matrix(m,p,n,noise,s);
A = sqrt(s/g)*(1/m)*A;
u = zeros(m,p,21);
v = zeros(4*n,p,21);
a = randi(4,p,n);
I = eye(4);
a = I(a,:);
v(:,:,1) = reshape(a',[4*n,p]);
%{
[~,~,v_0] = svd(A,'econ');
for i = 1:p
    if sum(v_0(:,i)) < 0
        v_0(:,i) = -v_0(:,i);
    end
end
v(:,:,1) = v_0(:,1:p);
%}
p = gene_Vsoftmax(v(:,:,1),5);
u(:,:,1) = 20*A*p*(p'*p)^(-1);
overlap = zeros(101,2);
for i = 2:101
    q = gene_Usoftmax(u(:,:,i-1),2);
    v(:,:,i) = sqrt(20)*(A'*q*(q'*q)^(-1) - gene_Vsoftmax(v(:,:,i-1),5)*gene_Udoftmax(u(:,:,i-1),2)');
    p = gene_Vsoftmax(v(:,:,i),5);
    u(:,:,i) = A*p*(p'*p)^(-1);
end


U_loss = squeeze(sum(sum((sqrt(2*m/(p*(p+1)))*gene_Usoftmax(u,10)-U).^2,1),2))/(m*p);
V_loss = squeeze(sum(sum(-V.*(log(sqrt(n)*gene_Vsoftmax(v,10))) - (1-V).*(log(sqrt(n)*gene_Vsoftmax(v,10))),1),2))/(p*4*n);

figure(1);
plot(0:100,U_loss,'DisplayName','U')
title('U loss')
figure(2);
plot(0:100,V_loss,'DisplayName','V')
title('V loss')
%{
hold on
figure(2);
plot(0:100,V_loss)
title('MSE of V estimate - s = ')
%}
hold on
end

