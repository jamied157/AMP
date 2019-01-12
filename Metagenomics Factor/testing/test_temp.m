function [Uoverlap] = test_temp(m,p,n,noise,s,g)
%TEST_TEMP tries different temperature values in the softmax for U and V
%and plots results for the overlap
Uoverlap = zeros(100,1);
%Uoverlap = zeros(100,1);
for i = 1:100
    [A,U,~] = ex_gene_matrix(m,p,n,noise,s);
    A = sqrt(s/g)*(1/m)*A;
    %u = zeros(m,p,21);
    u = sqrt(1/n)*ones(m,p,21);
    [u_0,~,~] = svd(A,'econ');
    for k = 1:p
        if sum(u_0(:,k)) < 0
            u_0(:,k) = -u_0(:,k);
        end
    end
    u(:,:,1) = u_0(:,1:p);
    for j = 1:100
        Uoverlap(j) = ((i-1)/i)*Uoverlap(j) + (1/i)*gene_Uoverlap(u(:,:,1),U,2*j+50);
    end
    disp(i)
end

