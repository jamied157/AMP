function [ V2overlap ] = test_initialisation( m , p , n )
% TEST_INTIALISATION: tests how we intitialise the algorithm between three
% choices - svd, E(priors), svd + E(priors), random draw from priors

mse = zeros(4,1);
overlap = zeros(4,1);
V2overlap = 0;
for i = 1:50
    s = 1;
    [A,~,V] = ex_gene_matrix(m,p,n,'m',s);
    %E(priors)
    v = zeros(4*n,p);
    mse(1) =  mse(1)*(i-1)/i +  sum(sum((sqrt(n)*gene_Vsoftmax(v,1)-V).^2,1),2)/(p*4*n)/i;
    overlap(1) = overlap(1)*(i-1)/i + gene_Voverlap(v,V)/i;
    %svd intialisation
    [~,~,v_0] = svd(A,'econ');
    for j = 1:p
        if sum(v_0(:,j)) < 0
            v_0(:,j) = -v_0(:,j);
        end
    end
    v = v_0(:,1:p);
    mse(2) =  mse(2)*(i-1)/i + sum(sum((sqrt(n)*gene_Vsoftmax(v,1)-V).^2,1),2)/(p*4*n)/i;
    overlap(2) = overlap(2)*(i-1)/i + gene_Voverlap(v,V)/i;
    V2overlap = V2overlap*(i-1)/i + gene_Voverlap(v(:,2),V(:,2))/i;
    %SVD + E(priors)
    v = [v(:,1),zeros(4*n,p-1)];
    mse(3) = mse(3)*(i-1)/i + sum(sum((sqrt(n)*gene_Vsoftmax(v,1)-V).^2,1),2)/(p*4*n)/i;
    overlap(3) = overlap(3)*(i-1)/i + gene_Voverlap(v,V)/i;
    %draw from priors
    v = randi(4,p,n);
    I = eye(4);
    v = I(v,:);
    v = reshape(v',[4*n,p]);
    mse(4) = mse(4)*(i-1)/i + sum(sum((sqrt(n)*gene_Vsoftmax(v,1)-V).^2,1),2)/(p*4*n)/i;
    overlap(4) = overlap(4)*(i-1)/i + gene_Voverlap(v,V)/i;
    disp(i)
end
%Plot the results

bar(mse)
title('Mean Squared Error')

figure(2);
bar(overlap)
title('Mean Overlap')
end

