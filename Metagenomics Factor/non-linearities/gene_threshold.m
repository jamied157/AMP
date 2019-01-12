function [ V ] = gene_threshold( V )
%GENE_THRESHOLD Applies the hard threshold non linearity to V
[nu,r] = size(V);
I = eye(4);
for i = 1:r
    for j = 0:(nu/4-1)
        V(4*j + 1:4*(j+1),i) = I(:,argmax(V(4*j + 1:4*(j+1),i)));

    end
end
end

