function [ A ] = GOE( n )
%GOE generates a Gaussian orthogonal ensemble of size nxn
A=zeros(n,n);
for i=1:n
    A(i,i) = sqrt(2/n)*randn;
    for j=i+1:n
        A(i,j)= sqrt(1/n)*randn;
        A(j,i) = A(i,j);
    end
end
end


