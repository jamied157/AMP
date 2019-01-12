function [ y , gamma ] = pm_factor( A , e , lambda)
%A simple application of approximate message passing applied to
%factorising a simple one rank matrix from the point mixture in section 2
n = size(A,1);
%We first compute the principal eigenvector
[v,~] = eigs(A,1);
if sum(sign(v)) > 0
    v = -v;
end
%next we setup our iteration functions
%lambda = .5*(abs(ell) + sqrt(abs(ell)^2 -4));

y_old = sqrt(n*(lambda^2)*((lambda^2)-1))*v;
gamma_old = (lambda^2) - 1;
y = A*lambda*pm_nl(y_old,gamma_old,e);

for i = 1:200
    gamma = lambda^2*(1-mmse_est(gamma_old,e));
    y_older = y_old;
    y_old = y;
    y = pm_iterate(A,y_old,y_older,gamma,gamma_old,lambda, e);
    gamma_old = gamma;
end

y = pm_nl(y,gamma,e);
end



