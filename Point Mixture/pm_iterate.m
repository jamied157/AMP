function [ y ] = pm_iterate( A, y , y_old, gamma, gamma_old , lambda,e )
%AMP_ITERATE
%the conditional expectation of the sum of the delta
%mixture and gaussian noise as in the amp matrix factorisation paper (eqn 2.11)


B = lambda*mean(pm_dnl(y,gamma,e));
c_term = lambda*pm_nl(y_old,gamma_old,e);
p_term = lambda*pm_nl(y,gamma,e);

y = A*p_term - B*c_term;
end

