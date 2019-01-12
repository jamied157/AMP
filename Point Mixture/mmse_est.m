function [ mmse ] = mmse_est( g , e )
%MMSE_EST gives an estimate of the minimum mean square error as given in
%equation (2.7) of the amp factorisation paper with the mixture 
%delta function as the distritbution (example 2.4?)
a_plus = sqrt((1-e)/e);
a_minus = 1/a_plus;
mmse = zeros(1,length(g));
G = randn(10000,length(g));
for i = 1:length(g)
    mmse(i) = mean(e*(a_plus-pm_nl(g(i)*a_plus+sqrt(g(i))*G(:,i),g(i),e)).^2 ...
        + (1-e)*(-a_minus-pm_nl(-g(i)*a_minus+sqrt(g(i))*G(:,i),g(i),e)).^2);
end

end

