function [ f ] = pm_nl( y , g , e)
%PM_NL is the non linearity for the 2 point mixure 
%as given in section 2.4
a_plus = sqrt((1-e)/e);
a_minus = 1/a_plus;

f = (-a_minus*exp(- a_minus*y - (a_minus^2*g)/2)*(1 - e) ...
    + a_plus*e*exp(a_plus*y - (a_plus^2*g)/2))./(e*exp(a_plus*y... 
    - (a_plus^2*g)/2) + exp(- a_minus*y - (a_minus^2*g)/2)*(1 - e));
%For small values of e, f can sometimes be NaN however it is usually
%obvious what the value of f should be in this case, the code below fixes
%this
j = 1:1:length(y);
j = j(isnan(f));

for i = j
    if y(i) < 0
        f(i) = -a_minus;
    else
        f(i) = a_plus;
    end
end

end

