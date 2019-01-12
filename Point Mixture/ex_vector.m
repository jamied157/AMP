function [ x ] = ex_vector( n , e)
%ex_vector Gives a column n-vector of the form given in the 2.4 example of the
%amp factorisation notes
x = zeros(n,1);
for i = 1:n
    r = rand;
    if r < e
        x(i) = sqrt((1-e)/e);
    else
        x(i) = -sqrt(e/(1-e));
    end
end

end

