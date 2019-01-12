function [ o ] = pm_overlap( y,x )
%OVERLAP finds how close two vectors are to one another
o = dot(y,x)/(norm(x)*norm(y));


end

