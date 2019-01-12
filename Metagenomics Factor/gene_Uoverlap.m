function [output] = gene_Uoverlap(u,U,b)

u = gene_Usoftmax(u,b);
output = sum(sum(U.*u))/(norm(u,'fro')*norm(U,'fro'));
end

