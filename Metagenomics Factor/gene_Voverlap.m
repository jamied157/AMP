function [output] = gene_Voverlap(v,V)
v = gene_Vsoftmax(v,b);
output = sum(sum(V.*v))/(norm(v,'fro')*norm(V,'fro'));
end

