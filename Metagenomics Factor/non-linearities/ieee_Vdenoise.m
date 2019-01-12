function [T] = ieee_Vdenoise(sigma,T)
%IEEE_VDENOISE is the denoises for the V matrix taken from the ieee paper
[p,n] = size(T);n=n/4;
T0 = exp(-(-2*T+1)./(2*sigma));
%T1 = exp(-((T - 1).^2)./(2*sigma));
for i = 0:n-1
    for j = 1:p
        t0 = T0(j,4*i+1:4*i+4);
        %t1 = T1(j,4*i+1:4*i+4);
        t = (1 + t0.*(sum(t0) - t0)).^(-1);
        T(j,4*i+1:4*i+4) = t;
    end
end
end

