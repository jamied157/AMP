function [output] = var_normal_denoise(Z,W,sigma)
%IEEE_UDENOISE is the normal denoising term listed in the ieee paper
[~,p] = size(W);

output =  (sigma*Z)./(p*sigma + Z);
end
