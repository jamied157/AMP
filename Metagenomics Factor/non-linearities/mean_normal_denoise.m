function [output] = mean_normal_denoise(Z,W,sigma,mu)
%IEEE_UDENOISE is the normal denoising term listed in the ieee paper
[~,p] = size(W);

output =  (mu.*Z + W.*sigma*sqrt(p))./(Z + p*sigma);
end

