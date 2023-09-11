function [ BW,ATFR ] = BandWidth( TF,dt )
%BANDWITHD Summary of this function goes here
%   Detailed explanation goes here
n = max(size(TF));
df = 1./(dt*n);
f = df*[fix(n/2)-1:-1:0];
ATFR = abs(TF(fix(n/2)+1:end,:)).^1.5;
miu = (f * ATFR) ./ (sum(ATFR) );% + eps
f_miu = repmat(f',1,n)-repmat(miu,fix(n/2),1);
sigma = sqrt(sum(f_miu.^2 .* ATFR) ./ (sum(ATFR) ));% + eps
BW = [miu - sigma ;miu + sigma];
end

