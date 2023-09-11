function [ G ] = adapt_winopt( s )
%ADAPT_WINOPT Summary of this function goes here
%   Detailed explanation goes here
N = length(s);
l = fix(N/10);% N/u : u depend on signal
s = padarray(s,[0,l]);
L = 2*l;
H = hankel(s);H = H(1:L,1:N);
alpha = [1:1:100];% :2:
for i=1:length(alpha)
w = gausswin(L,alpha(i));w = w/sum(w);
Tt = sparse(diag(w)) * H;
ZT = fft(Tt);ZT = ZT(1:fix(L/2),:);
CF(i,:) = sum(log10(abs(ZT) + eps));
end
[aa,P] = min(CF);
for i=1:N
w = gausswin(L,alpha(P(i)));w = w/sum(w);
W(:,i) = w;
end
G = fft(W .* H,N);
end

