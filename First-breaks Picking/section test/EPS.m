function [ smooted_s ] = EPS( s,l )
%% --------->>> WINMTX Summary of this function goes here
%   Detailed explanation goes here
N = length(s);
s = padarray(s,[0 l-1]);
xidx = (1:l)';
yidx = 0:(N+l-2);
h = xidx(:,ones(N+l-1,1)) + yidx(ones(l,1),:);
H = s(h);
smooted_s(N) = 0;
for i = 1:N
    [~,rr] = min(std(H(:,i:i+l-1)));
    smooted_s(i) = mean(H(:,i-1+rr));
end
end

