function [ fb ] = FistBreak( s )
%FISTBREAK Summary of this function goes here
%   Detailed explanation goes here
%% PRIMARY PROCESS
s = s(:)';
n0 = length(s);
s = s(1:2*fix(n0/2));
l = fix(n0/20);
s = [zeros(1,l) s];
n0 = length(s);
%% ARRIVAL TIME VICINITY EXTRACTION
w1 = gausswin(l,5);% selected parameter
c = conv(w1,abs(s));
d = c(l:end)./(c(1:n0) + eps);
d(1:l) = 0;
[~, m1] = max(d);
x = s(m1-l+1:m1+l);% selected parameter
n1 = length(x);
w2 = tukeywin(n1,.4);w2 = w2(:)';% selected parameter
x = x .* w2;
%% TF PRODUCTION
X = repmat(x',1,n1);
W = winmtx( n1,8 );% selected parameter
STFT = fft(W.*X);
TFR = abs(STFT(n1/2+1:end,:));
%% CHARACTERISTIC FUNCTION DETERMINATION
CFy = sum(TFR,2);
freq = CFy  > 4/5 * max(CFy );% selected parameter
CFx = sum(TFR(freq,:));
Smooth_CFx = EPS(CFx,20); % selected parameter
Diff_smooth_CFx =  diff(Smooth_CFx);
[~,m2] = max(Diff_smooth_CFx);
m = CFx(m2+1)-CFx(m2);
zc = round( (m*m2-CFx(m2)) / m );
fb = round(m1-n1+zc+1); % if n1 = 2*l
end

%__________________________________________________________________________
function [ W ] = winmtx( N,sd )
%% -------->>> STFT Window Matrix Function
%   Detailed explanation goes here
W = zeros(N);
if length(sd) == 1
    W = toeplitz(fftshift(1/(sqrt(2*pi)*sd) * ...
        exp(-1/2*((-fix(N/2):fix(N/2)-1)/sd).^2)));
else
    for i = 1:N
        W(:,i) = circshift(1/(sqrt(2*pi)*sd(i)) * ...
            exp(-1/2*((-fix(N/2):fix(N/2)-1)/sd(i)).^2)',fix(N/2)+i);
    end
end
end

%__________________________________________________________________________
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



