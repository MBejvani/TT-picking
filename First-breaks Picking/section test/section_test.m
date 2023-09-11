clc, clear,close
%%
load real_shot
[xd,yd] = size(y);
normal_y = y * sparse(diag(1./max(y)));
yn =   normal_y    +   .05*randn(size(normal_y));
%%
tic
l1 = fix(xd/20);%========================
for i = 1:yd
    
    s = yn(:,i);
%     s = normal_y(:,i);
    extra = l1;
    ex_s = [.5*ones(extra,1);s];
    xd1 = length(ex_s);
    %% ARRIVAL TIME VICINITY EXTRACTION
    w1 = gausswin(l1,5);
    c = conv(w1,abs(ex_s));
    d = c(2*l1:end)./(c(l1+1:xd1) + eps);
    [~, m1] = max(d);
    M1(i) = m1;
    l2 = fix(xd/10);%=======================
    a1 = m1-2*l2;if a1<1;a1=1;end
    b1 = m1+2*l2;if b1>xd;b1=xd;end
    x = s(a1:b1);
    n1 = length(x);
    w2 = tukeywin(n1,.4);
    x = x .* w2;
%% TF PRODUCTION
X = repmat(x,1,n1);
W = winmtx( n1,8 );
STFT = fft(W.*X);
TFR = abs(STFT(fix(n1/2)+1:end,:));
%% CHARACTERISTIC FUNCTION DETERMINATION
CFy = sum(TFR,2);
freq = CFy  > 3/4 * max(CFy );freq(1)=1;
CFx = sum(TFR(freq,:));
Smooth_CFx = EPS(CFx,20); 
Diff_smooth_CFx =  diff(Smooth_CFx);
tr = 2*std(Diff_smooth_CFx);
p = find(Diff_smooth_CFx>tr);
m2 = p(1);
[~,m2] = max(Diff_smooth_CFx);
m = CFx(m2+1)-CFx(m2);
zc = round( (m*m2-CFx(m2)) / m );
fb(i) = round(a1+zc+1); 
end

plotseismogram(gain_y);hold;
line([1:yd],fb,'Marker','s','MarkerFaceColor',[1 1 1],...zz
    'MarkerSize',3,'LineWidth',1.5,'Color',[1 .1 .1])
line([1:yd],M1,'Marker','s','MarkerFaceColor',[1 .1 1],...zz
    'MarkerSize',3,'LineWidth',1.5,'Color',[.1 .1 1])

toc
	


