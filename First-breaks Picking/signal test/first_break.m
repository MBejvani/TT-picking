close all
clear all ;clc
%% LOADING DATA
load real_shot
y = y(:,30);
n0 = length(y);
dt = .04;
t0 = dt*(0:n0-1);
y = y/max(abs(y)); % unnecessary code
%% ADDING NOISE
yn =   y    +  .2 * randn(size(y));
SNR = 10 * log10(sum(y.^2)./sum((y-yn).^2));
%% ARRIVAL TIME VICINITY EXTRACTION
l = fix(n0/20);
w1 = gausswin(l,5); % selected parameter
c = conv(w1,abs(yn));
d = c(l:end)./(c(1:n0) + eps);
d(1:l) = 0;
[~,m1] = max(d);
b1 = fix(m1-2*l+1);
b2 = fix(m1+2*l);
x = yn(b1:b2); % selected parameter
n1 = length(x);
w2 = tukeywin(n1,.4); % selected parameter
% x = x.*w2;
t1 = dt*(b1:b2);
f1 = linspace(0,1/dt/2,n1/2);
%% TF PRODUCTION
X = repmat(x,1,n1);
W = winmtx_1( n1,8 ); % selected parameter
STFT = fft(W.*X);
TFR = abs(STFT(fix(n1/2)+1:end,:));
%% CHARACTERISTIC FUNCTION DETERMINATION
abs_fft = abs(fft(x));
abs_fft = abs_fft(fix(n1/2)+1:end,:);
CFy = sum(TFR,2);
amount = 4/5 * max(CFy); % selected parameter
freq = CFy > amount;
freq_ = find(diff(freq));
CFx = sum(TFR(freq,:));
CFx = CFx/max(CFx); % unnecessary code
Smooth_CFx = EPS(CFx,20); % selected parameter
Diff_smooth_CFx =  diff(Smooth_CFx);
[xm2,m2] = max(Diff_smooth_CFx);
m = CFx(m2+1)-CFx(m2);
y_line = m*(1:n1) - m*m2 + CFx(m2);
zc = round( (m*m2-CFx(m2)) / m ); % y=0 , x:
fb = round(m1-n1/2+zc+1);
user_pick = 209;
%% figures
figure(5)
plot(t0,yn,'k');hold;axis tight;box on;xlabel Time(sec);ylabel Amplitude
plot(dt*[b1 b1],[min(yn) max(yn)],'r-.',...
    dt*[b2 b2],[min(yn) max(yn)],'r-.');
figure(6)
plot(t0,d,'k',dt*m1,d(m1),'sb');hold;ylim([0 1.1*max(d)]);box on;xlabel Time(sec);ylabel Amplitude
plot(dt*[b1 b1],[0 1.2*max(d)],'r-.',...
    dt*[b2 b2],[0 1.2*max(d)],'r-.');
figure(1)
plot(CFy,'. b-');hold;plot(abs_fft,'k:');plot([freq_(1) freq_(1)],[0 1.2*max(CFy)],'k',...
    [freq_(2) freq_(2)],[0 1.2*max(CFy)],'k');
plot([0 fix(n1/2)],[amount amount],'-');
axis tight;set(gca,'XTick',[]);ylabel Amplitude
figure(2)
imagesc(t1,f1,TFR);ylabel Frequency(Hz)
set(gca,'XTick',[])
set(gca,'YTick',f1(1:14:end))
set(gca,'YTickLabel',{num2str(fliplr(round(f1(1:14:end)))')})
hold;plot(dt*[b1 b2],f1([freq_(1) freq_(1)]),'w-',...
    dt*[b1 b2],f1([freq_(2) freq_(2)]),'w-');
figure(7)
plot(t1,CFx,'-b.');axis tight;xlabel Time(sec)

figure(3)
plot(t1,CFx,':k');hold;plot(t1,Smooth_CFx,'-b');plot(t1,y_line,'-r');...
    plot(dt*([1,n1]+b1-1),[0,0],'k',dt*(m2+b1-1),CFx(m2),'sk',dt*(zc+b1-1),0,'ok','LineWidth',1,...
    'MarkerSize',6);ylim([-.1 1.1]);xlim([t1(1),t1(end-1)]);;xlabel Time(sec);
    legend('ch fun',...
    'smoothed ch fun','tangent line','zero level');
figure(8)
h = stem(t1(1:end-1),Diff_smooth_CFx,'k');set(h,'Marker','none'),...
    set(get(h,'BaseLine'),'LineStyle','-');hold,h1 = stem(dt*(m2+b1-1),xm2,'sk','fill');...
    set(h1,'MarkerFaceColor','red','MarkerSize',3);ylim([min(Diff_smooth_CFx)-.1,...
    max(Diff_smooth_CFx)+.1]);xlim([t1(1),t1(end-1)]);xlabel Time(sec);ylabel 'dd derivative'

figure(4)
subplot 211;plot(y);hold;plot([fb fb],[-1 1],'r-.',fb,y(fb),...
    'or');ylim([-1,1]);rectangle('Position',[fb-l,-.25,2*l,.5]);
subplot 212;plot(yn);hold;plot([fb fb],[-1 1],'r-.',fb,y(fb),...
    'or');ylim([-1,1]);rectangle('Position',[fb-l,-.25,2*l,.5]);

