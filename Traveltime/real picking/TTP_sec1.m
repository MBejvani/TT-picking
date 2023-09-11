clc, clear all, close all
%%
load sec1;
y = y(300:550-1,50:250);
[xd,yd] = size(y);
normal_y = y * sparse(diag(1./max(y)));
[y_agc] = gain(y,.004,'agc',.5,2);
for i = 1:yd
s = y(:,i)';
end
n = length(s);
t = 0:dt:(n*dt)-dt;
df = 1/max(t) ;
%% TF
S = repmat(s',1,n);
f = linspace(0,1/dt/2,n/2);
W = winmtx_1( n,5 );
TF = fft(W.*S);
%% BandWidth
ATFR = abs(TF(fix(n/2)+1:end,:));
%% TravelTime
dTF = fft(sqrt(-1) * sparse(diag([0:n-1]*dt)) * real(ifft(TF)));
dTF = dTF(fix(n/2)+1:end,:);
TF = TF(fix(n/2)+1:end,:);
Tau = (real(TF).*imag(dTF)-real(dTF).*imag(TF)) ./ (conj(TF).*TF+1e-10);% eps % ...
T = repmat(t,n/2,1);
% [r,c,v] = find(diff(sign(T-Tau),1,2)>0);
Tau(diff(sign(T-Tau),1,2)>0) = 0;
zpl = double(diff(sign(T-Tau)+1,1,2)>0);
%--------------------------------------------------------------------------
l1 = 10;
l2 = 40;
bn = fix((xd-1)/l1);
for j = 1:fix((xd-1)/l1)
    sum(zpl(j*l1:l1));
% p1 = [zpl(1:l,:);zpl(end-l:end,:)];
p1 = [zpl(end-l1:end,:)];
sp1 = sum(p1);
% sp1 = sp1-mean(sp1);
% ps = double((sp1>2).*sp0)>0;
ps = double(sp1>0);
% [lmaxima,indices] = localmax(sp1,[],false);
% ps = double(lmaxima>0);
PS(:,i) = ps';
end
%% graph result
figure(1)
imagesc(y_agc);
% mycmap = get(figure(1),'Colormap');
% save('MyColormaps','mycmap')
load('MyColormaps','mycmap')
set(figure(1),'Colormap',mycmap)
% mycmap = get(figure(1),'Colormap');
% save('MyColormaps','mycmap')
hold
plotseismogram (PS,1.2)
% % rectangle('Position',[5,350,45,350],'EdgeColor','k','LineWidth',5)
% % imagesc(ATFR);
% % imagesc(Tau);
% % imagesc(zpl);colormap gray



