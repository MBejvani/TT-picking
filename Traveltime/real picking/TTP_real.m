clc, clear all, close all
%%
% load real
load partial_real
[xd,yd] = size(y);
% for o=1:40
%     if o<=60
%     y(1:288-4.2*o,o)=0;
%     else
%     y(1:4.5*o-260,o)=0;
%     end
% end
normal_y = y * sparse(diag(1./max(y)));
% [y_agc] = gain(y,.004,'agc',.5,2);
for i = 1:yd
s = y(:,i)';
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
Tau = (real(TF).*imag(dTF)-real(dTF).*imag(TF)) ./ (conj(TF).*TF+eps);% eps % ...
T = repmat(t,n/2,1);
% [r,c,v] = find(diff(sign(T-Tau),1,2)>0);
Tau(diff(sign(T-Tau),1,2)>0) = 0;
% Tau(diff(sign(T-Tau)+1,1,2)>0) = 0;
zpl = double(diff(sign(T-Tau)+1,1,2)>0);
% zpl = zpl(1:fix(n/4),:);
%--------------------------------------------------------------------------
l = 10;
% p0 = conv2(ones(1,l),zpl);p0 = p0(:,fix(l/2):n-1+fix(l/2)-1);
p1 = zpl(end-3:end,:);
% sp0 = sum(p0);
sp1 = sum(p1);
% ps = double((sp1>2).*sp0)>0;
ps = double(sp1>0);
PS(:,i) = ps';
end
%% graph result
figure(1)
imagesc([1:yd],t,y_agc);
rectangle('Position',[15,.75,20,1.5],'EdgeColor','b','LineWidth',1)
rectangle('Position',[15,2.75,20,1],'EdgeColor','b','LineWidth',1)
% contourf(PS,10);
% set(gca,'YTick',[]);
ylabel TWT(s)
xlabel 'Trace #'
% mycmap = get(figure(1),'Colormap');
% save('MyColormaps','mycmap')
load('MyColormaps','mycmap')
set(figure(1),'Colormap',mycmap)
% mycmap = get(figure(1),'Colormap');
% save('MyColormaps','mycmap')
hold
plotseismogram (PS,.9,[1:yd],t)



