clc, clear all, close all
%%
load y
%% TF
S = repmat(y',1,n);
f = linspace(0,1/dt/2,n/2);
W = winmtx_1( n,10 );
TF = fft(W.*S);
%% BandWidth
[ BW,ATFR ] = BandWidth( TF,dt );
%% TravelTime
dTF = fft(sqrt(-1) * sparse(diag([0:n-1]*dt)) * real(ifft(TF)));
dTF = dTF(fix(n/2)+1:end,:);
TF = TF(fix(n/2)+1:end,:);
Tau = (real(TF).*imag(dTF)-real(dTF).*imag(TF)) ./ (conj(TF).*TF+10e-7);% eps
T = repmat(t,n/2,1);
% [r,c,v] = find(diff(sign(T-Tau),1,2)>0);
Tau(diff(sign(T-Tau),1,2)>0) = 0;
Tau(diff(sign(T-Tau)+1,1,2)>0) = 0;
zpl = double(diff(sign(T-Tau)+1,1,2)>0);
%%
BW = fix(((1/dt/2)-BW)./df);
BW(BW(:,:)>fix(n/2)) = fix(n/2);
BW(BW(:,:)<=0) = 1; 
bw = fix(mean(BW,2));
%--------------------------------------------------------------------------
tau = mean(Tau(min(bw):max(bw),:));
%%
dif = t-tau;
[out, index] = zcr(dif, 'p', 0);
d = [];
d(index) = 1;d(1,[1 n])=0;
%% graph result
figure;
imagesc(ATFR);hold;plot(BW(2,:),'w:');plot(BW(1,:),'w:');
figure
imagesc(Tau);hold;plot(BW(2,:),'k-');plot(BW(1,:),'k-');
figure
stem(t,d,'Marker','none');ylim ([-.15 1.2])