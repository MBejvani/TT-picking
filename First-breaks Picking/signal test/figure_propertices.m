axis([0 1 0 max(f)])
xlabel 'time/sec';ylabel 'freq/Hz'
print(figure(1),'-djpeg','-r500','picture1')

figure(1),imagesc(abs(TFR0(fix(N/2)+1:end,:)));colormap(jet)
set(gca,'XTick',[1,512,1024])
set(gca,'XTickLabel',{'0';'0.5';num2str(t(1024))})
set(gca,'YTick',[1,256,512])
xlabel 'time/sec';ylabel 'freq/Hz'
set(gca,'YTickLabel',{'512';'256';'0'})
print(figure(1),'-djpeg','-r500','picture2')

scrsz = get(0,'ScreenSize');
figure('Position',[scrsz(3)/4 scrsz(4)/4 scrsz(3)/2 scrsz(4)/2])
get(gcf,'position')
set(gcf,'position',[200   104   600   380]);
h22 = plot(scf);
zdir = [0 0 1];
rotate(h22,zdir,-90)
set(gca,'Units')
grid
set(gca,'XDir','reverse')
set(gca,'XGrid','on','YGrid','off')
grid(gca,'minor')
box on
% subplot('position',[.1,.4,.8,.5])