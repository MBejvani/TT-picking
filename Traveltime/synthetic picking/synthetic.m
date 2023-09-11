clc, clear all, close all
%%
dt = .004 ;
trng = 2 ; 
df = 1/trng ;
t = 0:dt:trng-dt ;
n = fix(trng/dt);
%% signal #1
r1 = [1:n]*0;r2 = r1;r3 = r2;
r1([100,180,230]) = [1.7 -1.2 -.5];
r2([160,250,310]) = [1.2 +1 .7];
r3([280,385,450,420]) = [-.9 .7 -.5 -.6];
r = r1+r2+r3;
fd1 = 25;fd2 = 15;fd3 = 10;
x1 = (1-2*(pi*fd1*(t-trng/2)).^2) .* exp(-(pi*fd1*(t-trng/2)).^2) ;
x2 = (1-2*(pi*fd2*(t-trng/2)).^2) .* exp(-(pi*fd2*(t-trng/2)).^2) ;
x3 = (1-2*(pi*fd3*(t-trng/2)).^2) .* exp(-(pi*fd3*(t-trng/2)).^2) ;
y1 = conv(r1,x1);y1 = y1(fix(n/2)+1:n+fix(n/2));
y2 = conv(r2,x2);y2 = y2(fix(n/2)+1:n+fix(n/2));
y3 = conv(r3,x3);y3 = y3(fix(n/2)+1:n+fix(n/2));
y = y1+y2+y3;
%% graph result
figure
G = [y;y3;y2;y1]';

tt=-.2:.001:.2;
fm=5;

[nz,nx]=size(G);

trmx= max(abs(G));
 amx=mean(trmx);  
 x=[1:nx]; z=[1:nz];
 scal =.7;

 % take the average as dx
	dx1 = abs(x(2:nx)-x(1:nx-1));
 	dx = median(dx1);

 dz=z(2)-z(1);
 xmx=max(max(G)); xmn=min(min(G)); 

 if scal == 0; scal=1; end;
 G = G * dx /amx; 
 G = G * scal;

%  fprintf(' PlotWig: data range [%f, %f], plotted max %f \n',xmn,xmx,amx);
 
% set display range 
x1=min(x)-2.0*dx; x2=max(x)+2.0*dx;
z1=min(z)-dz; z2=max(z)+dz;
 
set(gca,'NextPlot','add','Box','on', ...
  'XLim', [x1 x2], ...
  'YDir','reverse', ...
  'YLim',[z1 z2]);
 

	fillcolor = [0 0 0];
	linecolor = [0 0 0];
	linewidth = 1;

	z=z'; 	% input as row vector
	zstart=z(1);
	zend  =z(nz);

for i=1:nx,
   
  if trmx(i) ~= 0;    % skip the zero traces
	tr=G(:,i); 	% --- one scale for all section
  	s = sign(tr) ;
  	i1= find( s(1:nz-1) ~= s(2:nz) );	% zero crossing points
	npos = length(i1);


	%12/7/97 
	zadd = i1 + tr(i1) ./ (tr(i1) - tr(i1+1)); %locations with 0 amplitudes
	aadd = zeros(size(zadd));

	[zpos,vpos] = find(tr >0);
	[zz,iz] = sort([zpos; zadd]); 	% indices of zero point plus positives
	aa = [tr(zpos); aadd];
	aa = aa(iz);

	% be careful at the ends
		if tr(1)>0, 	a0=0; z0=1.00;
		else, 		a0=0; z0=zadd(1);
		end;
		if tr(nz)>0, 	a1=0; z1=nz; 
		else, 		a1=0; z1=max(zadd);
		end;
			
	zz = [z0; zz; z1; z0];
 	aa = [a0; aa; a1; a0];
		

	zzz = zstart + zz*dz -dz;

	patch( aa+x(i)+1 , zzz,  fillcolor);

	line( 'Color',[1 1 1],'EraseMode','background',  ...
         'Xdata', x(i)+[0 0]+1, 'Ydata',[zstart zend]); % remove zero line

%'LineWidth',linewidth, ...
%12/7/97 'Xdata', x(i)+[0 0], 'Ydata',[z0 z1]*dz);	% remove zero line

	line( 'Color',linecolor,'EraseMode','background',  ...
	 'LineWidth',linewidth, ...
	 'Xdata', tr+x(i)+1, 'Ydata',z);	% negatives line

   else % zeros trace
	line( 'Color',linecolor,'EraseMode','background',  ...
	 'LineWidth',linewidth, ...
         'Xdata', [x(i) x(i)]+1, 'Ydata',[zstart zend]);
   end;
end

line(r*.7,t/dt)
% plot(t,y1+8.5),hold on
% plot(t,y2+6.5)
% plot(t,y3+4.5)
% plot(t,y+2.5,'k')
% stem(t,r,'Marker','none')
% % hold off
% set(gca,'YTick',[1:500])
% set(gca,'YTickLabel',{num2str(t)})
% 'reflectivity';'synthetic';'0.5 Hz';'1.5 Hz';'2.5 Hz'
% xlabel Time/sec
% text([.1,.1],[5.5,7.5],'+','FontSize',16)
% text([.1],[3.5],'=','FontSize',16)
%--------------------------------------------------------------------------
% clear -regexp ^x ^d ^a z t r f i l n;
clearvars -except y dt trng df t n
save y