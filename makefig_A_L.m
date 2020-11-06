%FILE: hyp_diam_nonconvex_v1.m
% 2020-10-07
clear
addpath bie fmm files
% 
% G = The doubly connected domain between the outer polygon and the inner
%     polygon
% E = The simply connected domain inside the inner polygon
% 
% This code to compute:
% cap   = the capacity of the domain G
% hd    = the hyp diameter of the set E with respect to the outer polygon
%
%
%%
% choose values of a, b, c (as in the figure: 2020-10-07)
a      =   0.2;
ver    = [ 0 ; 3 ; 3+i ; 2+i ; 2+a*i ; 1+a*i ; 1+i ; i];
alpha  =   0.5+0.5i; % alpha must be in the domain E (the domain between the
                  % two ploygons.
%%
n        =  2^13; 
t        = (0:2*pi/n:2*pi-2*pi/n).';
[et,etp] =  polygonp(ver,n/8);
%%   
[x,y] =  meshgrid(linspace(0,3,3000),linspace(0,1,1000));
z     =  x+i*y;
zv    =  z(:);
[in,on] = inpolygon(real(zv),imag(zv),real(ver),imag(ver));
zv(~in) = NaN+i*NaN; 
zv(on)  = NaN+i*NaN; 
zvb     = zv(abs(zv)>=0).';
%%
rzvb   = hypdist (et,etp,n,alpha,zvb,alpha);
rzv    = NaN(size(zv));
rzv(abs(zv)>=0) = rzvb;
rz     = NaN(size(z));
rz(:)  = rzv;
%%
hrad  = [0.5,1.5,4,10,16,21,23,23.5,23.75];
figure
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on;     
contour(real(z),imag(z),rz,hrad,'color','b','LineWidth',1.5)
plot(real(et),imag(et),'k','LineWidth',1.5);
% plot(real(ver),imag(ver),'ok');
plot(real(alpha),imag(alpha),'ok','MarkerSize',4,'MarkerFaceColor','k');
% plot(real(zvb),imag(zvb),'or');
% text(4.35,2.8,{'$\alpha$'},'Interpreter','LaTeX','FontSize',20)
axis equal
grid(gca,'minor')
grid on
set(gca, 'XMinorTick','on')
set(gca, 'YMinorTick','on')
ax=gca;
ax.GridAlpha=0.5;
ax.MinorGridAlpha=0.5;
box on
set(gca,'FontSize',22)
set(gca,'LooseInset',get(gca,'TightInset'))
axis([-0.2 3.2 -0.2 1.2])
print -dpdf HypCir_2.pdf
%%