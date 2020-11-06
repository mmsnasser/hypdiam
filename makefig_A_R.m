%FILE: hyp_diam_2cubes.m
% 2020-10-29
clear
addpath bie fmm files
% 
%%
% choose the value of a
a = 0.2;
h = 0.2;
% 
% vertices of the external polygon. The vertices are counterclockwise
% orriented
ver_out = [ 0 ; 3 ; 3+i ; 2+i ; 2+a*i ; 1+a*i ; 1+i ; i];
% vertices of the internal polygon. The vertices are clockwise orriented
ver_in  = 0.5+0.5i+h.*[-1-i ; -1+i ; 1+i ; 1-i];
alpha  =  1.5+0.1i; % alpha must be in the domain E (the domain between
                       % the two ploygons.
z2     =  0.5+0.5i;    % z2 must be inside the inner polygon
%%
n      =  2^13; 
t      = (0:2*pi/n:2*pi-2*pi/n).';
[eto,etop]=polygonp(ver_out,n/8);
[eti,etip]=polygonp(ver_in,n/4);
zet =  eti;
et  = [eto ; eti];
etp = [etop; etip];
%%   
figure
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on;     box on
k=1; crv = et((k-1)*n+1:k*n);crv(end+1) = crv(1); 
fill(real(crv),imag(crv),[0.85 0.85 0.85])
plot(real(crv),imag(crv),'k','LineWidth',1.5);
k=2; crv = et((k-1)*n+1:k*n);crv(end+1) = crv(1); 
fill(real(crv),imag(crv),[1 0.0 0.0])
plot(real(crv),imag(crv),'-b','LineWidth',1.5);
% plot(real(alpha),imag(alpha),'pk','MarkerSize',8,'MarkerFaceColor','k');
% plot(real(z2),imag(z2),'dk','MarkerSize',6,'MarkerFaceColor','k');  
% text(10.35,2,{'$z_2$'},'Interpreter','LaTeX','FontSize',20)
% text(10.35,imag(alpha),{'$\alpha$'},'Interpreter','LaTeX','FontSize',20)
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
print -dpdf HypSqu_2.pdf
%%
% cap   = the capacity of the domain G
tic
[~,cap] =  annq (et,etp,n,alpha,z2,'b');
toc
% hd    = the hyp diameter of the set E with respect to the outer polygon
tic
hd      = hypdiam(eto,etop,n,z2,zet.');
toc
'        h       cap         hd    '
myvalues=[h cap hd ]
% eval(['print -dpdf figure' num2str(ppp)])
% fig
%%