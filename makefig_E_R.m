% FILE: makefigE.m
% hypdia20201028.m
%Based on:  hyp_Reuleaux_triangle_v_n3.m
% USES: widemarg.m
%Output: EuHRtri2.pdf  :   Euclidean and hyperbolic Reuleaux triangle
close all
clear
addpath bie fmm files
%%
% Computing the hyperbolic distance
rhoB = @(x,y)(2*asinh(abs(x-y)./sqrt((1-abs(x).^2).*(1-abs(y).^2))));
%
%%
n      =  3*2^8; 
t      = (0:2*pi/n:2*pi-2*pi/n).';
%%
r     = 0.5;
ver    =  r.*exp(-i*2*pi*(0:2).'/3);
[etr,etrp,hy_diam,Eucent1,Eurad1]=hypReuleaux(r,n/3);
%
% parametrization of ring domain border by the unit disk and the hyp 
% Reuleaux triangle in the unit disk.
et     =  exp(i.*t);
etp    =  i.*exp(i.*t);
et     = [et;etr];
etp    = [etp;etrp];
% Choose alpha in the ring domain and z2 interior to the triangle
alpha  =(0.4+0.6*r).*exp(i*pi/3);
z2     = 0;
%
[Eucent1,Eurad1] = HypDisk(ver(1),rhoB(ver(1),ver(2)));
circ1 = Eucent1+Eurad1.*exp(2i*pi.*linspace(0,1,100));
%
xet     =  r.*exp(-i.*t);
%%  
figure
hold on;     box on
crv = circ1; 
k=1; crv = et((k-1)*n+1:k*n);crv(end+1) = crv(1); 
plot(real(crv),imag(crv),'k','LineWidth',1.5);
k=2; crv = et((k-1)*n+1:k*n);crv(end+1) = crv(1); 
plot(real(crv),imag(crv),'-r','LineWidth',1.5);
crv = xet;crv(end+1) = crv(1); 
plot(real(ver),imag(ver),'or','MarkerFaceColor','r');
hydia=rhoB(r,r*exp(i*2*pi/3));
circ2=tanh(hydia/4)*exp(i*2*pi*linspace(0,1,100));
plot(real(circ2), imag(circ2),'-.b','LineWidth',1.5);
eucRtri=r*sqrt(3)*exp(i*5*pi*linspace(1,1.4,100)/6)+r*ones(1,100);
plot( real(eucRtri), imag(eucRtri),'k--','LineWidth',1.5);
plot( real(exp(i*2*pi/3)*eucRtri), imag(exp(i*2*pi/3)*eucRtri),'k--','LineWidth',1.5);
plot( real(exp(i*4*pi/3)*eucRtri), imag(exp(i*4*pi/3)*eucRtri),'k--','LineWidth',1.5);
axis equal
grid on
% widemarg(gcf)
axis([-1.25 1.25 -1.15 1.15])
set(gca,'FontSize',20)
set(gca,'LooseInset',get(gca,'TightInset'))
%
print -dpdf EuHRtri_2.pdf
%%
