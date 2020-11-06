%FILE: maketableA.m
% For hypdia20201028.tex
% Basedon:Capacity_special_sets_table.m 
%Based on: hyp_Reuleaux_triangle_v6.m
% 2020-09-21....2020-10-06
% USES: hyR_tri.m, R_tri.m, widemarg.m
% Output: maketableA   and hypVsEuc2.pdf
clear
addpath bie fmm files
% rhoB is the hyp distance in the unit disk
rhoB=@(x,y)(2 *asinh(abs(x-y)./sqrt((1-abs(x).^2).*(1-abs(y).^2))));
%%
% The next formulas give capacities of special sets with hyp diam = t
% Jung's upper bound: from example15BLoop.m 2020-09-07
% based on Dekster's paper  
capJung = @(t) 2*pi./log((1+sqrt(1+(4/3)*sinh(t/2).^2))./(sqrt(4/3)*sinh(t/2)));
capDisk = @(t) 2*pi./log(1/tanh(t/4));
capSeg2 = @(t) 2*pi./mu(tanh(t/2));
%%
n      =  3*2^8; 
t      = (0:2*pi/n:2*pi-2*pi/n).';
% parametrization of the unit circle
etc    =  exp(i.*t);
etcp   =  i.*exp(i.*t);
%%
rv     = [0.05:0.10:0.95].';

for k=1:length(rv)
    r = rv(k);  
    ver    =  r.*exp(-i*2*pi*(0:2).'/3);
    %
    % alpha is a point in the annulus domain for both the hyp Reuleaux
    % triangle and the Euclidean Reuleaux triangle
    alpha  =(0.4+0.6*r).*exp(i*pi/3);
    z2     = 0;    
    %
    % 
    % parametrization of the hyp Reuleaux triangle
    %FOR ThinkPad:   [etr,etrp,hy_diam(k,1)]=hyR_tri(r,n/3);
    [ethr,ethrp]=hyR_tri(r,n/3);
    % parametrization of the whole boundary of the annulus domian between 
    % the unit circle and the hyp Reuleaux triangle
    eth     = [etc  ; ethr];
    ethp    = [etcp ; ethrp];
    %
    % capHRT   = the capacity of the hyp Reuleaux triangle
    [~,capHRT(k,1)] =  annq (eth,ethp,n,alpha,z2,'b');
    %
    % 
    % parametrization of the Euclidean Reuleaux triangle
    [eter,eterp]=R_tri(r,n/3);
    % parametrization of the whole boundary of the annulus domian between 
    % the unit circle and the Euclidean Reuleaux triangle
    ete     = [etc  ; eter];
    etep    = [etcp ; eterp];
    %
    % capHRT   = the capacity of the Euclidean Reuleaux triangle
    [~,capERT(k,1)] =  annq (ete,etep,n,alpha,z2,'b');
    %
    % hy_diam(k,1) is equal to rhoB(ver(1), ver(2))
    hy_diam(k,1)=rhoB(ver(1), ver(2));
    ec_rad  =  tanh(hy_diam(k,1)/4);
    ec_diam =  2*ec_rad;
    xt      =  ec_rad.*exp(-i.*t);    
    %
    % capRing = the capacity of the ring
    capRing(k,1)   = 2*pi/log(1/ec_rad);    
    %
    % capSeg = the capacity of the segment [-ec_rad,ec_rad]
    capSeg(k,1)   = capSeg2(hy_diam(k,1)); %OLD: 2*pi/mu(2*ec_rad./(1+ec_rad^2));
    Jung(k,1)=capJung(hy_diam(k,1));
    %
end
%%


%%
figure
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on; box on
t=hy_diam(1:8,1);
y1=capRing(1:8,1); % capERT capRing capHRT Jung
y2=capERT(1:8,1);
y3=capHRT(1:8,1);
y4=Jung(1:8,1);
y5=y2./y1;
y6=y5(t<4);
plot(t,y2./y1,'-sk','LineWidth',1.5)
plot(t,y3./y1,'-db','LineWidth',1.5)
plot(t,y4./y1,'-or','LineWidth',1.5)
plot(t,y1./y1,'--+k','LineWidth',1.5,'MarkerSize',8)
% plot(t(t<4),y6,':sr','LineWidth',1.5,'MarkerSize',8)
%
xlabel('hyperbolic diameter {$t$}','Interpreter','latex');
% ylabel('??','FontSize',15,'Interpreter','latex');
legend({'capERtri{$/b_1(t)$}','capHRtri{$/b_1(t)$}','capJung{$/b_1(t)$}','{$b_1(t)/b_1(t)$}'},...
        'Interpreter','latex','Location','southwest');

grid on
grid('minor')
set(gca, 'XMinorTick','on')
set(gca, 'YMinorTick','on')
ax=gca;
ax.GridAlpha=0.5;
ax.MinorGridAlpha=0.5;
% widemarg(gcf)
set(gca,'FontSize',20)
set(gca,'LooseInset',get(gca,'TightInset'))
axis([0.0 4  0.7 1.2])
print -dpdf hypVsEuc_2.pdf