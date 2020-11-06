%FILE: maketableA.m
% For hypdia20201028.tex
% Basedon:Capacity_special_sets_table.m 
%Based on: hyp_Reuleaux_triangle_v6.m
% 2020-09-21....2020-10-06
% USES: hyR_tri.m, R_tri.m, widemarg.m
% Output: maketableA   and hypVsEuc2.pdf
close all
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
%
for k=1:length(rv)
    r = rv(k);  
    ver    =  r.*exp(-i*2*pi*(0:2).'/3);
    %
    % alpha is a point in the annulus domain for both the hyp Reuleaux
    % triangle and the Euclidean Reuleaux triangle
    alpha  =(0.4+0.6*r);%.*exp(i*pi/3);
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
fname= ['maketableA' datestr(now,29) '.txt'];
delete(fname)
diary(fname)
disp([ datestr(now,31)    '  FILE: ' fname])
%disp(datestr(now))
disp(['      r       h-diam    capSeg    capERtri  capDisk   capHRtri  capJung'] )
disp([rv hy_diam capSeg capERT capRing capHRT Jung])
diary off
%%