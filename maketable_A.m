%FILE: hyp_diam_2cubes.m
% 2020-10-29
clear
addpath bie fmm files
% 
%%
% choose the value of a
a  = 0.2;
% vertices of the external polygon. The vertices are counterclockwise
% orriented
ver_out = [ 0 ; 3 ; 3+i ; 2+i ; 2+a*i ; 1+a*i ; 1+i ; i];
alpha   =  1.5+0.1i; % alpha must be in the domain E (the domain between
                       % the two ploygons.
z2      =  0.5+0.5i;    % z2 must be inside the inner polygon
%
n      =  2^13; 
t      = (0:2*pi/n:2*pi-2*pi/n).';
[eto,etop]=polygonp(ver_out,n/8);
%
hv = [0.1,0.2,0.3,0.4,0.45];
% 
%%
myvalues = [];
for kk=1:length(hv)
    h = hv(kk);
    % vertices of the internal polygon. The vertices are clockwise orriented
    ver_in  = 0.5+0.5i+h.*[-1-i ; -1+i ; 1+i ; 1-i];
    %
    [eti,etip]=polygonp(ver_in,n/4);
    zet =  eti;
    et  = [eto ; eti];
    etp = [etop; etip];
    %   
    % cap   = the capacity of the domain G
    tic
    [~,cap] =  annq (et,etp,n,alpha,z2,'b');
    toc
    % hd    = the hyp diameter of the set E with respect to the outer polygon
    tic
    hd      = hypdiam(eto,etop,n,z2,zet.');
    toc
    myvalues=[myvalues;h hd  cap];
end
%%
' h        hd       cap    '
myvalues
%%