function   [et,etp] = R_tri(rad,ns)
% polygonp.m
% Nasser, June 10, 2019
% This function compute the discretization of the parametrization of the
% Reuleaux triangle with the vertices "ver' where "ns" is the graded points
% on each side of the triangle. The number of graded points for the whole 
% triangle is "n=number of sides*ns". The graded mesh points are computed 
% by the function "deltw.m".
% 
% 
r      =  rad*sqrt(3);
n      =   3*ns;
t      =  (0:2*pi/n:2*pi-2*pi/n).';
[s,sp] =   deltw(t,3,3);
for j=1:3
    sv{j}  =  s((j-1)*n/3+1:j*n/3);
end
%
verT   =  rad.*exp(-i*2*pi*(0:2).'/3);
ver    = [verT;verT];
for j=1:3
    alp     = (j-1)*2*pi/3;
    bet     =     j*2*pi/3;
    a       =  Arg(ver(j)-ver(j-1+3));
    b       =  Arg(ver(j+1)-ver(j-1+3));
    sc{j}   = ((b-a).*sv{j}+a*bet-b*alp)/(bet-alp);
    scp{j}  = ((b-a).*ones(size(sv{j})))/(bet-alp);
    etv{j}  =  ver(j-1+3)+r.*exp(i.*sc{j});
    etvp{j} =  r.*scp{j}.*exp(i.*sc{j});
end
eto = []; etopo = [];
for j=1:3
    eto((j-1)*n/3+1:j*n/3,1)     =  etv{j};
    etopo((j-1)*n/3+1:j*n/3,1)   =  etvp{j};
end
et  =  eto; etp =  etopo.*sp;
end
%