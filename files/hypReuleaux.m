function   [et,etp,hyprad,Eucent,Eurad]=hypReuleaux(r,ns)
% hyR_tri.m
% Nasser, September 21, 2020
% Based on, Matti, Reuleaux_triangle.m, 2020-09-18
% 
% 
verT   =  r.*exp(-i*2*pi*(0:2).'/3);
ver    = [verT;verT];
% 
% Computing the hyperbolic distance
rhoB = @(x,y)(2*asinh(abs(x-y)./sqrt((1-abs(x).^2).*(1-abs(y).^2))));
% The hyperbolic rad 
hyprad=rhoB(verT(1),verT(2));
% The Euclidean center 
[Eucent,Eurad] = HypDisk(r,hyprad);
cv = [exp(2i*pi/3) ; 1 ; exp(-2i*pi/3)].*Eucent;
%
% Find the parametrization of the arc C_k
n      =  3*ns;
t      = (0:2*pi/n:2*pi-2*pi/n).';
[s,sp] =  deltw(t,3,3);
for j=1:3
    sv{j}  =  s((j-1)*n/3+1:j*n/3);
end
%
for j=1:3
    alp     = (j-1)*2*pi/3;
    bet     =     j*2*pi/3;
    a       =  Arg(ver(j)-cv(j));
    b       =  Arg(ver(j+1)-cv(j));
    sc{j}   = ((b-a).*sv{j}+a*bet-b*alp)/(bet-alp);
    scp{j}  = ((b-a).*ones(size(sv{j})))/(bet-alp);
    etv{j}  =  cv(j)+Eurad.*exp(i.*sc{j});
    etvp{j} =  Eurad.*scp{j}.*exp(i.*sc{j});
end
eto = []; etopo = [];
for j=1:3
    eto((j-1)*n/3+1:j*n/3,1)     =  etv{j};
    etopo((j-1)*n/3+1:j*n/3,1)   =  etvp{j};
end
et  =  eto; etp =  etopo.*sp;
end
%