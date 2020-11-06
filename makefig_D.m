% FILE: makefig5.m
% hypdia20201028.m
% Based on :cmpRtri.m
% USES: widemarg.m
% Output: cmpRtri2.pdf
% From Mathematica: 
% 2020-10-22
% capJung[t_] := 
% 2 Pi/Log[(1 + 
%       Sqrt[1 + (Sqrt[4/3] Sinh[t/2])^2])/ (Sqrt[4/3] Sinh[t/2])]
%Plot[{ capJung[t]/( 2 Pi/Log[1/Tanh[t/4]])}, {t, 0.000001, 15}, 
% GridLines -> Automatic]
clear 
close all
uuu=sqrt(4/3);
capJung= @(t) 2*pi./log((1+sqrt(1+(uuu^2) * sinh(t/2).^2))./(uuu * sinh(t/2)));
capDisk=@(t) 2*pi./log(1./tanh(t/4));
t=0.001:0.001:10;
plot(t,(2/sqrt(3))*ones(size(t)),'-b',t,capJung(t)./capDisk(t),'k-','LineWidth',1.5)
grid on
grid('minor')
% widemarg(gcf)
set(gca,'FontSize',20)
set(gca,'LooseInset',get(gca,'TightInset'))
axis([-0.5 10.5  1.0 1.17])
print -dpdf cmpRtri_2
