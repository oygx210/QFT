function qbodeplot
%QBODEPLOT Open bode plot window. (Utility)

% Author: Craig Borghesani
% 9/5/93
% Copyright (c) 2003, Terasoft, Inc.

f=gcf;
bthan=get(f,'userdata');
infmat=get(bthan(16),'userdata');
cont = get(bthan(19),'userdata');
lomat = get(bthan(20),'userdata');
nom = get(bthan(2),'userdata');
uL0 = get(bthan(8),'userdata');

if isempty(cont),
   cont = get(bthan(3),'userdata');
   lomat=get(bthan(1),'userdata');
end
T=get(bthan(13),'userdata');

w=lomat(1,:);
[z,p,k] = cnt2zpk(cont,T);

if T == 0,
   OpenLoopZPK = zpk(z,p,k);
else
   OpenLoopZPK = zpk(z,p,k,T);
end

if isempty(nom),
   OpenLoopZPK = OpenLoopZPK * frd(uL0, w);

else
   OpenLoopZPK = OpenLoopZPK * nom;

end

ClosedLoop1 = OpenLoopZPK / (1 + OpenLoopZPK);
ClosedLoop2 = 1 / (1 + OpenLoopZPK);

BodeFigure = findobj(allchild(0),'name','Bode Plots');
if isempty(BodeFigure),
   figure('name','Bode Plots');

else
   figure(BodeFigure);

end

bode(OpenLoopZPK, 'r', ClosedLoop1, 'g', ClosedLoop2, 'b');
legend('L_0','L_0/(1+L_0)','1/(1+L_0)');

