function [pr,pc1,pc0]=basis2b(np0,dp0,bdb,w,ic,zeta,fac);
% sets up the poles in the basis for optimization
%

% 8.26.96 lpbasis2a.m
% 8.30.96 revision: match autot3b.m
% 8.30.96 revision: allow do-loop for selecting wbase

pc0 = []; pc1 = []; pr = [];

% allocate no. of real and complex basis poles
%
ipc0=length(dp0)-length(np0);  % pre-weight to insure properness, at least
% one real, one complex, one real....
dev=fix(ic/3);
rdev=ic-dev*3;
ipr=0; ipc1=0;
if dev>0,
   ipr=dev; ipc1=dev;
end
if rdev==1,
   ipr=ipr+1;
elseif rdev==2,
   ipc1=ipc1+1;
end

% find high-gain (> 0db) bound type closest to 0db
gain=bdb(1,:);
ig=find(gain>-200);
igain=find(gain(ig)==min(gain(ig)));
g1=min(gain(ig));
wbase=w(igain)*fac;

% find freq for lowest value below bound
[ir,ij]=size(bdb);
il=(ir-2)/2;
gain=bdb(il+1:2*il,ij);
ig=find(gain<200);
mngain=min(gain(ig));
wmn=bdb(ir-1,ij);
% assuming slope of -20 db/dec in this nhbd
dg=mngain;
whigh=wbase*10^(abs(dg/20));


index=1:.2:20;

% assign real poles
if ipr,
   pr=index(1:ipr)*wbase*1.05;
end
% assign complex poles with a zero upstairs
if ipc1,
   pc1=index(1:ipc1)*wbase;
   pc1=[zeta*ones(length(pc1),1),pc1'];
end

if ipc0,
    pc0=[1:ipc0]*whigh*5;
end

%pc1=[   6.4499e-001  1.5139e+001
%   6.2035e-001  4.3021e+001
%  5.0775e-001  5.7381e+002
%  5.9077e-001  2.2271e+003
%7.0004e-001  7.3283e+003];

%pc0=[1.3739e-004,7.1000e-001,6.0643e+000,5.8056e+001,1.7194e+002,2.1521e+002];

%pr= 1.8330e+3 ;

