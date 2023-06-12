function bdb = qrobust(w,wbd,W,uP,vP,R,nom,uC,vC,ctype,ph_r)
% QROBUST Robust Stability Bound for R ~= 0. (Utility Function)
%         QROBUST determines the stability bound for the case when R is not
%         equal to zero for non-parametric plant cases.

% Author: Craig Borghesani
% Date: 10/5/92
% Revision: 2/18/96 12:28 AM V1.1 updates
% Copyright (c) 2003, Terasoft, Inc.


m=length(ph_r); nbd=length(w);
pbds=qsubset(wbd,w);
myeps=1e-16;

% preallocating matrix for more efficient execution
bmag=[ones(m,length(pbds))*myeps;ones(m,length(pbds))*1/myeps];

[rmp,cmp]=size(uP); [rmc,cmc]=size(uC);

% declaring sizes of replicating matricies
lo2=max(rmp,rmc);

if repltest,
 u=ones(lo2,1); v=ones(1,m);
else
 u=ones(1,lo2); v=ones(1,m);
end

[rW,cW]=size(W);
if rW == 1, W = W(u,:); end
if cW == 1, W = W(:,v); end

if rmp~=rmc, val=min(rmp,rmc); else val=1; end

j=1; pct=1; cct=1; Wct=1;
while j<=length(pbds),

 if cmp>1, pct=pbds(j); end
 if cmc>1, cct=pbds(j); end

 phi=ph_r-(vP(nom(1),pct)+vC(nom(2),cct));

 cnt=1;
 if rmp>rmc, rp=1:rmp; rc=cnt;
 elseif rmc>rmp, rp=cnt; rc=1:rmc;
 else rp=1:rmp; rc=1:rmc; end

 while cnt<=val,

%%%%%% V5 code
% Reason: V5 does "the right thing" with replicating matrices.  in V4.2,
% if a vector was 10 elements long and you indexed it with a boolean
% matrix of the same length, the identical matrix was returned.  in V5,
% a matrix of the first element is returned.  this is consistent
% behavior and the if-statement below makes sure that the
% vectors to be replicated are of length 1; which they will be if
% length(rp) == 1

  Wt=W(rp,pct); mP=uP(rp,pct); pP=vP(rp,pct);
  if length(rp) == 1,
     Wt=Wt(u); mP=mP(u); pP=pP(u);
  end

  mC=uC(rc,cct); pC=vC(rc,cct);
  if length(rc) == 1,
     mC=mC(u); pC=pC(u);
  end

  psi=pP+pC; psi=psi(:,v)+phi(u,:);

  A=mP.^2 .*mC.^2 .*(1-1 ./Wt.^2); A=A(:,v);
  B=2*mP(:,v).*mC(:,v).*(cos(psi));
  C=ones(lo2,m);

  [g1,g2]=quadrtic(A,B,C); size_g1 = size(g1);
  cbdb=[[[g1';g2'],bmag(:,j)];ones(2,size_g1(1)+1)];
  [abvblw,state(j)]=sectbnds(cbdb,0);
  bmag(:,j)=abvblw(1:2*m);
  if rmp>rmc, rc=rc+1; cnt=rc;
  elseif rmc>rmp, rp=rp+1; cnt=rp; else cnt=2; end
 end

 z=find(bmag(:,j)~=myeps & bmag(:,j)~=1/myeps & ...
        bmag(:,j)~=-248 & bmag(:,j)~=248 & ...
        bmag(:,j)~=-302 & bmag(:,j)~=302);
 if length(z), bmag(z,j)=bmag(z,j)*uP(nom(1),pct)*uC(nom(2),cct); end

 j=j+1;
end

bmag(bmag~=-248 & bmag~=248 & bmag~=-302 & bmag~=302) = ...
   20*log10(bmag(bmag~=-248 & bmag~=248 & bmag~=-302 & bmag~=302));
bdb=[bmag; w(pbds); ones(1,length(pbds))];
