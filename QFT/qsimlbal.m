function [Gr,hsv,Tr,Tl]=qsimlbal(G,P,Q)
% QSIMLBAL Minimal PQ-balanced realization of G. (Utility Function)
%          [GR,HSV,TR,TL] = qsimlbal(G,P,Q) finds the minimal PQ-balanced
%          realization of G.
%          G  : system in mu-format with
%          P,Q: Gramians related to the realization of G
%          GR : minimal PQ-balanced realization of order nr
%          HSV: PQ-Hankel singular values
%          TR : PQ-balancing transformation
%          TL : TL*TR=I, TL=inv(TR) iff both P and Q are of rank n.
%          Let n be the order of G and nr be the order of Gr then
%          TL*P*TL'=TR'*Q*TR=diag(HSV(1:nr))

% Author: Pepijn Wortelboer
% 8/31/93
% Copyright (c) 2003, Terasoft, Inc.


[A,B,C,D]=qunpack(G);
n=length(A);

p=qsym2def(P)';
q=qsym2def(Q)';

[nrp,ncp]=size(p);
[nrq,ncq]=size(q);
nr=min(nrp,nrq);

%[Tr,hsv,Tl]=svd0(p*q');
[U,s,V] = svd(p*q',0);
if min(size(s))>1
  s=diag(s);
end
i=find(s>0);
Tr=U(:,i);
hsv=s(i);
Tl=V(:,i);

Tr=p'*Tr;
Tl=q'*Tl;
[hsv,index]=sort(-hsv);
hsv=-hsv;
Tr=real(Tr(:,index));
Tl=real(Tl(:,index));
vsh=sign(diag(Tr(1:nr,1:nr)))./sqrt(hsv);

if ~length(vsh),
 Gr = 0;
else
 Tr=Tr.*vsh(:,ones(1,n)).';
 Tl=Tl.*vsh(:,ones(1,n)).';
 Gr=qsnsys([Tl'*[A*Tr B];C*Tr D],nr);
end
