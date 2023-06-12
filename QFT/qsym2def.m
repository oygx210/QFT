function [pp,pn]=qsym2def(P)
% QSYM2DEF [pp,pn] = qsym2def(P)
% factorizes a symmetric matrix:  P= pp*pp' - pn*pn'

% Author: Pepijn Wortelboer
% 8/31/93
% Copyright (c) 2003, Terasoft, Inc.


[nr,nc]=size(P);
if nr~=nc
  error('P not square');
else
  n=nr;
end
% P=U*S*V'=u1*s1*v1'+u2*s2*v2'+...
[U,S,V]=svd(P);
%
s=diag(S);
if0=find(s==0);
% P=Un*Sn*Un'
% Un=U*diag(teken)
[nr,nc]=size(U.'.*V');
X=[zeros(nr,1) U.'.*V'];
teken=sum(X')';
U=U.*teken(:,ones(1,n)).';
teken(if0)=zeros(length(if0),1);
ifp=find(teken>0);
ifn=find(teken<0);
if length(ifp)>0
  Up=U(:,ifp);
  sp=sqrt(s(ifp));
  pp=sp(:,ones(1,n)).'.*Up;
end
if length(ifn)>0
  Un=U(:,ifn);
  sn=sqrt(s(ifn));
  pn=sn(:,ones(1,n)).'.*Un;
end
