function [E,Xtf2e,Xinv] = qf2e(F,nrcp)
% QF2E   Transform from LBC to XBC form.
%        QF2E transforms the LBC form (diag(L),Bl,Cl) compressed in
%        F=[L Bl Cl'] to an XBC form (diax(X),Bx,Cx) compressed in
%        E=[X Bx Cx']
%        F : complex modal realization in compact form (LBC)
%        nrcp : [number of real poles, number of complex conjugate pole
%               PAIRS]
%        def: calculated by the program
%        E : real modal realization in compact form (XBC)
%        Xtf2e: compact notation of unitary transformation matrix:
%        Tf2e = diax(Xtf2e) , Te2f=Tf2e'

% Author: Pepijn Wortelboer
% 8/31/93
% Copyright (c) 2003, Terasoft, Inc.


[n,nn]=size(F);
if n==0
  return
end
if nargin==1
  nr=length(find(imag(F(:,1))==0));
  ncp=(n-nr)/2;
else
  nr=nrcp(1);
  ncp=nrcp(2);
end
Xtf2e=zeros(n,2);
if ncp>0
  Xtf2e(1:ncp,:)=ones(ncp,1)*[j 1]/sqrt(2);
  Xtf2e(ncp+nr+1:n,:)=ones(ncp,1)*[1 -j]/sqrt(2);
end
if nr>0
  Xtf2e(ncp+1:ncp+nr,1)=ones(nr,1);
end
Xte2f1=conj(Xtf2e(:,1));
Xte2f2=conj(Xtf2e(n:-1:1,2));
E=[real(F(:,1)),imag(F(:,1)),Xte2f1(:,ones(1,nn-1)).*F(:,2:nn)+...
                             Xte2f2(:,ones(1,nn-1)).*F(n:-1:1,2:nn)];
E=real(E);

X=Xtf2e;
[n,nc]=size(X);
if nc~=2
  error('X not a two-column matrix')
end
Xinv=zeros(n,2);
if n>1
  nr=length(find(imag(X(:,1))==0));
  ncp=(n-nr)/2;
  Xinv(:,1)=conj(X(:,1));
  Xinv(:,2)=conj(X(n:-1:1,2));
else
  nr=X(1);
  ncp=X(2);
  if ncp>0
    Xinv(1:ncp,:)=ones(ncp,1)*[j 1]/sqrt(2);
    Xinv(ncp+nr+1:n,:)=ones(ncp,1)*[1 -j]/sqrt(2);
  end
  if nr>0
    Xinv(ncp+1:ncp+nr,1)=ones(nr,1);
  end
end
