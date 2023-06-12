function gf = qrpk2gf(R,P,K)
% QRPK2GF SISO system to complex modal system realization. (Utility Function)
%         GF = QRPK2GF(R,P,K) makes a complex modal system realization for
%         a siso system in Residue form (see RESIDUE)
%         R : residues
%         P : poles
%         K : direct feedthrough term (D)
%         GF: complex modal form (scaled, sorted)

% Author: Pepijn Wortelboer
% 8/31/93
% Copyright (c) 2003, Terasoft, Inc.


i=find(R==0);
if length(i)>0
 R(i)=[];
 P(i)=[];
end
n=length(P);
S=sqrt(abs(R));
F=[P S conj(R./S)];
cr=[100*eps;1e-8];
X=F; v=imag(F(:,1));
nv=size(v); nX=size(X);
if nv(1)==nX(1) & nv(2)<=1
  [w,l]=sort(v);
  F=X(l,:);
elseif  nv(2)==nX(2) & nv(1)==1
  [w,l]=sort(v);
  F=X(:,l);
else
  return;
end
L=F(:,1);
p=atan2(real(L),imag(L));
ireal=sort([find(abs(abs(p)-.5*pi)<cr(1))
            find(real(L)==0&imag(L)==0)]);
nr=length(ireal);
ncp=(n-nr)/2;
if ncp>0
  ij=[1:ncp];
  ji=[n:-1:n-ncp+1];
  F(ji,1:3)=conj(F(ij,1:3));
end
if nr>0
  index=ncp+1:ncp+nr;
  X=F(index,:); v=real(F(index,1));
  nv=size(v); nX=size(X);
  if nv(1)==nX(1) & nv(2)<=1
    [w,l]=sort(v);
    F(index,:)=X(l,:);
  elseif  nv(2)==nX(2) & nv(1)==1
    [w,l]=sort(v);
    F(index,:)=X(:,l);
  else
    return;
  end
end

if (~length(K))
  gf(n+1,1) = 0;
else
  gf(n+1,1)=K;
end
gf(n+1,[2 3])=[nr+j*(ncp+1e-300) -Inf];
gf(1:n,:)=F;
