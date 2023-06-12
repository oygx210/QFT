function [S,n] = qdf2sn(D,F)
% QDF2SN  [S,n] = qdf2sn(D,F)
%         transforms the complex (D,F) form
%                 to the complex (S,n) form
%         D : direct feedthrough matrix
%         F : [L Bl Cl'], the so-called LBC form. L is a column
%             vector with the system poles.
%         S : [Al Bl;Cl D], with Al=diag(L)
%         n : number of states
% If the system is already in (S,n) form, it is transferred to the
% output directly

% Author: Pepijn Wortelboer
% 8/31/93
% Copyright (c) 2003, Terasoft, Inc.


if size(F)==[1 1]
  S=D;
  n=F;
  return
end
[no,ni]=size(D);
[n,nn]=size(F);
if nn==2+no+ni
%  disp('input probably in (D,E) form instead of (D,F) form')
%  disp('Conversion achieved with F=e2f(F)')
%  F=e2f(F);
   E=F;
%   [n,nn]=size(E);
   if n==0
     return
   end
   nr=length(find(E(:,2)==0));
   ncp=(n-nr)/2;
% first construct Xtf2e
   Xtf2e=zeros(n,2);
   if ncp>0
     Xtf2e(1:ncp,:)=ones(ncp,1)*[j 1]/sqrt(2);
     Xtf2e(ncp+nr+1:n,:)=ones(ncp,1)*[1 -j]/sqrt(2);
   end
   if nr>0
     Xtf2e(ncp+1:ncp+nr,1)=ones(nr,1);
   end
% Xte2f can be derived easily since Te2f=Tf2e'
   Xtf2e1=Xtf2e(:,1);
   Xtf2e2=Xtf2e(:,2);
   F1=[E(:,1)+j*E(:,2), Xtf2e1(:,ones(1,nn-2)).*E(:,3:nn)+...
                    Xtf2e2(:,ones(1,nn-2)).*E(n:-1:1,3:nn)];
   ij=[1:ncp nr+ncp+1:n];
   ji=[n:-1:nr+ncp+1 ncp:-1:1];
   F1(ij,:)=(F1(ij,:)+conj(F1(ji,:)))/2;
   F=F1;

end
if n==0
  S=D;
else
  L=[]; Bl=[]; Cl=[];
  if n>0
    L =F(:,1);
    Bl=F(:,2:2+ni-1);
    Cl=F(:,2+ni:nn)';
  end
  S = [diag(L) Bl;Cl D];
end
