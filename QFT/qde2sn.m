function [S,n] = qde2sn(D,E)
% QDE2SN  [S,n] = qde2sn(D,E)
%         transforms the (D,E) form to the (S,n) form
%         D : direct feedthrough matrix
%         E : [X Bx Cx'], the so-called XBC form. X is a two column
%             matrix with real and imaginary part of the system poles.
%         S : [Ax Bx;Cx D], with Ax=diax(X) in the form of a cross.
%         n : number of states
% If the system is already in (S,n) form, it is transferred to the
% output directly

% Author: Pepijn Wortelboer
% 8/31/93
% Copyright (c) 2003, Terasoft, Inc.


X=[];  Bx=[];   Cx=[];
if size(E)==[1 1]
  S=D;
  n=E;
  return
end
[no,ni]=size(D);
[n,nn]=size(E);
if n==0
  S=D;
else
  X=E(:,1:2);
  Bx=E(:,3:2+ni);
  Cx=E(:,3+ni:nn)';
  Z=X;
  [nr,nc]=size(Z);
  if nc==2
    Y=diag(Z(:,2));
    nrc=ceil(nr/2);
    if nr~=2*nrc
      Y(nrc,nrc)=0;
    end
    Y=Y(:,nr:-1:1)+diag(Z(:,1));
  elseif nr==nc
    Y=[diag(Z) diag(Z(:,nc:-1:1))];
  else
    return;
  end
  S = [Y Bx;Cx D];
end
