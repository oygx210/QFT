function sys = qsnsys(a,b,c,d,Cx,Cv,Ca)
% QSNSYS Transforms a system to the (S,n) mu-format. (Utility Function)
%        Input of state-space matrices:     sys=QSNSYS(A,B,C,D)
%              of the SN form:              sys=QSNSYS(S,n)
%                 the DE form:              sys=QSNSYS(D,E)
%              or the DF form:              sys=QSNSYS(D,F)
%              of mechanical matrices:      sys=QSNSYS(M,L,K,F,Cx,Cv,Ca)
%              or of mu-formatted systems:  sys=QSNSYS(sys)
%                                               QSNSYS(syse)
%                                               QSNSYS(sysf)
%                                               QSNSYS(sysm)
%
%        [S,n]=UNPEK(sys) returns the (S,n) form

% Author: Pepijn Wortelboer
% 8/31/93
% Copyright (c) 2003, Terasoft, Inc.


if nargin>4
  M=a; L=b; K=c; F=d;
  [n,nu]=size(F);
  A=[zeros(n,n) eye(n);-M\K -M\L];
  B=[zeros(n,nu);M\F];
  C=[Cx-Ca*(M\K) Cv-Ca*(M\L)];
  D=Ca*(M\F);
  a=[A B;C D];
  b=2*n;
elseif nargin==4
  n=length(a);
  a=[a b;c d];
  b=n;
elseif nargin==2
  [ny,nu]=size(a);
  [n,nn]=size(b);
  if nn==1+nu+ny   % DF form
    [a,b]=qdf2sn(a,b);
  else             % DE form
    [a,b]=qde2sn(a,b);
  end
elseif nargin==1
  [nr,nc]=size(a);
  if a(nr,nc)~=-Inf
%    help snsys
%    error('sys is not a system in mu-format')
     return;
  end
  n=a(nr,nc-1);
  if n~=0      % DF or DE form
    [a,b]=qunpack(a);
    [ny,nu]=size(a);
    [n,nn]=size(b);
    if nn==1+nu+ny   % DF form
      [a,b]=qdf2sn(a,b);
    else             % DE form
      [a,b]=qde2sn(a,b);
    end
  else
    if a(1,nc)<0
      [M,L,K,F,Cx,Cv,Ca]=qunpack(a);
      [n,nu]=size(F);
      A=[zeros(n,n) eye(n);-M\K -M\L];
      B=[zeros(n,nu);M\F];
      C=[Cx-Ca*(M\K) Cv-Ca*(M\L)];
      D=Ca*(M\F);
      a=[A B;C D];
      b=2*n;
    else
      sys=a;
      return
    end
  end
end
[nr,nc]=size(a);
sys=zeros(nr+1,nc+1);
sys(nr+1,nc+1)=-Inf;
if nr>0
  sys(1:nr,1:nc)=a;
  sys(1,nc+1)=b;
end
