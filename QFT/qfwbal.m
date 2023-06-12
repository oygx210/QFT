function [sys,hsv]=qfwbal(r,p,k,fs,zorr)
% QFWBAL Computes the frequency weighted balanced realization. (Utility)
%        [SYS,HSV] = QFWBAL(R,P,K,FS,ZOOR)
%        R  : residues
%        P  : poles
%        D  : D-term
%        FS : packed frequency function definition: piecewise linear with
%        (m0 + m1 * w) between w1 and w2 stored as
%        [m0 m1 w1 w2 i]  i=1 for input, i=2 for output, i=3 for both
%        Multiple rows are allowed, [1 0 0 Inf i] for unit weight
%        default: [1 0 0 Inf 3]
%        SYS: frequency weighted balanced realization
%        HSV: frequency weighted Hankel singular values

% Author: Pepijn Wortelboer
% 8/31/93
% Copyright (c) 2003, Terasoft, Inc.


if nargin<5
  zorr='r';
end
if length(fs)==0
  fs=[1 0 0 Inf 3];
end
n=length(p);
[nfs,nc]=size(fs);
if zorr=='z'
  g=zpk2schr(r,p,k);
  [A,B,C,D]=qunpack(g);
  T=eye(n);
  [T,A]=rsf2csf(T,A);
  B=T\B;
  C=C*T;
  L=diag(A);
end
if zorr=='r'
  g = qrpk2gf(r,p,k);
  [D,F,nrcp]=qunpack(g);
  if sum(nrcp)==0
    sys=qsnsys(D,0);
    hsv=[];
    return
  end
  [n,nn]=size(F);
  ni=1;
  L=[]; B=[]; C=[];
  if n>0
    L =F(:,1);
    B=F(:,2:2+ni-1);
    C=F(:,2+ni:nn)';
  end
  n=length(L);
  [E,Vi,V] = qf2e(F,nrcp);
  g=[E;D cumsum(nrcp.*[1 2]) -Inf];
end


P=zeros(n); Q=P;
for ii=1:nfs
  if fs(ii,[2 3 4])==[0 0 Inf]
    if zorr=='r'
      X=-L(:,ones(1,n))-L(:,ones(1,n))';
      if fs(ii,5)==1 | fs(ii,5)==3
        P=P+fs(ii,1)*(B*B')./X;
      end
      if fs(ii,5)==2 | fs(ii,5)==3
        Q=Q+fs(ii,1)*(C'*C)./conj(X);
      end
    else
      if fs(ii,5)==1 | fs(ii,5)==3
        P=P+qlyaps('ul',A,A',fs(ii,1)*B*B');
      end
      if fs(ii,5)==2 | fs(ii,5)==3
        Q=Q+qlyaps('lu',A',A,fs(ii,1)*C'*C);
      end
    end
  else
    if zorr=='r'
      nL=length(L);
      Li=L(:,ones(1,nL)); Lj=L(:,ones(1,nL))';
      E1=ones(nL);
      oma=fs(ii,3); omb=fs(ii,4); pl=fs(ii,1:2);
      X1= (-j*(Li+Lj)).\(pl(1)*E1-j*Li*pl(2))...
          .*log((j*omb*E1-Li)./(j*oma*E1-Li))+...
          (j*(Li+Lj)).\(pl(1)*E1+j*Lj*pl(2))...
          .*log((j*omb*E1+Lj)./(j*oma*E1+Lj));
      oma=-fs(ii,4); omb=-fs(ii,3); pl=[1 -1].*fs(ii,1:2);
      X2= (-j*(Li+Lj)).\(pl(1)*E1-j*Li*pl(2))...
          .*log((j*omb*E1-Li)./(j*oma*E1-Li))+...
          (j*(Li+Lj)).\(pl(1)*E1+j*Lj*pl(2))...
          .*log((j*omb*E1+Lj)./(j*oma*E1+Lj));
      X=X1+X2;
      if fs(ii,5)==2 | fs(ii,5)==3
        Q=Q+(C'*C).*conj(X)/2/pi;
      end
      if fs(ii,5)==1 | fs(ii,5)==3
        P=P+(B*B').*X/2/pi;
      end
    else
      Sp=-j*logm((j*fs(ii,3)*eye(n)-A)\(j*fs(ii,4)*eye(n)-A));
      Sm=-j*logm((-j*fs(ii,4)*eye(n)-A)\(-j*fs(ii,3)*eye(n)-A));
      Ap=eye(n)*fs(ii,1)-j*A*fs(ii,2);
      Am=eye(n)*fs(ii,1)+j*A*fs(ii,2);
      if fs(ii,5)==2 | fs(ii,5)==3
        Xp=qlyaps('lu',A',A,Ap'*C'*C);
        Xm=qlyaps('lu',A',A,Am'*C'*C);
        Q=Q+(Xp'*Sp+Sp'*Xp+Xm'*Sm+Sm'*Xm)/2/pi;
      end
      if fs(ii,5)==1 | fs(ii,5)==3
        Xp=qlyaps('ul',A,A',B*B'*Ap');
        Xm=qlyaps('ul',A,A',B*B'*Am');
        P=P+(Xp*Sp'+Sp*Xp'+Xm*Sm'+Sm*Xm')/2/pi;
      end
    end
  end
end
if zorr=='r'
  P=xxxpep(xxxpep(V,P,'x#'),Vi,'#x');
  Q=xxxpep(xxxpep(V,Q,'x#'),Vi,'#x');
  P=real(P);
  Q=real(Q);
  [sys,hsv]=qsimlbal(qsnsys(g),P,Q);
else
  P=real(T*P*T');
  Q=real(T*Q*T');
  [sys,hsv]=qsimlbal(g,P,Q);
end
