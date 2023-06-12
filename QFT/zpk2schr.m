function sys=zpk2schr(Z,P,k)
% ZPK2SCHR Zeros-poles-gain to real Schur system realization form.
%          SYS = ZPK2SCHR(Z,P,K) converts zeros-poles-gain format to real
%          Schur system realization.
%          Z: vector with zeros
%          P: vector with poles
%          K: gain
%                   (s-z1)(s-z2)(s-z3)(s-z4)
%          sys = ------------------------------ * k
%                (s-p1)(s-p2)(s-p3)(s-p4)(s-p5)

% Author: Pepijn Wortelboer
% 4/1/94
% Copyright (c) 2003, Terasoft, Inc.


AA=[];BB=[];CC=[];DD=k;
Z=cplxpair(Z);
P=cplxpair(P);
nZ=length(Z);
nP=length(P);
while nP>0
  p=P(1);
  if imag(p)==0
    np=1;
    P(1)=[];
  elseif P(2)==conj(p)
    np=2;
    P(1:2)=[];
  else
    error('problem with pole pair sorting');
  end
  if nZ>0
    z=Z(1);
    if imag(z)==0
      if np==2 & nZ>=2
        z2=Z(2);
        nz=-2;  % complex pole pair versus two real zeros
        Z(1:2)=[];
      else
        nz=1;
        Z(1)=[];
      end
    elseif Z(2)==conj(z)
      nz=2;
      Z(1:2)=[];
      if np==1
        p2=P(1);
        np=-2;  % complex zero pair versus two real poles
        P(1)=[];
      end
    else
      error('problem with zero pair sorting');
    end
  else
    nz=0;
  end
  if abs(np)==2
    if np==-2
      rez=real(z);
      zz=z*z';
      A=[p 1;0 p2];
      C=[zz+p*p-p*2*rez,-2*rez+p+p2];
      D=1;
      B=[0;1];
    else
      rep=real(p);
      imp=imag(p);
      pp=p*p';
      spp=sqrt(pp);
      A=[rep -imp;imp rep];
      if abs(nz)==2
        if nz==2
          rez=real(z);
          imz=imag(z);
          zz=z*z';
          D=1;
          C=[2*(rep-rez) zz/imp+rep*(rep/imp)-imp-2*rez*(rep/imp)];
          B=[1;0];
        else
          zz=z*z2;
          D=1;
          C=[2*rep-z-z2 zz/imp+rep*(rep/imp)-imp-(z+z2)*(rep/imp)];
          B=[1;0];
        end
      elseif nz==1
        B=[1;0];
        C=[1 (rep-z)/imp];
        D=0;
      elseif nz==0
        B=[1;0];
        C=[0 1/imp];
        D=0;
      end
    end
  elseif np==1
    A=p;
    if nz==1
      B=1;
      C=p-z;
      D=1;
    elseif nz==0
      B=1;
      C=1;
      D=0;
    end
  end
  nZ=length(Z);
  nP=length(P);
    % To keep states of G1 on top:
    % [ AA BB*C  | BB*D ]
    % [  0   A   | B    ]
    % -----------+-------
    % [ CC DD*C  | DD*D ]
  n=length(AA);
  if isempty(BB),
     BB = zeros(0,size(C,1));
  end;%if
  AA=[AA BB*C;zeros(abs(np),n) A];
  BB=[BB*D;B];
  CC=[CC DD*C];
  DD=DD*D;
end
sys=qsnsys(AA,BB,CC,DD);
