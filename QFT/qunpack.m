function   [S,n,C,D,Cx,Cv,Ca]=qunpack(sys)
% QUNPACK    [S,n]=unpack(sys)
%
%   Unpacks a SYSTEM matrix in mu-format into its equivalent
%   (S,n), (D,E), (D,F), (D,H)  or (M,L,K,F,Cx,..) form
%     [S,n]=qunpack(sysSn);
%     [D,E,nrcp]=qunpack(sysDE);   % with nrcp [#real-poles,
%     [D,F,nrcp]=qunpack(sysDF);   % #complex-conjugate-pole-pairs]
%     [D,H,nrcp]=qunpack(sysDH);
%
%     [M,L,K,F,Cx,..]=qunpack(sysMLK);   % for mechanical systems
%                                        stored with SOSYS
%     [A,B,C,D]=qunpack(sys);    % transforms any system to ABCD form
%

% Author: Pepijn Wortelboer
% 8/31/93
% Copyright (c) 2003, Terasoft, Inc.


[nr,nc]=size(sys);
if sys(nr,nc)~=-Inf
  error('sys should be a mu-format system');
elseif nr==1 & nc==1
  S=[]; n=0;
  return
end
N=sys(nr,nc-1);
if N==0  % MLK or Sn form
  n=sys(1,nc);
  if n<0   % MLK form: mechanical system
    n=-n;
    nu=nc-3*n-1;
    ny=nr-n-1;
    i=1:n;
    M=sys(i,i);
    L=sys(i,i+n);
    K=sys(i,i+2*n);
    F=sys(i,(3*n+1):(nc-1));
    Ca=sys(n+1:nr-1,i);
    Cv=sys(n+1:nr-1,i+n);
    Cx=sys(n+1:nr-1,i+2*n);
    if nargout==4
      K=-M\K;
      L=-M\L;
      S=[zeros(n) eye(n);K L]; %A
      n=[zeros(n,nu);M\F];     %B
      C=[Cx+Ca*K Cv+Ca*L];
      D=zeros(ny,nu);
    else
      S=M;
      n=L;
      C=K;
      D=F;
    end
  else      % Sn form
    S=sys(1:nr-1,1:nc-1);
    if nargout==4
      D=S(n+1:nr-1,n+1:nc-1); C=S(n+1:nr-1,1:n);
      B=S(1:n,n+1:nc-1); S=S(1:n,1:n); n=B;
    end
  end
else        % DE, DF or DH form
  if imag(N)>0  % DF or DH form
    nrcp=[real(N) round(imag(N))];
    nr_cp=nrcp(1)+nrcp(2);
    n=nrcp(1)+2*nrcp(2);
    C=nrcp;
    if any(imag(sys(1:nr_cp,1))<0)  % DF form
      systype='DF';
      nu=nc-1-nr+n;
      S=sys(n+1:nr,1:nu);
      n=sys(1:n,1:nc);
    else                            % DH form
      systype='DH';
      nu=nc-1-nr+nr_cp;
      S=sys(nr_cp+1:nr,1:nu);
      n=sys(1:nr_cp,1:nc);
    end
    if nargout==4
      if all(systype=='DH')
        n=h2f(n,nu);
      end
      D=S;
      F=n; ni=nu;
      S=[]; n=[]; C=[];
      [n1,nn]=size(F);
      if n1>0
        S =F(:,1);
        n=F(:,2:2+ni-1);
        C=F(:,2+ni:nn)';
      end
      S=diag(S);
    end
  else                % DE form
    n=sys(nr,nc-1);
    nrcp=sys(nr,nc-2);
    C=[nrcp (n-nrcp)/2];
    nu=nc-2-nr+n;
    S=sys(n+1:nr,1:nu);
    n=sys(1:n,1:nc);
    if nargout==4
      D=S;
      E=n; ni=nu;
      S=[]; n=[]; C=[];
      [n1,nn]=size(E);
      if n1>0
        S=E(:,1:2);
        n=E(:,3:2+ni);
        C=E(:,3+ni:nn)';
      end
      Z=S;
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
      S=Y;
    end
  end
end
