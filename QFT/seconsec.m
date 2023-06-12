function [contnew,flag]=seconsec(delph,delmag,w,T)
% SECONSEC Identify a super-second term. (Utility Function)
%          SECONSEC identifies a rational transfer function G(s) that takes
%          on the complex value p at s=jw, w not 0.
%
%          note: G(0) is set to 1 regardless of p(0).

% Author: Yossi Chait
% 12/06/94
% Integrated into toolbox 12/07/94 CWB
% Copyright (c) 2003, Terasoft, Inc.


% 1st try 2nd/2nd term

contnew=[];
flag=0;

epsq=0.01;
if T > 0,                             % discrete-time
 wm=exp(-w*T);  wm1=exp(-(epsq+w)*T);
 p1=qcpqft(1,[1,-wm],w,T);
 p2=qcpqft(1,conv([1,-wm],[1,-wm1]),w,T);
else                                      % continuous-time
 wm=w;  wm1=epsq+w;
 p1=qcpqft(1,[1/wm,1],w);
 p2=qcpqft(1,conv([1/wm,1],[1/wm1,1]),w);
end
s=[p1;p2];

rs=real(s);  is=imag(s);
A=[[1;0],[rs';is']];
if T > 0, A=[[1,1/(1-wm),1/(1-wm)/(1-wm1)];A];
else A=[ones(1,3);A]; end

p=delmag*exp(i*(delph*pi/180));
rp=real(p); ip=imag(p);
B=[1;rp';ip'];

x=A\B;

% form num/den of result
te = [1, (1-2*(T > 0))*wm];  te1 = [1, (1-2*(T > 0))*wm1];
if T > 0,
 nn2=[x(1)*conv(te,te1)+[0,x(2)*te1]+[0,0,x(3)]]; dd2=conv(te,te1);
 zero2=roots(nn2);
else
 nn2=[x(1)*conv(te,te1)+[0,x(2)*wm*te1]+[0,0,x(3)*wm*wm1]]; dd2=conv(te,te1);
 zero2=roots(nn2);
end
zero2=conj(zero2');
% now try simple 2nd order solution
if T > 0,
 [zeta,wn,state]=dsecond(delph,delmag,w,T);
 if state == 0,
   a=zeta*wn; b=wn*sqrt(1-zeta^2);
   nn1=(1-2*cos(b*T)*exp(-a*T)+exp(-2*a*T))*[1,0,0];
   dd1=[1,-2*cos(b*T)*exp(-a*T),exp(-2*a*T)];
 end

else
 [zeta,wn,state]=csecond(delph,delmag,w);
 nn1=wn^2; dd1=[1,2*zeta*wn,wn^2];
end

if delph>0 & state == 0,                               % pole or zero type?
 temp=dd1; dd1=nn1; nn1=temp;
end

% assign "best" solution
nn=[]; dd=[]; in=[];
if T > 0,                             % find rhp zeros
 in=find(abs(zero2) >=1);
else
 in=find(real(zero2) >=0);
end
if state~=0 & length(in)~=0,
 flag=2;                                  % detected unstable or rhp terms!
elseif state~=0,                          % simple 2nd is "bad"
 nn=nn2; dd=dd2;
elseif length(in)~=0,                     % 2/2 is "bad"
  nn=nn1; dd=dd1;
elseif state==0 & length(in)==0,          % both are "ok"
  nn=nn1; dd=dd1;
  if zeta<0.5,
   nn=nn2; dd=dd2;
  end
end

if flag==0,
 contnew = cntpars(tf(nn, dd, T));
end
