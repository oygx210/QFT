function [zeta,wn,qflag]=dsecond(dph,m,w,T)
% DSECOND Discrete-time second order. (Utility Function)
%         DSECOND computes the zeta and natural frequency when given the
%         desired change in magnitude and phase.

% Author: Yossi Chait
% 3/22/94
% Copyright (c) 2003, Terasoft, Inc.


qflag=0;

if dph>0,
 dph=-dph; m=1/m;
end
p=dph*pi/180;
c1=cos(w*T);  s1=sin(w*T);
c2=cos(p+w*T);  s2=sin(p+w*T);
c3=cos(p+2*w*T);  s3=sin(p+2*w*T);
c4=cos(p);  s4=sin(p) ;

n=-(c3-c1/m)*(s2-s1/m)+(c2-c1/m)*(s3-s1/m);
d=(c4-c1/m)*(s2-s1/m)-(c2-c1/m)*(s4-s1/m);
a=-log(n/d)/2/T;

cbt=(s3-s1/m+(n/d)*(s4-s1/m))/(2*sqrt(n/d)*(s2-s1/m));
b=acos(cbt)/T;

zeta=sqrt(a^2/(a^2+b^2));
wn=b/sqrt(1-zeta^2);

if imag(zeta)~=0 | imag(wn)~=0 | zeta==0,
 zeta=[]; wn=[];
 qflag = 1;
else
 a=zeta*wn; b=wn*sqrt(1-zeta^2);
 rts=roots([1 -2*exp(-a*T)*cos(b*T) exp(-2*a*T)]);
 if any(abs(rts)>=1),
  qflag = 2;
 end
end
