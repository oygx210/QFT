function [cab,z,p]=ldlgcplx(phi,wd,w,T)
% LDLGCPLX Lead/lag frequency response. (Utility Function)
%          LDLGCPLX computes the frequency response of a lead/lag component.
%          The function is also used to return the zero and pole
%          value of a specific phase and frequency for CNTDISP.

% Author: Craig Borghesani
% Date: 9/6/93
% Revised: 2/17/96 9:32 AM V1.1 update
% Copyright (c) 2003, Terasoft, Inc.


%%%%%% V5 change to accomodate nargin change
nargval = nargin;
cab = [];

if nargval==4,
 if T == 0,
  nargval=3;
 end
else
   T = 0;
end

if nargval==3,
 alpha=(1-sin(abs(phi*pi/180)))/(1+sin(abs(phi*pi/180)));
 a=wd*sqrt(alpha);
 b=wd/sqrt(alpha);
 if nargout==1,
  ca=rlroot(a,w,[sign(phi), T]); cb=rlroot(b,w,[-sign(phi), T]);
  cab=ca.*cb;
 else
  if sign(phi)<0,
   z=num2str(b,4);
   p=num2str(a,4);
  else
   z=num2str(a,4);
   p=num2str(b,4);
  end
 end
else
 x=cos(wd*T)+sin(phi*pi/180);
 y=sin(phi*pi/180)*cos(wd*T)+1;
 a=-(-y+sqrt(y^2-x^2))/x;
 b=-((sin(phi*pi/180)-a)/(-a*sin(phi*pi/180)+1));
 ac=-log(a)/T; bc=-log(b)/T;
 if nargout==1,
  ca=rlroot(ac,w,[1 T]); cb=rlroot(bc,w,[-1 T]);
  cab=ca.*cb;
 else
  z=num2str(a,8);
  p=num2str(b,8);
 end
end
