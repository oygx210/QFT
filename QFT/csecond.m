function [zeta,wn,flag]=csecond(dph,dmag,wo)
% CSECOND Continuous-time second order. (Utility Function)
%         CSECOND computes the zeta and natural frequency when given the
%         desired change in magnitude and phase.

% Author: Yossi Chait
% 3/21/94
% Copyright (c) 2003, Terasoft, Inc.


if all([dph,dmag]),
 flag=0;
 if dph>0,
  dph=-dph; dmag=1/dmag;
 end
 wn=wo/sqrt(1-cos(pi/180*dph)/dmag);
 zeta=(wo^2-wn^2)*tan(pi/180*dph)/2/wo/wn;
 if imag(wn)~=0 | zeta<=0,

%%%% what the hell was this line here for??
%  wn=[]; zeta=[];
  if imag(wn)~=0,
   flag=1;
  elseif zeta<=0,
   flag=2;
  end
 end
elseif dph==0,
 flag=0;
 zeta = 0.707;
 wn = (dmag^2*wo^2*(2*zeta^2-1) + dmag*wo^2*sqrt(dmag^2*(2*zeta^2-1)^2 - dmag^2 + 1))/abs(dmag^2-1);
 wn = sqrt(wn);
elseif dmag==0,
 flag=0;
 zeta = 0.707;
 sign=1;
 P = tan(abs(dph)*pi/180);
 if abs(dph) > 90, sign = -1; end
 wn = wo*(zeta + sign*sqrt(zeta^2 + P^2))/P;
end
