function dcgain=cntdcgn(cont,T)
% CNTDCGN Controller D.C. Gain. (Utility Function)
%         CNTDCGN determines the dc-gain of a controller matrix.

% Author: Craig Borghesani
% Date: 10/27/93
% Revised: 2/18/96 12:08 AM V1.1 updates
% Copyright (c) 2003, Terasoft, Inc.


nargval = nargin;

if nargval==2,
 if T == 0,
  nargval=1;
 end
end

dcgain=1;
if nargval==1,
 for p=1:length(cont(:,1)),
  if cont(p,4)==1,
   dcgain=dcgain/cont(p,1);
  elseif cont(p,4)==2,
   dcgain=dcgain*cont(p,1);
  elseif cont(p,4)==3,
   if cont(p,2)~=0, dcgain=dcgain/cont(p,2)^2; end
  elseif cont(p,4)==4,
   if cont(p,2)~=0, dcgain=dcgain*cont(p,2)^2; end
  elseif cont(p,4)==5,
   phi=cont(p,1); wd=cont(p,2);
   alpha=(1-sin(abs(phi*pi/180)))/(1+sin(abs(phi*pi/180)));
   a=wd*sqrt(alpha); b=wd/sqrt(alpha);
   if sign(phi)<0, dcgain=dcgain*a/b;
   else dcgain=dcgain*b/a; end
  end
 end
else
 for p=1:length(cont(:,1)),
  if cont(p,4)==1,
   dcgain=dcgain/(1-real(exp(-cont(p,1)*T)));
  elseif cont(p,4)==2,
   dcgain=dcgain*(1-real(exp(-cont(p,1)*T)));
  elseif cont(p,4)==3,
   a=cont(p,1)*cont(p,2);
   b=cont(p,2)*sqrt(1-cont(p,1)^2);
   dcgain=dcgain/(1-2*exp(-a*T(1))*cos(b*T(1))+exp(-2*a*T(1)));
  elseif cont(p,4)==4,
   a=cont(p,1)*cont(p,2);
   b=cont(p,2)*sqrt(1-cont(p,1)^2);
   dcgain=dcgain*(1-2*exp(-a*T(1))*cos(b*T(1))+exp(-2*a*T(1)));
  elseif cont(p,4)==5,
   phi=cont(p,1); wd=cont(p,2);
   x=cos(wd*T)+sin(phi*pi/180); y=sin(phi*pi/180)*cos(wd*T)+1;
   a=-(-y+sqrt(y^2-x^2))/x;
   b=-((sin(phi*pi/180)-a)/(-a*sin(phi*pi/180)+1));
   dcgain=dcgain*(1-b)/(1-a);
  end
 end
end
