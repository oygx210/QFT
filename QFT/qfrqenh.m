function wnew = qfrqenh(wn,zeta,wold,T)
% QFRQENH Frequency enhancement. (Utility Function)
%         QFRQENH enchances the present frequency vector when a second order
%         pole or zero is added.  This is accomplished by adding frequencies
%         around the damped natural frequency.

% Author: Craig Borghesani
% Date: 9/5/93
% Revised: 2/17/96 9:34 AM V1.1 updates
% Copyright (c) 2003, Terasoft, Inc.


%%%%%% V5 change to accomodate nargin change
nargval = nargin;

if nargval==4,
 if T == 0,
  nargval=3;
 end
end

%%%%%% V5 change because sometimes zeta == [] and the zeta == 0 statement
%%%%%% was causing an annoying warning.
if ~length(zeta),
   wnew = wold;
   return;
end

if zeta==0,
 zeta=1e-6;
end

wd=abs(wn)*sqrt(1-zeta^2);
if nargval==3,
 wpre=logspace(log10(wd/2),log10(.89*wd),25);
 wwd=logspace(log10(0.9*wd),log10(1.1*wd),50);
 wpost=logspace(log10(1.11*wd),log10(2*wd),25);
 wnew=sort([wpre wd wwd wpost wold]);
else
 ws=pi/T;
 x=wd/2; y=0.89*wd;
 if x>ws, wpre=[];
 elseif y>ws, y=ws; wpre=logspace(log10(x),log10(y),25);
 else wpre=logspace(log10(x),log10(y),25); end
 x=0.9*wd; y=1.1*wd;
 if x>ws, wwd=[];
 elseif y>ws, y=ws; wwd=logspace(log10(x),log10(y),50);
 else wwd=logspace(log10(x),log10(y),50); end
 x=1.11*wd; y=2*wd;
 if x>ws, wpost=[];
 elseif y>ws, y=ws; wpost=logspace(log10(x),log10(y),25);
 else wpost=logspace(log10(x),log10(y),25); end
 wnew=sort([wpre wd wwd wpost wold]);
end
