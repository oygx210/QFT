function cab=ntchcplx(zta1,zta2,wn,w,T);
% NTCHCPLX Notch frequency response. (Utility Function)
%          NTCHCPLX computes the frequency response of a notch component.

% Author: Craig Borghesani
% Date: 9/6/93
% Revised: 2/17/96 9:44 AM V1.1 updates
% Copyright (c) 2003, Terasoft, Inc.


nargval = nargin;

if nargval==5,
 if T == 0,
  nargval=4;
 end
else
 T = 0;
end

if nargval==4,
 ca=cproot(zta1,wn,w,[1, T]);
 cb=cproot(zta2,wn,w,[-1, T]);
else
 ca=cproot(zta1,wn,w,[1 T]);
 cb=cproot(zta2,wn,w,[-1 T]);
end
cab=ca.*cb;
