function cpl=qcntbode(cont,w,T)
% QCNTBODE Controller frequency response. (Utility Function)
%          QCNTBODE computes the frequency response of a the elements stored
%          in a controller matrix.

% Author: Craig Borghesani
% Date: 9/3/93
% Revised: 2/17/96 9:33 AM V1.1 updates
% Copyright (c) 2003, Terasoft, Inc.
%       $Revision: 1.5 $

% Because all the CAD functions use this function, the following code
% detemines whether it is a continuous or discrete CAD function
% calling.  T=[] in all the continuous CAD environments

%%%%%% V5 change to accomodate nargin change
nargval = nargin;

if nargval==3,
 if T == 0,
  nargval=2;
 end
else
 T = 0;
end

[r,c]=size(cont);
cpl=ones(1,length(w));
for k=1:r,
 if nargval==2,
  if cont(k,4)==0,
   cp=cont(1,1);
  elseif (cont(k,4)==0.7 & (cont(k,1)~=0 | cont(k,2)~=0)),
   cp=cintegtr(cont(k,1),w) .* cintegtr(-cont(k,2),w);
  elseif (cont(k,4)==1 | cont(k,4)==2),
   cp=rlroot(cont(k,1),w,[(cont(k,4)==2)-(cont(k,4)==1), T]);
  elseif (cont(k,4)==3 | cont(k,4)==4),
   cp=cproot(cont(k,1),cont(k,2),w,[(cont(k,4)==4)-(cont(k,4)==3), T]);
  elseif cont(k,4)==5,
   cp=ldlgcplx(cont(k,1),cont(k,2),w);
  elseif cont(k,4)==6,
   cp=ntchcplx(cont(k,1),cont(k,2),cont(k,3),w);
  elseif cont(k,4)==7,
   cp=complexlead(cont(k,1), cont(k,2), cont(k,3), w);
  else cp=1; end
 else
  if cont(k,4)==0,
   cp=cont(1,1);
  elseif (cont(k,4)==1 | cont(k,4)==2),
   cp=rlroot(cont(k,1),w,[(cont(k,4)==2)-(cont(k,4)==1) T]);
  elseif (cont(k,4)==3 | cont(k,4)==4),
   cp=cproot(cont(k,1),cont(k,2),w,[(cont(k,4)==4)-(cont(k,4)==3) T]);
  elseif (cont(k,4)==0.5 & cont(k,1)~=cont(k,2)),
   ca=1; cb=1;
   if cont(k,1)~=0, ca=dintegtr(cont(k,1),w,T,-1); end
   if cont(k,2)~=0, cb=dintegtr(cont(k,2),w,T,1); end
   cp=ca.*cb;
  elseif (cont(k,4)==0.6 & (cont(k,1)~=0 | cont(k,2)~=0)),
   cp=cintegtr(cont(k,1),w,T) .* cintegtr(-cont(k,2),w,T);
  elseif cont(k,4)==5,
   cp=ldlgcplx(cont(k,1),cont(k,2),w,T);
  elseif cont(k,4)==6,
   cp=ntchcplx(cont(k,1),cont(k,2),cont(k,3),w,T);
  elseif cont(k,4)==7,
   cp=complexlead(cont(k,1), cont(k,2), cont(k,3), w, T);
  else cp=1; end
 end
 cpl=cpl.*cp;
end
