function [contldlg,msg]=zp2ldlg(contzp,T)
% ZP2LDLG Zero-pole to lead/lag. (Utility Function)
%         ZP2LDLG converts zero-pole pairs to lead/lag terms.

% Author: Craig Borghesani
% Date: 10/10/93
% Revised: 2/17/96 11:26 PM V1.1 updates
% Copyright (c) 2003, Terasoft, Inc.


nargval = nargin;

if nargval==2,
 if T == 0,
  nargval=1;
 end
end

contldlg=[]; msg = [];
p=sort(contzp(find(contzp(:,4)==1),1));
z=sort(contzp(find(contzp(:,4)==2),1));
lz=length(z);
lp=length(p);

if nargval==1,
 if lz==lp & all([z(:);p(:)]>0),
  for k=1:lz,
   w = sqrt(z(k)*p(k));
   alph = p(k)/z(k);
   delph = ((p(k)>z(k))-(p(k)<z(k)))*abs(asin((1-alph)/(alph+1))*180/pi);
   contldlg=[contldlg;delph,w,NaN,5];
  end
 elseif any([z(:);p(:)]<0),
  if any(z<0), msg=3;
  else msg=4; end
 else
  if lz>lp, msg=1;
  else msg=2; end
 end
elseif nargval==2
 if lz==lp & all(abs([z(:);p(:)])>0),
  for k=1:lz,
   z = real(exp(-z(k)*T)); p = real(exp(-p(k)*T));
   delph=asin((z-p)/(1-z*p))*180/pi;
   c=(z-p)/(1-z*p); d=(z^2+1)/(2*z);
   w=acos((1-c*d)/(d-c))/T;
   contldlg=[contldlg;delph,w,NaN,5];
  end
 elseif any(abs([z(:);p(:)])<0),
  if any(abs(z)<0), msg=3;
  else msg=4; end
 else
  if lz>lp, msg=1;
  else msg=2; end
 end
end
