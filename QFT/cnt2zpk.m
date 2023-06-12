function [z,p,k]=cnt2zpk(contr,T)
% CNT2ZPK Controller matrix to zero/pole/gain format. (Utility Function)
%         CNT2ZPK converts a controller matrix into the zero/pole/gain
%         format.

% Author: Craig Borghesani
% 10/10/93
% Copyright (c) 2003, Terasoft, Inc.
%      $Revision: 1.4 $

nargval = nargin;

if nargval==2,
 if T == 0,
  nargval=1;
 end
else
   T = 0;
end

z=[]; p=[];
if nargval==1,
 contr=cntcvrt(contr);
 k=contr(1,1);
 for v=1:length(contr(:,1)),
  if contr(v,4)==0.7,
   if contr(v,1) > 0,
    p=[p;zeros(contr(v,1),1)];
   end

   if contr(v,2) > 0,
    z=[z;zeros(contr(v,2),1)];
   end

  elseif contr(v,4)==1,
   k=k*(contr(v,1));
   p=[p;-contr(v,1)];
  elseif contr(v,4)==2,
   k=k/(contr(v,1));
   z=[z;-contr(v,1)];
  elseif contr(v,4)==3,
   k=k*(contr(v,2)^2);
   p=[p;roots([1,2*contr(v,1)*contr(v,2),contr(v,2)^2])];
  elseif contr(v,4)==4,
   k=k/(contr(v,2)^2);
   z=[z;roots([1,2*contr(v,1)*contr(v,2),contr(v,2)^2])];
  end
 end
else
 nc=0; dc=0;
 k=contr(1,1);
 contr=cntcvrt(contr,T);
 for v=1:length(contr(:,1)),
  if contr(v,4)==0.5,
   p=[p;ones(contr(v,1),1)];
   if contr(v,1)==3,
      z=[z;-1];
   end
   z=[z;ones(contr(v,2),1)];
   if contr(v,2)==3,
      p=[p;-1];
   end
   if contr(v,1)>0,
      nc=nc+1;
      k = k * T^(contr(v,1));
   end
   if contr(v,2)>0,
      dc=dc+1;
      k = k / T^(contr(v,2));
   end
%  elseif contr(v,4)==0.6,
%   if contr(v,1)>0,
%    p=[p;zeros(contr(v,1),1)];
%   else
%    z=[z;zeros(-contr(v,1),1)];
%   end
  elseif contr(v,4)==1,
   k=k*(1-real(exp(-contr(v,1)*T)));
   p=[p;real(exp(-contr(v,1)*T))];
   nc=nc+1;
  elseif contr(v,4)==2,
   k=k/(1-real(exp(-contr(v,1)*T)));
   z=[z;real(exp(-contr(v,1)*T))];
   dc=dc+1;
  elseif contr(v,4)==3,
   a=contr(v,1)*contr(v,2);
   b=contr(v,2)*sqrt(1-contr(v,1)^2);
   k=k*(1-2*exp(-a*T(1))*cos(b*T(1))+exp(-2*a*T(1)));
   p=[p;roots([1 -2*exp(-a*T(1))*cos(b*T(1)) exp(-2*a*T(1))])];
   nc=nc+2;
  elseif contr(v,4)==4,
   a=contr(v,1)*contr(v,2);
   b=contr(v,2)*sqrt(1-contr(v,1)^2);
   k=k/(1-2*exp(-a*T(1))*cos(b*T(1))+exp(-2*a*T(1)));
   z=[z;roots([1 -2*exp(-a*T(1))*cos(b*T(1)) exp(-2*a*T(1))])];
   dc=dc+2;
  end
 end

 %t1 = (nc-dc) + sign(contr(3,1))*(nc-dc) + 1;
 t1 = (nc-dc) - contr(3,1);
 if t1>0, z=[z;zeros(t1,1)];
 elseif t1<0, p=[p;zeros(-t1,1)]; end

% if (nc-dc)>0,
%  z=[z;zeros(nc-dc,1)];
% else
%  p=[p;zeros(dc-nc,1)];
% end

end
