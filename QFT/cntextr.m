function [num,den] = cntextr(cont,T)
% CNTEXTR Controller extraction. (Utility Function)
%         CNTEXTR extracts the controller out of the controller matrix
%         created within the IDE and places it into num/den format.

% Author: Craig Borghesani
% 9/3/93
% Date: 7/1/03 2:09PM : Updated gain modified by discrete
%                       integ/diff.
% Copyright (c) 2003, Terasoft, Inc.
%      $Revision: 1.7 $

% Because all the CAD functions use this function, the following code
% detemines whether it is a continuous or discrete CAD function
% calling.  T=[] in all the continuous CAD environments

nargval = nargin;

if nargval==2,
 if T == 0,
  nargval=1;
 end
else
   T = 0;
end

% Convert Notch, Lead/Lag, and Complex Lead/Lag terms into pole and zero terms
[cont,lc]=cntcvrt(cont,T);

num=1;nc=0;den=1;dc=0;
if nargval==1,
 for p=2:lc,
  if cont(p,4)==1,
   den=conv(den,[1 cont(p,1)]);
   cont(1,1)=cont(1,1)*cont(p,1);
  elseif cont(p,4)==2,
   num=conv(num,[1 cont(p,1)]);
   cont(1,1)=cont(1,1)/cont(p,1);
  elseif cont(p,4)==3,
   den=conv(den,[1 2*cont(p,1)*cont(p,2) cont(p,2)^2]);
   if cont(p,2)~=0, cont(1,1)=cont(1,1)*cont(p,2)^2; end
  elseif cont(p,4)==4,
   num=conv(num,[1 2*cont(p,1)*cont(p,2) cont(p,2)^2]);
   if cont(p,2)~=0, cont(1,1)=cont(1,1)/cont(p,2)^2; end
  end
 end
else
 for p=2:lc,
  if cont(p,4)==1,
   pz=real(exp(-cont(p,1)*T));
   den=conv(den,[1 -pz]);
   cont(1,1)=cont(1,1)*(1-pz);
   nc=nc+1;

  elseif cont(p,4)==2,
   zz=real(exp(-cont(p,1)*T));
   num=conv(num,[1 -zz]);
   cont(1,1)=cont(1,1)/(1-zz);
   dc=dc+1;

  elseif cont(p,4)==3,
   a=cont(p,1)*cont(p,2);
   b=cont(p,2)*sqrt(1-cont(p,1)^2);
   den=conv(den,[1 -2*exp(-a*T(1))*cos(b*T(1)) exp(-2*a*T(1))]);
   cont(1,1)=cont(1,1)*(1-2*exp(-a*T(1))*cos(b*T(1))+exp(-2*a*T(1)));
   nc=nc+2;

  elseif cont(p,4)==4,
   a=cont(p,1)*cont(p,2);
   b=cont(p,2)*sqrt(1-cont(p,1)^2);
   num=conv(num,[1 -2*exp(-a*T(1))*cos(b*T(1)) exp(-2*a*T(1))]);
   cont(1,1)=cont(1,1)/(1-2*exp(-a*T(1))*cos(b*T(1))+exp(-2*a*T(1)));
   dc=dc+2;

  elseif cont(p,4)==0.5,

   if cont(p,1)==1,
      den=conv(den,[1 -1]);
      nc=nc+1;

   elseif cont(p,1)==2,
      den=conv(den,[1 -2 1]);
      nc=nc+1;

   elseif cont(p,1)==3,
      den=conv(den,[1 -3 3 -1]); num=conv(num,[1 1]);
      nc=nc+1;

   end
   cont(1,1) = cont(1,1) * T^(cont(p,1));

   if cont(p,2)==1,
      num=conv(num,[1 -1]);
      dc=dc+1;

   elseif cont(p,2)==2,
      num=conv(num,[1 -2 1]);
      dc=dc+1;

   elseif cont(p,2)==3,
      num=conv(num,[1 -3 3 -1]); den=conv(den,[1 1]);
      dc=dc+1;

   end
   cont(1,1) = cont(1,1) / T^(cont(p,2));

  end
 end
end

% take care of integrators/differentiators/delays/predictors
if nargval==2,

% add the appropriate number of preds/delays to make things
% correct for the external world

 if length(num) > 1 | length(den) > 1,
%   t1 = (nc-dc) + sign(cont(3,1))*(nc-dc) + 1;
   t1 = (nc-dc) - (cont(3,1)-cont(3,2));
   if t1>0, num=conv(num,[1 zeros(1,t1)]);
   elseif t1<0, den=conv(den,[1 zeros(1,-t1)]); end
 end

else
 if (cont(2,1)-cont(2,2))>0, den=conv(den,[1 zeros(1,(cont(2,1)-cont(2,2)))]);
 else num=conv(num,[1 zeros(1,abs((cont(2,1)-cont(2,2))))]); end
end

% which ever is of the greatest order remains monotonic
if length(num)<=length(den),
 num=num*cont(1,1);
else
 den=den/cont(1,1);
end
