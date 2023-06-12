function [cont,lc]=cntcvrt(con_ar,T)
% CNTCVRT Convert lead/lag and notch terms to pole/zero terms. (Utlilty)
%         CNTCVRT converts Notch and Lead/Lag terms, which are stored in
%         zeta1, zeta2, and wn and phase and w, respectively, and convert
%         them into poles and zeros.

% Author: Craig Borghesani
% Date: 9/10/93
% Revised: 2/18/96 12:08 AM V1.1 updates
% Copyright (c) 2003, Terasoft, Inc.

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

c=1;
cont = [];
loc = find(con_ar(:,4) == 5 | con_ar(:,4) == 6 | con_ar(:,4) == 7);
for k=loc(:)',
   if con_ar(k,4)==5,
      phi=con_ar(k,1); wd=con_ar(k,2);
      if nargval==1,
         alpha=(1-sin(abs(phi*pi/180)))/(1+sin(abs(phi*pi/180)));
         a=wd*sqrt(alpha); b=wd/sqrt(alpha);
         if sign(phi)<0,
            cont=[cont;a,NaN,NaN,1]; cont=[cont;b,NaN,NaN,2];
         else
            cont=[cont;b,NaN,NaN,1]; cont=[cont;a,NaN,NaN,2];
         end

      else
         x=cos(wd*T)+sin(phi*pi/180); y=sin(phi*pi/180)*cos(wd*T)+1;
         a=-(-y+sqrt(y^2-x^2))/x;
         b=-((sin(phi*pi/180)-a)/(-a*sin(phi*pi/180)+1));
         ac=-log(a)/T;
         bc=-log(b)/T;
         cont=[cont;ac,NaN,NaN,2]; cont=[cont;bc,NaN,NaN,1];

      end
      c=c+2;

   elseif con_ar(k,4)==6,
      cont=[cont;con_ar(k,1),con_ar(k,3),NaN,4];
      cont=[cont;con_ar(k,2),con_ar(k,3),NaN,3];

   elseif con_ar(k,4)==7,
      [jk, contnew] = complexlead(con_ar(k,1),con_ar(k,2),con_ar(k,3),[],T);
      if T > 0,
         cont = [cont; contnew(4:end, :)];
      else
         cont = [cont; contnew(3:end, :)];
      end

   end
end
con_ar(loc,:)=[];
cont=[con_ar; cont];
lc=length(cont(:,1));
