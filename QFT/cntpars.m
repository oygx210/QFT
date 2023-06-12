function cont = cntpars(C0, cont_t)
% CNTPARS Parse controller. (Utility Function)
%         CNTPARS takes a transfer function and determine the D.C. gain
%         (K), real and complex poles, and real and complex zeros. This
%         information is then stored in what is called the controller
%         matrix. The controller matrix is used by all the IDEs.

% Author: Craig Borghesani
% Date: 9/3/93
% Revised: 2/17/96 9:31 AM V1.1 updates
%          4/16/03 11:30AM LTI support only.
% Copyright (c) 2003, Terasoft, Inc.

% Because all the CAD functions use this function, the following code
% detemines whether it is a continuous or discrete CAD function
% calling.  T=0 in all the continuous CAD environments

%%%%%% V5 change to accomodate nargin change
nargval = nargin;

T = C0.Ts;

% seperate LTI object into zeros, poles, and gain.
[z, p, k] = zpkdata(C0, 'v');
[num, den] = tfdata(C0, 'v');

% obtain D.C. gain (K) of controller
cont(1,1:4)=[NaN NaN NaN 0];
if T == 0,
   n=find(num~=0); d=find(den~=0);
   cont(1,1:4)=[num(n(length(n)))/den(d(length(d))) NaN NaN 0];
% integrators/differentiators
   cont(2,1:4)=[0 0 NaN 0.7];
   c=2;
else
   cont(1,1:4)=[k NaN NaN 0];
% integrators/differentiators
   cont(2,1:4)=[0 0 NaN 0.5];
% delay/predictors
   cont(3,1:4)=[0 0 NaN 0.6];
   c=3;
end

% separate and store poles and zeros
q=1;
delays=0;
preds=0;
nsum=1;
dsum=1;

while q<=length(z),
   c=c+1;
   if imag(z(q))==0,
      if T > 0 & abs(z(q)-1)>1e-10,
         nsum=nsum*(1-z(q));
         if z(q)~=0,
            z(q)=log(z(q))/T;
         end
      end
      cont(c,1:4)=[-z(q) NaN NaN 2];
      q=q+1;
   else
      if T > 0,
         z(q)=log(z(q))/T;
      end
      re=real(z(q)); im=imag(z(q));
      wn=sqrt(re^2+im^2); zta=-re/wn;
      cont(c,1:4)=[zta wn NaN 4];
      q=q+2;
      if T > 0,
         a=zta*wn; b=wn*sqrt(1-zta^2);
         nsum=nsum*(1-2*exp(-a*T)*cos(b*T)+exp(-2*a*T));
      end
   end
end

q=1;
while q<=length(p),
   c=c+1;
   if imag(p(q))==0,
      if T > 0 & abs(p(q)-1)>1e-10,
         dsum=dsum*(1-p(q));
         if p(q)~=0,
            p(q)=log(p(q))/T;
         end
      end
      cont(c,1:4)=[-p(q) NaN NaN 1];
      q=q+1;
   else
      if T > 0,
         p(q)=log(p(q))/T;
      end
      re=real(p(q)); im=imag(p(q));
      wn=sqrt(re^2+im^2); zta=-re/wn;
      cont(c,1:4)=[zta wn NaN 3];
      q=q+2;
      if T > 0,
         a=zta*wn; b=wn*sqrt(1-zta^2);
         dsum=dsum*(1-2*exp(-a*T)*cos(b*T)+exp(-2*a*T));
      end
   end
end

% if discrete, determine the number of integrators
if T > 0,
   y=find(abs(1+cont(:,1))<1e-10 & cont(:,4)==1); cont(y,:)=[];
   x=find(cont(:,1)==0 & cont(:,4)==2);
   yl=length(y);
   if yl,
      if length(x),
         cont(x(1),:)=[];
      else
         delays=delays+1;
      end
      if yl==3,
         x=find(-cont(:,1)==-1 & cont(:,4)==2);
         if length(x),
            cont(x,:)=[];
         else
            cont=[cont; -1 NaN NaN 1];
         end
      elseif yl>3,
         error('Cannot have more than 3 integrators');
      end
   end
   yy=find(abs(1+cont(:,1))<1e-10 & cont(:,4)==2); cont(yy,:)=[];
   x=find(cont(:,1)==0 & cont(:,4)==1);
   yyl=length(yy);
   if yyl,
      if length(x),
         cont(x(1),:)=[];
      else
         preds=preds+1;
      end
      if yyl==3,
         x=find(cont(:,1)==-1 & cont(:,4)==1);
         if length(x),
            cont(x,:)=[];
         else
            cont=[cont; -1 NaN NaN 2];
         end
      elseif yyl>3,
         error('Cannot have more than 3 differentiators');
      end
   end
   yd=find(cont(:,1)~=0 & cont(:,4)==1);
   x=find(cont(:,1)==0 & cont(:,4)==2);
   if length(x)>=length(yd),
      cont(x(1:length(yd)),:)=[];

   elseif length(x)<length(yd),
      delays=length(yd)-length(x)+delays;
      cont(x,:)=[];
   end
   yn=find(cont(:,1)~=0 & cont(:,4)==2);
   x=find(cont(:,1)==0 & cont(:,4)==1);
   if length(x)>=length(yn),
      cont(x(1:length(yn)),:)=[];
   elseif length(x)<length(yn),
      preds=length(yn)-length(x)+preds;
      cont(x,:)=[];
   end

% when our parser parses a controller for our
% second order term, it will see 2 terms at zero
% for each second order term

% find all complex pole elements
   yd=find(cont(:,2)~=NaN & cont(:,4)==3);

% find all z/1 terms
   x=find(cont(:,1)==0 & cont(:,4)==2);

% if we have more z/1 terms than we have
% complex poles, get rid of the correct
% amount
   if length(x)>length(yd),
      cont(x(1:length(yd)*2),:)=[];

% if we have less, meaning that someone did not pass in
% a second order in the required format, make sure
% and compensate for this by adding the correct number of
% delays
   elseif length(x)<=length(yd),
      delays=length(yd)*2-length(x)+delays;
      cont(x,:)=[];
   end

% find all complex zero elements
   yn=find(cont(:,2)~=NaN & cont(:,4)==4);

% find all 1/z terms
   x=find(cont(:,1)==0 & cont(:,4)==1);

% if we have more 1/z terms than we have
% complex zeros, get rid of the correct
% amount
   if length(x)>length(yn),
      cont(x(1:length(yn)*2),:)=[];

% if we have less, meaning that someone did not pass in
% a second order in the required format, make sure
% and compensate for this by adding the correct number of
% predictors
   elseif length(x)<=length(yn),
      preds=length(yn)*2-length(x)+preds;
      cont(x,:)=[];
   end
   x=find(cont(:,1)==0 & cont(:,4)==2);
   xxl=length(x)+preds;
   cont(x,:)=[];
   x=find(cont(:,1)==0 & cont(:,4)==1);
   cont(x,:)=[];
   xl=length(x)+delays;
   cont(2,1)=yl;
   cont(2,2)=yyl;
   cont(3,1)=xxl;
   cont(3,2)=xl;
   if cont(3,1) == cont(3,2),
      cont(3,1:2) = [0,0];
   elseif (cont(3,1) - cont(3,2)) >= 1,
      cont(3,1) = cont(3,1) - cont(3,2);
      cont(3,2) = 0;
   else
      cont(3,2) = cont(3,2) - cont(3,1);
      cont(3,1) = 0;
   end

   n=find(num~=0);
   d=find(den~=0);
   cont(1,1)=nsum/dsum*num(n(1))/den(d(1));

% modify gain based upon number of integ/diff
   cont(1,1) = cont(1,1) / T^(cont(2,1));
   cont(1,1) = cont(1,1) * T^(cont(2,2));

else  % Continuous
 xl=find(cont(:,1)==0 & cont(:,4)==2); cont(xl,:)=[];
 xll=find(cont(:,1)==0 & cont(:,4)==1); cont(xll,:)=[];
 cont(2,1)=length(xll);
 cont(2,2)=length(xl);
end

% merge two controller matrices
if nargin == 2,

%   cont(1,1) = cont(1,1) * cont_t(1,1);
   cont(1,1) = cont_t(1,1);
   if T > 0,
      cont(2, 1:3) = cont(2, 1:3) + cont_t(2, 1:3);
      cont(3, 1:3) = cont(3, 1:3) + cont_t(3, 1:3);
      cont = [cont; cont_t(4:end, :)];

   else
      cont(2, 1:3) = cont(2, 1:3) + cont_t(2, 1:3);
      cont = [cont; cont_t(3:end, :)];

   end

   if T > 0,
      if cont(3,1) == cont(3,2),
         cont(3,1:2) = [0,0];
      elseif (cont(3,1) - cont(3,2)) >= 1,
         cont(3,1) = cont(3,1) - cont(3,2);
         cont(3,2) = 0;
      else
         cont(3,2) = cont(3,2) - cont(3,1);
         cont(3,1) = 0;
      end
   end

end
