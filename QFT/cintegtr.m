function h=cintegtr(rt,w,T)
% CINTEGTR Continuous-time integrator. (Utility Function)
%          CINTEGTR computes the frequency response of continuous-time
%          integrators/differentiators or discrete-time delays/predictors.

% Author: Craig Borghesani
% Date: 9/2/93
% Revised: 2/17/96 9:30 AM V1.1 updates
% Copyright (c) 2003, Terasoft, Inc.


% Because all the CAD functions use the same function (ELEMENTS) to
% compute the specific element the user chooses, the following code
% detemines whether it is a continuous or discrete CAD function
% calling.  T=[] in all the continuous CAD environments

%%%%%% V5 change to accomodate nargin change
nargval = nargin;

if nargval==3,
 if T == 0,
  nargval=2;
 end
end

h=ones(1,length(w));

if nargval==2, % continuous integrator(s)/differentiator(s)
 s=sqrt(-1)*w(:)';

% avoid 'Divide by zero' error
 zero=find(s==0);
 if length(zero),
  s(zero)=ones(1,length(zero))*eps;
 end
 h=(1 ./(s.^rt));

else   % discrete delay(s)/predictor(s)
 z=exp(sqrt(-1)*w(:)'*T);
 zero=find(z==0);
 if length(zero),
  z(zero)=ones(1,length(zero))*eps;
 end
 for n=1:abs(rt),
  if rt<0,
   h=h.*z;
  else
   h=h./z;
  end
 end
end
