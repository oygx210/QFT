function cp = qcpqft(num,den,w,T)
% QCPQFT Continuous and discrete-time frequency responses. (Utility Function)
%        QCPQFT computes both the continuous and discrete bodes for the
%        IDE environments.

% Author: Craig Borghesani
% Date: 9/5/93
% Revised: 2/17/96 9:33 AM V1.1 updates
% Copyright (c) 2003, Terasoft, Inc.


%%%%%% V5 change to accomodate nargin change
nargval = nargin;

if nargval==4,
 if T == 0,
  nargval=3;
 end
end

[rn,cn]=size(num); [rd,cd]=size(den); [rw,cw]=size(w);
if rw>1,
 disp('QCPQFT is making w a row vector. Do the same in workspace.');
 w = w(:)';
end

i=sqrt(-1);
if nargval==3,
 p=i*w;
else
 p=exp(i*w*T);
end

mx = max(rn,rd);
q=1; r=1;
for h=1:mx,
 q=(rn>1)*h+(rn==1); r=(rd>1)*h+(rd==1);
 upper = polyval(num(q,:),p);
 lower = polyval(den(r,:),p);
 cp(h,:)=upper./lower;
end
