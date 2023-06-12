function cp = dfreqcp(T,num,den,w)
% DFREQCP Compute discrete bode of numerator/denominator matrices.
%         CP=DFREQCP(T,NUM,DEN,W) returns the discrete bode of the numerators
%         and denominators in NUM and DEN, respectively. T is the sampling
%         period in seconds and W is the frequency vector in radians/second.
%
%         Rows and columns of NUM and DEN correspond to cases and coefficients,
%         respectively.
%
%         Rows and columns of CP correspond to cases and frequencies,
%         respectively.

% Author: Craig Borghesani
% 4/27/94
% Copyright (c) 2003, Terasoft, Inc.


[rn,cn]=size(num);
[rd,cd]=size(den);
[rw,cw]=size(w);

if rw>1,
 disp('DFREQCP is making w a row vector. Do the same in workspace.');
 w = w(:)';
end

if rn~=rd & (rn~=1 & rd~=1),
 error('If NUM and DEN have unequal rows, one must be only 1 row');
end

i=sqrt(-1);
z=exp(i*w*T);

mx = max(rn,rd);
q=1; r=1;
for h=1:mx,
 q=(rn>1)*h+(rn==1); r=(rd>1)*h+(rd==1);
 upper = polyval(num(q,:),z);
 lower = polyval(den(r,:),z);
 cp(h,:)=upper./lower;
end
