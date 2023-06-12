function cp = freqcp(num,den,w,delay)
% FREQCP Compute continuous bode of numerator/denominator matrices.
%        CP=FREQCP(NUM,DEN,W,DELAY) returns the continuous bode of the
%        numerators and denominators in NUM and DEN, respectively.  W is
%        the frequency vector in radians/second.  DELAY is the pure
%        time delay and it must have as many rows as NUM or DEN or be
%        just one element.
%
%        CP=FREQCP(NUM,DEN,W) returns the continuous bode as above,
%        except DELAY defaults to zero (0).
%
%        Rows and columns of NUM and DEN correspond to cases and
%        coefficients, respectively.
%
%        Rows and columns of CP correspond to cases and frequencies,
%        respectively.

% Author: Craig Borghesani
% 4/27/94
% Copyright (c) 2003, Terasoft, Inc.


if nargin==3, delay=0; end

[rn,cn]=size(num);
[rd,cd]=size(den);
[rw,cw]=size(w);
[rdel,cdel]=size(delay);

if rw>1,
 disp('FREQCP is making w a row vector. Do the same in workspace.');
 w = w(:)';
end

if rn~=rd & (rn~=1 & rd~=1),
 error('If NUM and DEN have unequal rows, one must be only 1 row');
end

if rdel~=rn & rdel~=rd & rdel~=1,
 error('DELAY must have as many rows as either NUM or DEN or be only 1 element');
end

i=sqrt(-1);
s=i*w;

mx = max(rn,rd);
q=1; r=1;
for h=1:mx,
 q=(rn>1)*h+(rn==1); r=(rd>1)*h+(rd==1); d=(rdel>1)*h+(rdel==1);
 upper = polyval(num(q,:),s);
 lower = polyval(den(r,:),s);
 cp(h,:)=(upper./lower).*exp(-s*delay(d));
end
