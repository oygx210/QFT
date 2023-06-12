function ph=qatan4(cp)
% QATAN4 Four-quadrant arctangent. (Utility Function)
%        QATAN4 provides a 4 quadrant arctangent of a complex number.
%        Starting at 0 and going through -360 degrees.

% Author: Craig Borghesani
% 8/31/93
% Copyright (c) 2003, Terasoft, Inc.


ph = atan2(imag(cp),real(cp));
ph = ph - 2*pi*(ph>1e-5);
