function [unstable]=chkstab(a1,b1,c1,d1,a2,b2,c2,d2,T)
% CHKSTAB Check stability. (Utility Function)
%         CHKSTAB checks the closed loop stability of the nominal loop.

% Author: Yossi Chait
% 2/17/94
% Copyright (c) 2003, Terasoft, Inc.
%       $Revision: 1.5 $

% set unstable flag to 0
unstable = 0;

% form closed-loop A matrix
e=inv(1+d2*d1);
if isempty(a1) | isempty(b1) | isempty(c1),
   A1 = [];
else
   A1 = a1-b1*e*d2*c1;
end

if isempty(b1) | isempty(c2),
   A2 = [];
else
   A2 = -b1*e*c2;
end

if isempty(c1) | isempty(b2),
   A3 = [];
else
   A3 = b2*c1-b2*d1*e*d2*c1;
end

if isempty(a2) | isempty(b2) | isempty(c2),
   A4 = [];
else
   A4 = a2-b2*d1*e*c2;
end

A=[A1,A2;A3,A4];

eg=eig(A);

if T == 0,
 if any(real(eg) >= 0), unstable = 1; end
else
 if any(abs(eg) >= 1), unstable = 1; end
end
