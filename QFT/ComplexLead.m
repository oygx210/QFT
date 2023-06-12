function [h, contnew]=ComplexLead(delph,damp,w,lo,T)
%Finds the element
%               s^2+2*d*w1*s+w1^2
%               -------------------
%               s^2+2*d*w2*s+w2^2
%  whose maximum phase is delph(deg) damping(damp) and occures at w.
% T is the sampling time for a discrete design, if nargin==4 a continous design
% If discrete, the translation from continuous is done by c2dNotch.m
%

% Author: Oded yaniv
% 12/06/03
% Integrated into toolbox 12/07/94 CWB
% Copyright (c) 2003, Terasoft, Inc.


% 1st try 2nd/2nd term
if nargin == 4, T = 0; end
zeta    = damp; %My choice Oded 0.45 is recommended as default
contnew = [];
flag    = 0;
% form num/den of result
% now try simple 2nd order solution
if T ~= 0,
   [w1,w2]   = PhWp2ComplexL2(delph,w,zeta);
   %[nn1,dd1] = c2dNotch(d1,d2,w1,w2,T);
   [nn,dd] = c2dNotch(zeta,zeta,w1,w2,T);
else
   [w1,w2] = PhWp2ComplexL2(delph,w,zeta);
   nn     = [1,2*zeta*w1,w1^2];
   dd     = [1,2*zeta*w2,w2^2];
end
if flag==0,
   contnew = cntpars(tf(nn,dd,T));
   contnew(1,1) = 1;
   h = [];
   if ~isempty(lo),
      h = qcntbode(contnew,lo,T);
   end
end
