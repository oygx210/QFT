function [P]=multmpl(P1,P2,typ)
% MULTMPL Multiplication of two LTI/FRD arrays.
%         MULTMPL(P1,P2,FLAG) produces the user-specified multiplication of P1
%         and P2.  FLAG = 1 for correlated; FLAG = 2 for uncorrelated.
%
%         Correlated implies that the i'th model in P1 is matched
%         with the i'th model in P2.
%         Uncorrelated implies that the i'th model in P1 is matched
%         with all other models in P2.
%
%         See also ADDTMPL, CLTMPL.

% Author: Craig Borghesani
% 8/31/93, 7/3/03 10:02AM : v2.5 updates.
% Copyright (c) 2003, Terasoft, Inc.

if ~isa(P1,'lti'),
   error('MULTMPL only for LTI arrays.');
end

if nargin==2, typ=1; end

rp1=prod(size(P1)); rp2=prod(size(P2));

if rp1~=1 & rp2~=1 & rp1~=rp2 & typ==1,
   disp('You have uncorrelated data. Program setting FLAG=2'); typ=2;
end

rm=max(rp1,rp2);
if typ==1,
   for a=1:rm,
      x=(rp1>1)*a+(rp1==1); y=(rp2>1)*a+(rp2==1);
      P(1,1,a)=P1(1,1,x)*P2(1,1,y);
   end
else
   ab=1;
   for a=1:rp1, for b=1:rp2,
      P(1,1,ab) = P1(1,1,a) * P2(1,1,b);
      ab=ab+1;
   end; end
end

