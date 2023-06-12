function [P,G,H,F,sgn,typ]=clcpdef(P,G,H,F,sgn,typ)
% CLCPDEF Defaults for CLCP. (Utility Function)
%         CLCPDEF sets the defaults for whatever the user either passed in
%         as [] or did not specify at all for CLCP.

% Author: Craig Borghesani
% Date: 9/2/93
% Revised: 2/17/96 11:03 AM V1.1 updates
% Copyright (c) 2003, Terasoft, Inc.

if ~length(P), P=1; end
if ~length(G), G=1; end
if ~length(H), H=1; end
if ~length(F), F=1; end
if ~length(sgn), sgn=-1; end
if ~length(typ), typ=1; end

if ~isa(P, 'lti'),
   [rp,cp]=size(P); [rg,cg]=size(G); [rh,ch]=size(H); [rf,cf]=size(F);
   chk=[cp cg ch cf]; chk(find(chk==1))=[];
   if diff(chk)~=0 & length(chk)>1,
    error('Elements in clcp must have equal columns or be just one number');
   end

   % make all matrices that have only one element the same size of the other
   % matrices

   cm=max([cp cg ch cf]);

   if cp == 1, P=P(:,ones(1,cm)); end
   if cg == 1, G=G(:,ones(1,cm)); end
   if ch == 1, H=H(:,ones(1,cm)); end
   if cf == 1, F=F(:,ones(1,cm)); end
end;%if

