function bds = sisobnds(ptype,w,W,P,R,nom,C,loc,ph_d)
% SISOBNDS Compute Single-Input/Single-Output QFT bounds.
%          SISOBNDS(PTYPE,W,Ws,P,R,NOM,C,LOC,PHS) computes bounds for the
%          closed-loop configuration designated by PTYPE.  W is the set of
%          frequencies designating bounds to compute.  Ws is the performance
%          specification (scalar or LTI/FRD model), P is the frequency response
%          data of the plant (scalar or LTI/FRD model), R is the disk radius for non-
%          parametric uncertainty (scalar or LTI/FRD model), NOM designates the
%          nominal plant and controller array indices.  LOC specifies location of unknown
%          controller in the loop: 1 for G, 2 for H.  PHS specifies the phase
%          (degrees) grid for bound computations (0 to -360 every -5 deg is default).
%
%          PTYPE specifies the closed-loop configuration.
%          PTYPE=1:  |FPGH/(1+PGH)|<Ws
%          PTYPE=2:  |F/(1+PGH)|<Ws
%          PTYPE=3:  |FP/(1+PGH)|<Ws
%          PTYPE=4:  |FG/(1+PGH)|<Ws
%          PTYPE=5:  |FGH/(1+PGH)|<Ws
%          PTYPE=6:  |FPG/(1+PGH)|<Ws
%          PTYPE=7:  Ws1<|FPG/(1+PGH)|<Ws2
%          PTYPE=8:  |FH/(1+PGH)|<Ws
%          PTYPE=9:  |FH/(1+PGH)|<Ws
%
%          In PTYPE=7, Ws is [Ws1,Ws2]
%
%          See also GENBNDS, GRPBNDS, PLOTBNDS, SECTBNDS.

% Author: Craig Borghesani
% Date: 5/21/94, 4/10/03 2:45PM, 7/3/03 10:07AM
% Copyright (c) 2003, Terasoft, Inc.

if ptype < 1 | ptype > 9,
 error('SISOBNDS only accepts problem types between 1 and 9');
end

% removing wbd
wbd = [];

if nargin==4,
 [w,wbd,W,uP,vP,R,nom,uC,vC,loc,ph_r,info]=bndsdef(w,wbd,W,P,[],[],[],[],[],ptype);
elseif nargin==5,
 [w,wbd,W,uP,vP,R,nom,uC,vC,loc,ph_r,info]=bndsdef(w,wbd,W,P,R,[],[],[],[],ptype);
elseif nargin==6,
 [w,wbd,W,uP,vP,R,nom,uC,vC,loc,ph_r,info]=bndsdef(w,wbd,W,P,R,nom,[],[],[],ptype);
elseif nargin==7,
 [w,wbd,W,uP,vP,R,nom,uC,vC,loc,ph_r,info]=bndsdef(w,wbd,W,P,R,nom,C,[],[],ptype);
elseif nargin==8,
 [w,wbd,W,uP,vP,R,nom,uC,vC,loc,ph_r,info]=bndsdef(w,wbd,W,P,R,nom,C,loc,[],ptype);
elseif nargin==9,
 [w,wbd,W,uP,vP,R,nom,uC,vC,loc,ph_r,info]=bndsdef(w,wbd,W,P,R,nom,C,loc,ph_d,ptype);
else
 error('Improper number of inputs');
end

if nargout == 1,
 eval(['bds=sisobnd',int2str(ptype),'(w,wbd,W,uP,vP,R,nom,uC,vC,loc,ph_r,info);'],'qfterror(1,info)');
else
 eval(['sisobnd',int2str(ptype),'(w,wbd,W,uP,vP,R,nom,uC,vC,loc,ph_r,info);'],'qfterror(1,info)');
end
