function [bdb,ptype,phase,axs,pos,wbs,wbs2,coora,coorb] = qplotdef(bdb,ptype,...
    phase,pos)
% QPLOTDEF Set defaults for PLOTBNDS. (Utility Function)
%          QPLOTDEF sets all the defaults for variables that were either not
%          specified or passed as [] for PLOTBNDS.

% Author: Craig Borghesani
% 9/6/93
% Copyright (c) 2003, Terasoft, Inc.
%       $Revision: 1.5 $

% load user default values
defs = qftdefs;

[r,c]=size(bdb);
if r==0,
    error('No bound(s) to plot');
end

if any(bdb(r,:)==13)
    flag=2;
else
    flag=1;
end

if ~length(phase)
    phase = defs(1,1):defs(1,2):defs(1,3);
end

if ((r-2)/2)~=length(phase)
    error('Bound information does not match phase vector');
end

if flag<2
    np=length(ptype); p=[];
    if np~=0
        for l=1:np
            p=[p, find(ptype(l)==bdb(r,:))];
        end
    else
        p=1:c;
    end
else
    p=1:c;
end

bdb=bdb(:,p);
[axsd,wbs,ptype,wbs2]=qfindinf(phase,bdb,1);

axs=[axsd(1) axsd(2) axsd(3)-5 axsd(4)+5];

[coora,coorb] = wherebnd(bdb);

if ~length(pos),
    pos=[0.333,0.28,0.6620,0.6604];
end
