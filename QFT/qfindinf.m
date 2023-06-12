function [axs,wbds,prob,wbds2] = qfindinf(phase,bdb,flag)
% QFINDINF Read bound info vector. (Utility Function)
%          QFINDINF finds the information that is stored at the end of bound
%          vectors.  It is also used to determine the axis limits required
%          to accomodate all the bounds

% Author: Craig Borghesani
% Date: 9/5/93
% Revised: 2/16/96 10:49 AM V1.1 updates
% Copyright (c) 2003, Terasoft, Inc.


dbrange = [-250, 250];
[r, c]  = size( bdb );

%%%%%% V5 code
% Reason: undefined variable
wbds2 = [];

% find wbds and prob vectors
prob = bdb(  r, : ); wbds = bdb( r-1, : ); wtemp = wbds;
prob = sort( prob ); prob( find(diff(prob)==0) ) = [];
wbds = sort( wbds ); wbds( find(diff(wbds)==0) ) = [];

% determine necessary axis limits
if( length(phase) )
    ymin =  500;
    ymax = -500;
    r    = r-2;
    bdb  = bdb( 1:r, :);
    xmin = min( phase );
    xmax = max( phase );
    for k=1:c
        bdbk = bdb(:, k);
        z = find( bdbk > dbrange(1) & bdbk < dbrange(2) & abs(bdbk)~=248 );
        if( length(z) )
            ymin = min( [ymin; bdbk(z)] );
            ymax = max( [ymax; bdbk(z)] );
        else
            wbds2 = [ wbds2, wtemp(k) ];
        end
    end
    axs = [xmin xmax ymin ymax];
    % if length(wbds2), wbds2=sort(wbds2); wbds2(find(diff(wbds2)==0))=[]; end
    if( nargin==3 )
        ymin    = ymin - 5;
        ymax    = ymax + 5;
        axs     = [xmin xmax ymin ymax];
    end
end
