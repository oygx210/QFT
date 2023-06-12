function bds = grpbnds(bd1,bd2,bd3,bd4,bd5,bd6,bd7,bd8,bd9)
% GRPBNDS Group QFT bounds.
%         GRPBNDS(BD1,BD2,...,BD9) creates a single bound matrix containing
%         all the passed bounds.
%
%         See also PLOTBNDS, SECTBNDS, SISOBNDS, GENBNDS.

% Author: Craig Borghesani
% 9/6/93, 7/3/03 10:00AM v2.5 updates.
% Copyright (c) 2003, Terasoft, Inc.


bds=[];
for k=1:nargin,
 ch=int2str(k);
 eval(['bds=[bds, bd',ch,'];']);
end
