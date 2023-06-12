function sub = qsubset(x,y)
% QSUBSET Intersection of two vectors. (Utility Function)
%         QSUBSET determines intersection of two vectors and then returns
%         the locations of these intersections.

% Author: Craig Borghesani
% 9/6/93
% Copyright (c) 2003, Terasoft, Inc.


% making sure to put into row format...just for the hell of it.
x=x(:)'; y=y(:)';

if ~length(x),
    x=y;
end

if length(x)~=length(y),
    for k=1:length(x),
        s=find(x(k)==y);
        if ~length(s), error('Mismatch between w and wbd');
        else sub(k)=s(1); end
    end
else sub=1:length(x); end
