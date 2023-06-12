function rep=repltest
% REPLTEST Replicate test function. (Utility Function)
%          REPLTEST determines which direction Matlab replicates a single
%          number using the : operator.

% Author: Craig Borghesani
% 2/17/94
% Copyright (c) 2003, Terasoft, Inc.



% perform replicating test with one number
x=5;
v=ones(1,5);

% newx should be a row vector with dimensions 1x5
newx = x(v);

[r,c]=size(newx);
rep=1;

% if r > 1 then matlab is replicating incorrectly
if r>1, rep=0; end
