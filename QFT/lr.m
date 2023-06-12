function [li,ri] = lr(v)
% LR finds the (lower)left-most and (upper)right-most points of the hull
%
% 7.25.96
% yc

[x,i] = min(v(:,1));
[y,j] = min(v(i,2));
if i==1,
 li=1;
else
 li = i(j(1));
end

[x,i] = max(v(:,1));
[y,j] = max(v(i,2));
ri = i(j(1));


