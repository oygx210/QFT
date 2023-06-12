function [m,b0] = mb(v)
% mb computed the line parameters for the vertices wth
% the points in v being the edges of the convex hull
%

% 7.25.96
% yc

[ir,ix]=size(v);

% compute slope and constant
for k=1:ir-1,
 m(k,1)=(v(k+1,2)-v(k,2))/(v(k+1,1)-v(k,1));
 b0(k,1)=v(k,2)-m(k)*v(k,1);
end

% last point is connected to first point
m(ir,1)=(v(ir,2)-v(1,2))/(v(ir,1)-v(1,1));
b0(ir,1)=v(ir,2)-m(ir)*v(ir,1);
