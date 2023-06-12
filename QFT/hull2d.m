function v = hull2d(p)
%function v = hull2d(p)
%Find the convex hull of a set of 2D points.
%p is a np x 2 array; v is a nv x 2 array,
%with the last point being a duplicate of the first if nv>2.
%(That way, plot(v(:,1),v(:,2)) yields a plot of the polygon.)
%
% Author: Richard Schreier, Oregon State University. <schreier@ece.orst.edu>

%Find the (lower)left-most and (upper)right-most points.
[x,i] = min(p(:,1));
[y,j] = min(p(i,2));
l = [x(j(1)) y(1)];
li = i(j(1));

[x,i] = max(p(:,1));
[y,j] = max(p(i,2));
r = [x(j(1)) y(1)];

if size(p,1) <=2        % degenerate cases
    if l == r
        v = l;
    else
        v = [l;r];
    end
    return
end

%Split the points into those above and those below the l-r line.
i = leftof(p,l,r);
above = [p(i,:);l];
below = p(~i,:);

%Sort them in terms of increasing first coordinate.
[junk,i] = sort(above(:,1));
above = above(i,:);
[junk,i] = sort(below(:,1));
below = below(i,:);     % includes the l and r points

%Move along the underside, building the vertex list as we go.
a = below(1,:);
nb = size(below,1);
b = below(2,:);
v = [a;b];
nv = 2;

for i=3:nb
    p = below(i,:);
    while( ~leftof(p,a,b) ) % backtrack vertices
        nv = nv-1;
        v = v(1:nv,:);
        b = a;
        if nv>1
            a = v(nv-1,:);
        else
            break;
        end
    end
    v = [v;p];
    nv = nv+1;
    a = b;
    b = p;
end

%Move along the top side, continuing to build the vertex list.
na = size(above,1);
for i=na:-1:1
    p = above(i,:);
    while( ~leftof(p,a,b) ) % backtrack vertices
        nv = nv-1;
        v = v(1:nv,:);
        b = a;
        a = v(nv-1,:);
    end
    v = [v;p];
    nv = nv+1;
    a = b;
    b = p;
end
