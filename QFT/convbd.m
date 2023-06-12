function [index,li,ri,m,b0,ix] = convbd(bdb,phs)
% computes bound on T from those on L
% using convex hull in T space

% 8.26.96

clear i

[ij,ib]=size(bdb);
w=bdb(ij-1,:);

index=[]; li=[]; ri=[]; m=[]; b0=[];  in=0;  ix=[];

for j=1:ib,
% extract bound
   if length(phs),
      bdbo = cpparse2(bdb(:,j),phs);
   else
      bdbo = cpparse2(bdb(:,j));
   end

% compute convex hull
   bdbcvx = hull2d(bdbo);
% remove last point (duplicate of 1st point)
   ix1=length(bdbcvx);
   bdbcvx=bdbcvx(1:ix1-1,:);
   ix1=length(bdbcvx);

% find the (lower)left-most and (upper)right-most points.
   [li1,ri1] = lr(bdbcvx);
% compute line eqns for each vertex
   [m1,b01] = mb(bdbcvx);

% allocate to matrices
   index=[index;[1,ix1]+in]; in=in+ix1; 
   li=[li;li1];
   ri=[ri;ri1];
   m=[m;m1];
   b0=[b0;b01];
end

save bddata index li ri m b0 ix bdb w
