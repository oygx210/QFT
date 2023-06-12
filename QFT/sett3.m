function [A,b] = sett2(w,pr,pc0,pc1,m,b0,li,ri,ix)
% sets up the LP program
%

% 7.29.96
% 8.15.96   modify: allow complex poles in T expansion
% 8.16.96   modify: LP in terms of T (not Q)
% 8.19.96   modify: add 1/(s+p1)..(s+pn) weight on T to get relative degree in G
% 8.20.96   modify: use 1/(s/pi+1)  terms (dc gain=1) 
% yc

clear i

ipr=length(pr);
ipc0=length(pc0);
[ipc1,tmp]=size(pc1);

A=[]; b=[];

% set up A and b matrices in Ax \leq b
% 1st point is extreme left, so response must lie above

% 1st compute complex response of basis functions
tr=[];  tc1=[];
% real poles
if ipr,
   for h=1:ipr,
      tr=[tr,freqcp(1,[1/pr(h),1],w)];
   end
end
% complex poles with a zero
if ipc1,
   for h=1:ipc1,
      tc1=[tc1,freqcp([1],[1/pc1(h,2)^2,2*pc1(h,1)/pc1(h,2),1],w)];
   end
end
t=[tr,tc1*i*w,tc1];

% pre-weight T
W=1;
if ipc0,
   for h=1:ipc0,
      W=W * freqcp(1,[1/pc0(h),1],w);
   end
end
t=t * W;

rt=real(t);  it=imag(t);  

% vertices from li to ri (ie, response to stay above each)
for n=1:ri-1, 
   A=[A;-it+m(n)*rt];
   b=[b;-b0(n)]; 
end

% vertices from ri to last one (ie, response to stay below each)
for n=ri:ix, 
   A=[A;it-m(n)*rt];
   b=[b;b0(n)]; 
end
