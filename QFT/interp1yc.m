function  [A,b] = interp1yc(up,uz,A,b,pr,pc0,pc1);
% sets up the interpolation conditions
% for unstable and nmp OL roots
% assuming simple roots

% 8.26.96
% 8.29.96 modif: allow nmp
% yc

clear i

ipr=length(pr);
ipc0=length(pc0);
[ipc1,tmp]=size(pc1);

% unstable pole
if length(up),
   w=up/i;
   tr=[];  tc1=[];
   if ipr,
      for h=1:ipr,
         tr=[tr,freqcp(1,[1/pr(h),1],w)];
      end
   end
   if ipc1,
      for h=1:ipc1,
         tc1=[tc1,freqcp([1],[1/pc1(h,2)^2,2*pc1(h,1)/pc1(h,2),1],w)];
      end
   end
   t=[tr,tc1*i*w,tc1];
   W=1;
   if ipc0,
      for h=1:ipc0,
         W=W * freqcp(1,[1/pc0(h),1],w);
      end
   end
   t=t * W;

   A=[t;A];
   b=[1;b];
end

% nmp zero
if length(uz),
   tr=[];  tc1=[];
   w=uz/i;
   if ipr,
      for h=1:ipr,
         tr=[tr,freqcp(1,[1/pr(h),1],w)];
      end
   end
   if ipc1,
      for h=1:ipc1,
         tc1=[tc1,freqcp([1],[1/pc1(h,2)^2,2*pc1(h,1)/pc1(h,2),1],w)];
      end
   end
   t=[tr,tc1*i*w,tc1];
   W=1;
   if ipc0,
      for h=1:ipc0,
         W=W * freqcp(1,[1/pc0(h),1],w);
      end
   end
   t=t * W;

   A=[t;A];
   b=[0;b];
end


