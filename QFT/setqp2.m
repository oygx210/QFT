function [H,c] = setqp2(w,pr,pc1)
% sets up the QP factor for minimization of very HF gain


% 8.20.96
% yc

clear i

ipr=length(pr);
[ipc1,tmp]=size(pc1);
ip=ipr+2*ipc1;

tr=[];  tc1=[];

% minimize gain of T at whf=w

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

% NOTE: no need to include pre-weight
rt=real(t);  it=imag(t);

% For QP
RT=rt'*rt;   IT=it'*it;  
H=diag(rt.^2+it.^2);
for g=1:ip-1,
   H(g,(g+1):ip)=2*(RT(g,(g+1):ip)+IT(g,(g+1):ip));
end

% for LP
c=abs(rt)+abs(it);


