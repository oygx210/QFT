function [nc,dc] = autot3c(np0,dp0,w,ic,zeta,iter,optype,index,li,ri,m,b0,ix,bdb)
% convbd(ubdb,[]);
% run first [index,li,ri,m,b0,ix] = convbd(ubdb,phs)
% ic=3;zeta=.6;iter=4;[nc,dc]=autot3d(np0,dp0,w,ic,zeta,iter,optype);
% sets up the QP program
%

% 7.29.96
% 8.15.96  modify:  allow 2nd order basis expansion
% 8.16.96  modify:  LP in terms of T (not Q)
% 8.19.96  modify:  add 1/(s+p1)..(s+pn) weight on T to get relative degree in G
% 8.19.96  modify:  HF gain is T(jw) with w being a large number
% 8.20.96  modify:  use 1/(s/pi+1) terms (dc gain=1)
% 8.29.96  modify:  pre-compute convex hull of bounds
% 8.27.96  modify:  auto determination of location of basis poles
% 8.27.96  modify:  allow T(jw)=0 in case with nmpe OL zeros
% 8.29.96  modify:  compute [nc,dc] using [z,p,k] forms
% 8.30.96  modify:  reduce # of input arguments
% 8.30.96  modify:  allow iterations to minimize high-freq gain
% 9.9.96   modify:  allow for either LP or QP optimization

clear i

[ij,ib]=size(w);

% load pre-computed bounds data
%load bddata

i1=zeros(1,38);
i11=[1.25:.5:10.5];
i1(1,[1:2:37])=i11;
i2=zeros(1,38);
i22=[.95:-.05:.05];
i2(1,[2:2:38])=i22;
ifac=[1,i1+i2];
hfg=1e55;  hown='no';
for g=1:iter,
   fac=ifac(g);
% assign basis poles for optimization
   [pr,pc1,pc0]=basis2b(np0,dp0,bdb,w,ic,zeta,fac);
   ipr=length(pr);
   [ipc1,ij]=size(pc1);
   ipc0=length(pc0);
   ip=ipr+2*ipc1;

% set up A and b matrices in Ax \leq b
   A=[]; b=[];
   for j=1:ib,
      inx=[index(j,1):index(j,2)];
      in=index(j,2)-index(j,1)+1;
      [A1,b1] = sett3(w(j),pr,pc0,pc1,m(inx),b0(inx),li(j),ri(j),in);
      A=[A;A1]; b=[b;b1];
   end
% nmp and unstable OL plant constraints
   rz=roots(np0);
   rp=roots(dp0);
   inmp=find(real(rz) >= 0);
   iuns=find(real(rp) >= 0);
   [A,b] = interp1yc(rp(iuns),rz(inmp),A,b,pr,pc0,pc1);

% now set up minimization of hf gain
   whf=100*max(w);
% n/d of pre-weight
   ntc0=1; dtc0=1;
   if ipc0,
      for h=1:ipc0,
         dtc0=conv(dtc0,[1/pc0(h),1]);
      end
   end
   W=freqcp(ntc0,dtc0,whf);
   [H,c] = setqp2(whf,pr,pc1);
   H=H*abs(W)^2;
% solve the QP
   ieq=0;
   if iuns > 0
      ieq=1;
   end
   if inmp > 0,
      ieq=ieq+1;
   end

   if ieq,
      if optype==1
         c=[zeros(1,ip)];
         [x,temp,how]=qp(H,c,A,b,[],[],[],ieq,-1);
      else
         [x,temp,how]=lp(c,A,b,[],[],[],ieq,-1);
      end
   else
      if optype==1
         c=[zeros(1,ip)];
         [x,temp,how]=qp(H,c,A,b,[],[],[],[],-1);
      else
         [x,temp,how]=lp(c,A,b,[],[],[],[],-1);
      end
   end

% exctract controller
% 1st convert T to a num-den format

% n/d with real part only
   ntr=0; dtr=1;
   if ipr,
      [ntr,dtr]=residue(x(1:ipr).*(pr'),-pr,[]);
   end

% n/d with complex part only with
   ntc1=0; dtc1=1;
   if ipc1,
      for j=1:ipc1,
         [ntc1,dtc1]=addnd(ntc1,dtc1,[0,x(ipr+j),x(ipr+ipc1+j)],[1/pc1(j,2)^2,2*pc1(j,1)/pc1(j,2),1]);
      end
   end

% combine
   [nt1,dt1]=addnd(ntr,dtr,ntc1,dtc1);
   [nt,dt]=mulnd(nt1,dt1,ntc0,dtc0);
   c=dt(1); dt=dt/c; nt=nt/c;
   nt=real(nt); dt=real(dt);

   [z1,p1,k1]=tf2zp(nt,dt-nt);
   [zp0,pp0,kp0]=tf2zp(np0,dp0);
   inmp=find(real(zp0) >= 0);
   iuns=find(real(pp0) >= 0);
   if iuns > 0
      ip=find(abs(p1 - pp0(iuns)) < 1e-8);
      if ip > 0,
         if (pp0(iuns)~=0),
            kp0=kp0*pp0(iuns);
            k1=k1*p1(ip);
         end
         p1(ip)=[];  pp0(iuns)=[];
      end
   end
   if inmp > 0,
      iz=find(abs(z1 - zp0(inmp)) < 1e-8);
      if iz > 0,
         if (zp0(inmp)~=0)&(z1(iz)~=0),
            kp0=kp0/zp0(inmp);
            k1=k1/z1(iz);
         end
      end
      z1(iz)=[];  zp0(inmp)=[];
   end
   zc=[z1;pp0]; pc=[p1;zp0];  kc=k1/kp0;
   hfgn=20*log10(abs(freqcp(nt,dt,whf)));
%   [20*log10(sqrt(abs(x'*H*x))),whf,hfgn]
   if (hfgn<hfg) & strcmp(how,'ok'),
      zcn=zc; pcn=pc; kcn=kc;
      hfg=hfgn; hown=how;
   end
end
if hown=='ok',
%   disp('got a solution for ya...')
%   putqft('junk.shp',[],zcn,pcn,kcn);
   [nc,dc]=zp2tf(zcn,pcn,kcn);
else
%   disp('cannot find a solution ... try higher orders')
   nc=[]; dc=[];
end
