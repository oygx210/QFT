% QFTEX15 2x2 MIMO.

% Author: Y. Chait
% 10/6/94
% 4.28.03 revised for V2 (YC)
% Copyright (c) 2002 by Terasoft, Inc.

clc
clear
echo on

% Example #15 (2x2 mimo) describes application of QFT to a robust feedback
% design problem with a 2x2 parametric uncertain plant.  It illustrates
% applicability of various functions to mimo qft design while the toolbox
% is not fully mimo user-friendly

% Please refer to manual for more details....

% Strike any key to advance from one plot to another....
pause % Strike any key to continue
clc

%	Consider a continuous-time, mimo, negative unity feedback system

%                                   |V(s)            |D(s)
%                         ----      |      ----      |
%             ---->x---->|G(s)|---->x---->|P(s)|---->x---->
%             R(s) |      ----             ----        | Y(s)
%                  |               --                  |
%                  ---------------|-1|------------------
%                                  --

%       The 2x2 plant P(s) has parametric a uncertainty model:
%
%	   P(s) = pij(s), i,j=1,2
%          p11(s)=a/d(s); p12(s)=(3+0.5a)/d(s)
%          p21(s)=1/(s);  p22(s)=8/d(s)
%          d(s)=s^2+0.03as+10
%          a in [6,8]
%

pause % Strike any key to continue
clc

%	The performance specifications are: design a diagonal
%  controller G(s)=diag[g1(s),g2(s)] such that it achieves

%  1) Robust stability

%  2) Robust margins (via closed-loop magnitude peaks)
%           |1+Li(jw)| > 1/1.8, i=1,2, w>0
%           (Li is the i'th open-loop function with j'th loop closed)

%  3) Robust sensitivity rejection (S={sij}=(I+PG)^(-1))
%          |sij(jw)| < ai(w), w<10
%          ai = 0.01w,  i=j
%          ai = 0.005w, i not == j

pause % Strike any key to continue
echo off
clc

% PROBLEM DATA

% compute the boundary of the plant templates
disp('Computing outer plant (20 points)....')
aa = linspace(6,8,20);
c = 1;
for a = aa,
  den = [1,0.03*a,10];
  P(1,1,c) = tf(a,den); P(1,2,c) = tf(3+0.5*a,den); 
  P(2,1,c) = tf(1,den); P(2,2,c) = tf(8,den); 
  c = c + 1;
end

disp('Skip showing individual templates...')
disp(' ')

% 1st loop desiign
% BOUNDS
disp(' ')
disp('Starting with design of 1st loop...')
disp(' ')
disp('Computing bounds...')
disp(' ')

% specs:
W11 = 1.8; W22 = 1.8; % stability margins
Ws11 = tf(0.01*[1,0],1);  Ws22 = Ws11; % sensitivity for diagonal terms
Ws12 = tf(0.005*[1,0],1); Ws21 = Ws12; % sensitivity for off-diagonal terms

% margins (1,1): g2=0
disp(' ')
disp('bdb11a=genbnds(10,w,W11,a,b,c,d,P(1,1,1)); % margins (1,1): g2=0');
drawnow
w = [.01,.1,1,5,7.5,10,100];  % working frequencies
a = 1; b = 0; c = 1; d = P(1,1,:);
bdb11a = genbnds(10,w,W11,a,b,c,d,P(1,1,1));

% margins (1,1): g2=inf
disp(' ')
disp('bdb11b=genbnds(10,w,W11,a,b,c,d,P(1,1,1)); % margins (1,1): g2=inf');
detP = P(1,1,:)*P(2,2,:)-P(1,2,:)*P(2,1,:);
a = P(2,2,:); b = 0; c = P(2,2,:); d = detP;
bdb11b = genbnds(10,w,W11,a,b,c,d,P(1,1,1));

% performance (1,1)
disp(' ')
disp('bdbs11 = genbnds(10,wbd,Ws11,a,b,c,d,P(1,1,1)); % s11');
Pw = freqresp(P,w);
a = frd(abs(squeeze(Pw(2,2,:,:))),w)+multmpl(frd(abs(squeeze(Pw(2,1,:,:))),w),frd(abs(squeeze(freqresp(Ws21,w))),w));
wbd = w(1:length(w)-1);
bdbs11 = genbnds(10,wbd,Ws11,a,b,c,d,P(1,1,1));

% performance (1,2)
disp(' ')
disp('bdbs12 = genbnds(10,wbd,Ws12,a,b,c,d,P(1,1,1)); % s12');
a = frd(abs(squeeze(Pw(1,2,:,:))),w)+multmpl(frd(abs(squeeze(Pw(1,2,:,:))),w),frd(abs(squeeze(freqresp(Ws22,w))),w));
bdbs12 = genbnds(10,wbd,Ws12,a,b,c,d,P(1,1,1));

disp(' ')
disp('bdb1=grpbnds(bdb11a,bdb11b,bdbs11,bdbs12); %group bounds')
drawnow
bdb1 = grpbnds(bdb11a,bdb11b,bdbs11,bdbs12);
disp('plotbnds(bdb1); %show all bounds')
plotbnds(bdb1),title('All Bounds: 1st loop');
qpause;close(gcf);

disp(' ')
disp('ubdb1=sectbnds(bdb1); %intersect bounds')
drawnow
ubdb1 = sectbnds(bdb1);
disp('plotbnds(ubdb); %show bounds')
drawnow
plotbnds(ubdb1),title('Intersection of Bounds: 1st loop');
qpause;close(gcf);

% DESIGN
disp(' ')
disp('Design: c1')
disp('lpshape(wl,ubdb1,P(1,1,1),c1); %loop shaping')
drawnow
wl = sort([linspace(3,3.7,25),logspace(-2,3,200)]); % define a frequency array for loop shaping
nc1 = 2.0532e+6*[1,173.4,1.361077e+4,2.75809e+4,1.1270886e+5];
dc1 = [1,854.43,2.5363e+5,2.87192337e+7,7.5729623405e+8,0];
c1 = tf(nc1,dc1);
lpshape(wl,ubdb1,P(1,1,1),c1);  qpause;

disp('Design of 2nd loop....')
% BOUNDS
disp(' ')
disp('Computing bounds...')
disp(' ')

x = inv(1+frd(Pw(1,1,:,:),w)*frd(freqresp(c1,w),w));
temp = multmpl(frd(Pw(1,2,:,:),w),frd(Pw(2,1,:,:),w));
temp = multmpl(temp,c1);
P22e = frd(Pw(2,2,:,:),w) - multmpl(temp,x);

% margins (2,2)
disp(' ')
disp('genbnds(10,w,W22,a,b,c,d,P22e(1,1,1,:)); % margins (2,2)');
a = 1; b = 0; c = 1; d = P22e;
bdb22a = genbnds(10,w,W22,a,b,c,d,P22e(1,1,1,:));

% margins (2,1)
disp(' ')
disp('bdb21a=genbnds(10,w,W11,a,b,c,d,P22e(1,1,1,:)); % margins (2,1)');
a = x; b = multmpl(P(2,2,:),x); c = 1; d = P22e;
bdb21a = genbnds(10,w,W11,a,b,c,d,P22e(1,1,1,:));

% performance (2,1)
disp(' ')
disp('bdbs21=genbnds(10,wbd,Ws21,a,b,c,d,P22e(1,1,1,:)); % s21');
a = multmpl(P(2,1,:),x); b = 0; c = 1; d = P22e;
wbd=w(1:length(w)-1);
bdbs21 = genbnds(10,wbd,Ws21,a,b,c,d,P22e(1,1,1,:));

% performance (2,2)
disp(' ')
disp('bdbs22=genbnds(10,wbd,Ws22,a,b,c,d,P22e(1,1,1,:)); % s22');
a = 1; b = 0; c = 1; d = P22e;
bdbs22 = genbnds(10,wbd,Ws22,a,b,c,d,P22e(1,1,1,:));

disp(' ')
disp('bdb1=grpbnds(bdb22a,bdb21b,bdbs21,bdbs22); %group bounds')
drawnow
bdb2 = grpbnds(bdb22a,bdb21a,bdbs21,bdbs22);
disp('plotbnds(bdb2); %show all bounds')
plotbnds(bdb2),title('All Bounds: 2st loop');
qpause;close(gcf);

disp(' ')
disp('ubdb2=sectbnds(bdb2); %intersect bounds')
drawnow
ubdb2 = sectbnds(bdb2);
disp('plotbnds(ubdb2); %show bounds')
drawnow
plotbnds(ubdb2),title('Intersection of Bounds: 2st loop');
qpause;close(gcf);

% DESIGN
disp(' ')
disp('Design: G2')
disp('lpshape(wl,ubdb2,P22e(1,1,1),c2); %loop shaping')
drawnow
x = inv(1+P(1,1,:)*c1);
P22e = P(2,2,:)-P(1,2,:)*P(2,1,:)*x*c1;
nc2 = 6.086878e+5 *[1,173.4,1.36108e+4,2.75809e+4,1.1270886e+5];
dc2 = [1,718.65,1.785867e+5,1.6925058e+7,4.165129287e+8,0];
c2 = tf(nc2,dc2);
lpshape(wl,ubdb2,P22e(1,1,1),c2); qpause;


% ANALYSIS
disp(' ')
disp('Analysis....')

% redefine wl with less points
wl = logspace(0,3);

C = [c1,0;0,c2];
%Lw = freqresp(P*C,wl);
S = inv(eye(2,2)+P*C);
%Sw = inv(eye(2,2)+frd(Lw,wl));
Sw = freqresp(S,wl);
mSw = max(abs(Sw),[],4);

disp(' ')
disp('Computing margins in 1st and 2nd channels')
drawnow

subplot(211)
loglog(wl,W11*ones(1,length(wl)),'--',wl,squeeze(mSw(1,1,:))); grid; title('1st channel')
set(gca,'xlim',[min(wl),max(wl)]);
subplot(212)
loglog(wl,W22*ones(1,length(wl)),'--',wl,squeeze(mSw(2,2,:))); grid; title('2nd channel')
set(gca,'xlim',[min(wl),max(wl)]);
qpause;close(gcf);

disp(' ')
disp('Performance')
drawnow

mWs11 = abs(squeeze(freqresp(Ws11,wl)));
mWs22 = mWs11;
mWs12 = abs(squeeze(freqresp(Ws12,wl)));
mWs21 = mWs12;

wl = logspace(-2,1);  % performance bandwidth
ind = find(wl<=10);
subplot(221)
loglog(wl(ind),mWs11(ind),'--',wl(ind),squeeze(mSw(1,1,ind))); grid; title('s11')
set(gca,'xlim',[min(wl),max(wl(ind))]);
subplot(222)
loglog(wl(ind),mWs12(ind),'--',wl(ind),squeeze(mSw(1,2,ind))); grid; title('s12')
set(gca,'xlim',[min(wl),max(wl(ind))]);
subplot(223)
loglog(wl(ind),mWs21(ind),'--',wl(ind),squeeze(mSw(2,1,ind))); grid; title('s21')
set(gca,'xlim',[min(wl),max(wl(ind))]);
subplot(224)
loglog(wl(ind),mWs22(ind),'--',wl(ind),squeeze(mSw(2,2,ind))); grid; title('s22')
set(gca,'xlim',[min(wl),max(wl(ind))]);
qpause;close(gcf);
