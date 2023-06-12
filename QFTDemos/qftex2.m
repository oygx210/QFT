% QFTEX2 QFT Tracking.

% Author: Y. Chait
% 2/15/93
% 4.18.03 revised for V2 (YC)
% Copyright (c) 2002 by Terasoft, Inc.
%       $Revision: 1.6 $

clc
clear
echo on

% Example #2 (2 DOF) describes the application of QFT to a 2
% degree-of-freedom feedback design problem with a parametrically
% uncertain plant and QFT tracking specification.

% Please refer to manual for more details....

% Strike any key to advance from one plot to another....
pause % Strike any key to continue
clc

%	Consider a continuous-time, siso, negative unity feedback system

%                                       |V(s)            |D(s)
%            ----             ----      |      ----      |
%      ----->|F(s)|--->x---->|G(s)|---->x---->|P(s)|---->x---->
%      R(s)  ----      |      ----             ----      | Y(s)
%                      |             --                  |
%                      -------------|-1|------------------
%                                    --

%       The plant P(s) has a parametric uncertainty model:
%	            |    k                                |
%	    P(s)  = { ------- :  k in [1,10], a in [1,10] }
%	            | s (s+a)                             |

pause % Strike any key to continue
clc

%	The performance specifications are: design a controller
%  G(s) and a pre-filter F(s) such that they achieve

%  1) Robust stability

%  2) Robust margins (via closed-loop magnitude peaks)
%         | P(jw)G(jw) |
%         |------------| < 1.2,  w>0
%         |1+P(jw)G(jw)|

%  3) Robust tracking (related to tracking step responses)
%                |       P(jw)G(jw) |
%         a(w) < |F(jw) ------------| < b(w),  w<10
%                |      1+P(jw)G(jw)|

%                |           120              |
%         a(w) = |----------------------------|
%                |(jw)^3+17(jw)^2+828(jw)+120 |

%                |  0.6584(jw+30)    |
%         b(w) = |-------------------|
%                |(jw)^2+4(jw)+19.752|

pause % Strike any key to continue
echo off
clc

% PROBLEM DATA

% computing the boundary of plant templates
disp('Computing plant templates (14 plants at 6 frequencies)....')
drawnow
c = 1;	a = 1;
for k = [1,2.5,6,10],
   P(1,1,c) = tf(k*a,[1,a,0]);    c = c + 1;
end
a = 10;
for k = [1,2.5,6,10],
   P(1,1,c) = tf(k*a,[1,a,0]);    c = c + 1;
end
k = 1;
for a = [1.5,3,6],
   P(1,1,c) = tf(k*a,[1,a,0]);    c = c + 1;
end
k = 10;
for a = [1.5,3,6],
   P(1,1,c) = tf(k*a,[1,a,0]);    c = c + 1;
end

% Define working frequencies
w = [.1,.5,1,2,10,15,100];

disp(' ')
disp('plottmpl(w,P); %show templates')
drawnow
plottmpl(w,P), title('Plant Templates')
qpause;close(gcf);

% BOUNDS
disp(' ')
disp('Computing bounds...')
disp(' ')
disp('bdb1=sisobnds(1,w,W1,P,R); %margins');
drawnow
W1 = 1.2; % robust margins weight
bdb1 = sisobnds(1,w,W1,P);
disp('plotbnds(bdb1); %show bounds');
plotbnds(bdb1),title('Robust Margins Bounds');
qpause;close(gcf);
%
disp(' ')
disp('bdb7=sisobnds(7,wbd7,W7,P); %tracking bounds');
drawnow
wbd7 = w(1:5);
mu = tf(0.6584*[1,30],[1,4,19.752]);
ml = tf(120,[1,17,82,120]);
W7 = [mu;ml]; % tracking weight
bdb7 = sisobnds(7,wbd7,W7,P);
disp('plotbnds(bdb7); %show bounds');
plotbnds(bdb7),title('Robust Tracking Bounds');
qpause;close(gcf);

disp(' ')
disp('bdb=grpbnds(bdb1,bdb7); %grouping bounds')
drawnow
bdb=grpbnds(bdb1,bdb7);
disp('plotbnds(bdb); %show all bounds')
plotbnds(bdb),title('Margins and Tracking Bounds');
qpause;close(gcf);
%
disp(' ')
disp('ubdb=sectbnds(bdb); %intersect bounds')
drawnow
ubdb=sectbnds(bdb);
disp('plotbnds(ubdb); %show bounds')
drawnow
plotbnds(ubdb),title('Intersection of Bounds');
qpause;close(gcf);

% DESIGN
disp(' ')
disp('Design')
disp('lpshape(wl,ubdb,P(1,1,nompt),G); %loop shaping')
drawnow
nompt=1;
%nL0=nump(nompt,:); dL0=denp(nompt,:); % nominal plant
del0=0;
wl=logspace(-2,4,200);
nc0 = [3.0787e+6,3.5365e+8,3.8529e+8];
dc0 = [1.0,1.5288e+3,1.0636e+6,4.2810e+7];
G = tf(nc0, dc0);
lpshape(wl,ubdb,P(1,1,nompt),G);  qpause;

disp(' ')
disp('pfshape(7,wl(ind),W7,P,[],G,[],F); %pre-filter shaping')
drawnow
nf0=1; df0=[1/4^2,2*0.7/4,1];
F = tf(nf0, df0);
ind=find(wl<=10); % tracking performance up to w=15
pfshape(7,wl(ind),W7,P,[],G,[],F); qpause;

% ANALYSIS
disp(' ')
disp('Analysis....')
disp(' ')
disp('chksiso(1,wl,W1,P,[],G); %margins spec')
drawnow
chksiso(1,wl,W1,P,[],G);
qpause;close(gcf);
%
disp(' ');
disp('chksiso(7,wl(ind),W7,P,[],G,[],F); %tracking spec')
drawnow
chksiso(7,wl(ind),W7,P,[],G,[],F);
qpause;close(gcf);
