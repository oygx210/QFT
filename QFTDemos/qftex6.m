% QFTEX6 Missile.

% Yossi Chait
% 2/15/93
% 4.18.03 revised for V2 (YC)
% Copyright (c) 2002 by Terasoft, Inc.

clc
clear
echo on

% Example #6 (missile) describes the application of QFT to a missile
% at different operating points in a flight envelop;  A plant model with
% both parametric and non-parametric uncertainty and different
% performance specifications at different operating points.

% Please refer to manual for more details....

% Strike any key to advance from one plot to another....
pause % Strike any key to continue
clc

%	Consider a continuous-time, siso feedback system

%                                              |V(s)
%                   ----    ----      ----     |      ----
%        ----->x-->|G(s)|--|A(s)|----|M(s)|--->x---->|P(s)|---->x---->
%        R(s)  |    ----    ----      ----            ----     | Y(s)
%              |             -----                             |
%              -------------|-D(s)|----------------------------
%                            -----

%       The plant P(s), simplified dynamics evaluated at three points
%       on the flight envelop, has a mixed parametric and non-parametric
%       uncertainty model:
%	          |                                                         |
%      P(s) = {P (s)(1+delta (s)): delta (s) stable, |delta (jw)|<r (w) }
%	          | i           i           i                  i       i    |
% where
%                           a1*s+a2
%	    P (s) =  ---------------------
%	     i       b1*s^3+b2*s^2+b3*s+b4

%      i=1:  a1=335, a2=237, b1=20.7, b2=39, b3=257, b4=-9.5,  r=0.1
%      i=2:  a1=315, a2=227, b1=19.7, b2=37, b3=247, b4=-9.0,  r=0.05
%      i=3:  a1=345, a2=247, b1=23.7, b2=36, b3=267, b4=-10.5, r=0.075

pause % Strike any key to continue
clc

%       The amplifier A(s) is:
%                     1
%           A(s) = -------
%                  0.01s+1

%       The servo motor M(s) is:
%                            1
%           M(s) = ----------------------
%                  107(0.001s^2++0.13s+1)

%       The rate gyro nD(s) is:
%                     27*40^2*s
%         D(s) = ----------------
%                s^2+1.2*40s+40^2

pause % Strike any key to continue
clc

%	The performance specifications (reflecting different requirements
%  at different operating points) are: design a controller
%  G(s) such that it achieves

%  1) Robust stability

%  2) Robust margins (via closed-loop magnitude peaks)
%         | P(jw)M(jw)A(jw)G(jw)D(jw) |             i=1: W1=1.3
%         |---------------------------| < W1, w>0,  i=2: W1=1.2
%         |1+P(jw)M(jw)A(jw)G(jw)D(jw)|             i=3: W1=1.25

%  3) Robust input disturbance rejection
%         |Y(jw)|                    i=1: W2=0.04
%         |-----| < W2, w in [1,8],  i=2: W2=0.036
%         |V(jw)|                    i=3: W2=0.038

pause % Strike any key to continue
echo off
clc

% PROBLEM DATA

% compute parametric plant templates
a1 = 335; a2 = 237; b1 = 20.7; b2 = 39; b3 = 257; b4 = -9.5;
P(1,1,1) = tf([a1,a2],[b1,b2,b3,b4]);
a1 = 315; a2 = 227; b1 = 19.7; b2 = 37; b3 = 247; b4 = -9.0;
P(1,1,2) = tf([a1,a2],[b1,b2,b3,b4]);
a1 = 345; a2 = 247; b1 = 23.7; b2 = 36; b3 = 267; b4 = -10.5;
P(1,1,3) = tf([a1,a2],[b1,b2,b3,b4]);

%s ervomotor
M = tf(1/107,[0.001,0.13,1]);

% amplifier
A = tf(1,[0.01,1]);

% rate sensor
S = tf(27*[1,0],[1/40^2,1.2/40,1]);

% Compute frequency response
w = [1,3,5,10,30];  % working frequencies

disp(' ')
disp('plottmpl(w,P); % show templates')
drawnow
plottmpl(w,P), title('Plant Templates  (parametric part w/o hardware)')
qpause;close(gcf);

% BOUNDS
disp(' ')
disp('Computing bounds...')
disp(' ')

disp('bdb1=sisobnds(1,w,W1,L); % margins')
drawnow
R = [.1;.05;.075]; % uncertainity disk radius
W1 = [1.3;1.2;1.25]; % margin weight
% include hardware
L = P*M*A*S;
bdb1 = sisobnds(1,w,W1,L);
disp('plotbnds(bdb1); % show bounds')
drawnow
plotbnds(bdb1);
title('Robust Stability Bounds');
qpause;close(gcf);

disp(' ')
disp('bdb3=sisobnds(3,wbd3,W3,P,[],[],H); % input disturbance rejection')
drawnow
wbd3=[1,3,5,10]; % performance frequency range
W3 = [0.04,.036,.038]'; % performance weight
% include hardware response
H = M*S*A;
bdb3 = sisobnds(3,wbd3,W3,P,[],[],H);
disp('plotbnds(bdb3); % show bounds')
drawnow
plotbnds(bdb3);
title('Robust Input Disturbance Rejection Bounds');
qpause;close(gcf);

disp('bdb=grpbnds(bdb1,bdb3); % group bounds')
drawnow
bdb = grpbnds(bdb1,bdb3);
disp('plotbnds(bdb); %s how all bounds')
drawnow
plotbnds(bdb),title('All Bounds');
qpause;close(gcf);

disp(' ')
disp('ubdb=sectbnds(bdb); % intersect bounds')
drawnow
ubdb = sectbnds(bdb);
disp('plotbnds(ubdb); % show bounds')
plotbnds(ubdb);
title('Intersection of Bounds');
qpause;close(gcf);

% DESIGN
disp(' ')
disp('Design')
disp('lpshape(wl,ubdb,L0,C0); % loop shaping')
disp(' ')
drawnow
% define a frequency array for loop shaping
wl = logspace(-3,3.3,100);
% obtain nominal plant from plant set
L0 = L(1,1,1);
% define initial controller
zc = [-2.6176e+001+ 3.2764e+001i
      -2.6176e+001- 3.2764e+001i
      -8.9383e+000+ 1.0620e+001i
      -8.9383e+000- 1.0620e+001i
      -3.7734e+002
      -1.4442e+002
      -7.0175e+001
      -6.1421e+000];
pc = [ 1.6848e+001
      -5.1293e+003+ 5.2329e+003i
      -5.1293e+003- 5.2329e+003i
      -1.2702e+003+ 1.6936e+003i
      -1.2702e+003- 1.6936e+003i
      -1.1330e+003
      -8.1836e+002+ 3.5211e+001i
      -8.1836e+002- 3.5211e+001i
      -3.8289e-001];
kc = 1.7204e+014;
C0 = zpk(zc,pc,kc);
lpshape(wl,ubdb,L0,C0); qpause

% ANALYSIS
disp(' ')
disp('Analysis....')

disp(' ')
disp('chksiso(1,wl,W1,L,[],C0); % margins spec')
drawnow
chksiso(1,wl,W1,L,[],C0);
qpause;close(gcf);

disp(' ')
ind = find((wl>=0.1) & (wl<=10));  % performnce range
disp('chksiso(3,wl(ind),W3,P,[],C0,H); % input disturbance rejection spec')
drawnow
chksiso(3,wl(ind),W3,P,[],C0,H);
qpause;close(gcf);
