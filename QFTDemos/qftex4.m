% QFTEX4 Classical design for fixed plant.

% Author: Y. Chait
% 2/15/93
% 4.18.03 revised for V2 (YC)
% Copyright (c) 2002 by Terasoft, Inc.

clc
clear
echo on

% Example #4 (classical design) describes the application of QFT to a
% feedback design problem with a fixed plant and both gain margin
% and bandwidth specifications.

% Please refer to manual for more details....

% Strike any key to advance from one plot to another....
pause % Strike any key to continue
clc

%	Consider a continuous-time, siso, negative unity feedback system

%                                   |V(s)            |D(s)
%                         ----      |      ----      |
%             ---->x---->|G(s)|---->x---->|P(s)|---->x---->
%             R(s) |      ----             ----        | Y(s)
%                  |               --                  |
%                  ---------------|-1|------------------
%                                  --

%       The plant model P(s) is fixed:
%               10
%      P(s) = ------
%             s(s+1)

pause % Strike any key to continue
clc

%	The performance specifications are: design a controller
%  G(s) such that it achieves

%  1) Stability

%  2) Gain margin of 1.8

%  3) Zero steady state error for velocity reference commands

%  4) Bandwidth limitation
%         | P(jw)G(jw) |
%         |------------| < 0.707,  w>10
%         |1+P(jw)G(jw)|

pause % Strike any key to continue
echo off
clc

% PROBLEM DATA

% define fixed plant
P0 = tf(10,[1,1,0]);

% BOUNDS
disp(' ')
disp('Computing bounds...')
disp(' ')
disp('bdb1=sisobnds(1,w,W1,P0); %margins')
drawnow
w = [10,100];
W1 = 1.2;  % margin weight
bdb1 = sisobnds(1,w,W1,P0);

disp(' ')
disp('bdb6=sisobnds(6,wbd6,W6,P0); %bandwidth')
drawnow
wbd6 = 10; % performance frequency range
W6 = 0.707; % performance weight
bdb6 = sisobnds(6,wbd6,W6,P0);

disp(' ')
disp('bdb=grpbnds(bdb1,bdb6); %grouping bounds')
drawnow
bdb = grpbnds(bdb1,bdb6);
disp('plotbnds(bdb); %show all bounds')
drawnow
plotbnds(bdb),title('All Bounds');
qpause;close(gcf);

disp(' ')
disp('ubdb=sectbnds(bdb); %intersect bounds')
drawnow
ubdb = sectbnds(bdb);
disp('plotbnds(ubdb); %show bounds')
drawnow
plotbnds(ubdb),title('Intersection of Bounds');
qpause;close(gcf);

% DESIGN
disp(' ')
disp('Design')
disp('lpshape(wl,ubdb,P0,C0); %loop shaping')
drawnow
wl =[logspace(-1,.999,75),logspace(1,2,25)];  % define a frequency array for loop shaping
C0 = tf(0.46*[1,1.4,1],conv([1,0],[1/16^2,1.2/16,1]));
lpshape(wl,ubdb,P0,C0);  qpause;

% ANALYSIS
disp(' ')
disp('Analysis....')

disp(' ')
disp('chksiso(1,wl,W1,P0,R,C0); %margins spec')
drawnow
chksiso(1,wl,W1,P0,0,C0);
qpause;close(gcf);

disp(' ')
disp('chksiso(6,wl(ind),W6,P0,0,C0); %bandwidth spec')
drawnow
ind = find(wl>=10); % performance range
chksiso(6,wl(ind),W6,P0,0,C0);
qpause;close(gcf);
