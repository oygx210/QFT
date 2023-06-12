% QFTEX3 Non-parametric uncertainty.

% Author: Y. Chait
% 10/10/94
% 4.18.03 revised for V2 (YC)
% Copyright (c) 2002 by Terasoft, Inc.

clc
clear
echo on

% Example #3 (non-parametric uncertainty) describes the application of
% QFT to a feedback design problem with a non-parametrically uncertain
% plant and several robust performance specifications.

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

%       The plant P(s) has a non-parametric uncertainty model:
%	           |   10                   delta(s) stable                |
%      P(s) = {--------- (1+delta(s)): |delta(jw)|<a(w)               }
%	           |s(0.1s+1)               a(w)=.9(jw/.91+1)/((jw/1.001+1)|

pause % Strike any key to continue
clc

%	The performance specifications are: design a controller
%  G(s) such that it achieves

%  1) Robust stability

%  2) Robust margins (via closed-loop magnitude peaks)
%         | P(jw)G(jw) |
%         |------------| < 1.2,  w>0
%         |1+P(jw)G(jw)|

%  3) Robust sensitivity reduction
%         |     1      |
%         |------------| < 0.089w^2,  w<5
%         |1+P(jw)G(jw)|

pause % Strike any key to continue
echo off
clc

% PROBLEM DATA

%  plant templates
P0 = tf(10,conv([1,0],[.1,1])); % nominal plant
R = tf(.9*[1/.91,1],[1/1.001,1]); % uncertainty model

% BOUNDS
disp(' ')
disp('Computing bounds...')
disp(' ')
disp('bdb1=sisobnds(1,w,W1,P0,R); %stability')
drawnow
W1 = 1.2;  % define weight
w = [0.05,0.5,2,5,90]; % working frequency range
bdb1 = sisobnds(1,w,W1,P0,R);
disp('plotbnds(bdb1); %show bounds')
plotbnds(bdb1),title('Robust Stability Bounds');
qpause;close(gcf);
%
disp(' ')
disp('bdb2=sisobnds(2,wbd2,W2,P0,R); %sensitivity')
drawnow
wbd2 = [0.05,0.5,2,5]; % performance frequency range
W2 = tf(0.089*[1,0,0],1); % performance weight
bdb2 = sisobnds(2,wbd2,W2,P0,R);
disp('plotbnds(bdb2); %show bounds')
drawnow
plotbnds(bdb2),title('Robust Sensitivity Bounds');
qpause;close(gcf);

disp(' ')
disp('bdb=grpbnds(bdb1,bdb2); %grouping bounds')
drawnow
bdb=grpbnds(bdb1,bdb2);
disp('plotbnds(bdb); %show all bounds')
drawnow
plotbnds(bdb),title('All Bounds');
qpause;close(gcf);

disp(' ')
disp('ubdb=sectbnds(bdb); %intersect bounds')
drawnow
ubdb=sectbnds(bdb);
disp('plotbnds(ubdb); %show bounds')
plotbnds(ubdb),title('Intersection of Bounds');
qpause;close(gcf);

% DESIGN
disp(' ')
disp('Design')
disp('lpshape(wl,ubdb,P0,C0); %loop shaping')
drawnow
wl = logspace(-2,8,300);  % define a frequency array for loop shaping
nc0 = [5.032245193701399e+015
       1.047651224010599e+021
       2.783298073926963e+025
       2.091612032833357e+029
       3.180392975934304e+032
       8.265557006790552e+034
       2.627506261665202e+036
       4.929810064633756e+036
       2.260724023982798e+036]';
dc0 = [1.000000000000000e+000
       1.123537042440132e+006
       2.734200460173087e+012
       7.574140703600440e+017
       4.919207207185055e+022
       7.266655789062400e+026
       2.669203151502317e+030
       1.910241138970628e+033
       1.932735079737017e+035
                            0]';
C0 = tf(nc0,dc0);
lpshape(wl,ubdb,P0,C0);  qpause;

% ANALYSIS
disp(' ')
disp('Analysis....')

disp(' ')
disp('chksiso(1,wl(ind),W1,P0,R,C0)); %stability spec')
drawnow
ind = find(wl<=5*10^6);  % performance range
chksiso(1,wl(ind),W1,P0,R,C0);
qpause;close(gcf);
%
disp(' ')
disp('chksiso(2,wl(ind),W2,P0,R,C0); %sensitivity spec')
drawnow
ind = find(wl<=5);  % performance range
chksiso(2,wl(ind),W2,P0,R,C0);
qpause;close(gcf);
