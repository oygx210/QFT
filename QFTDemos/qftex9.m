% QFTEX9 Uncertain flexible mechanism.

% Author: Y. Chait
% 2/15/93
% 4.19.03 revised for V2 (YC)
% Copyright (c) 2002 by Terasoft, Inc.


clc
clear
echo on

% Example #9 (flexible servomechanism) describes the application of QFT
% to a feedback design problem with a parametrically uncertain flexible
% mechanical system and several robust performance specifications.

% Please refer to manual for more details....

% Strike any key to advance from one plot to another....
pause % Strike any key to continue
clc

%	Consider a continuous-time, siso feedback system

%          U(s) --             --  V(s)
%          <---|Wu|<--|   |<--|Wv|<---
%               --    |   |    --
%                     |   |           ----       ----
%          ----->x--->x-->x----------|G(s)|----->|P(s)|---->x---->
%          R(s)  |                    ----       ----       | Y(s)
%                |             --                           |
%                -------------|-1|--------------------------x<--|Wu|<--
%                              --                                    N(s)

%       The plant P(s) has a parametric uncertainty model:
%	           |         b1*s+b2           a1,a2,a3,a4,b1 and b2   |
%	   P(s)  = { --------------------- :   functions of c          }
%	           | a1*s^3+a2*s^2+a3*s+a4     c in [0.0111,0.0195] N/m|

pause % Strike any key to continue
clc

%	The performance specifications are: design a controller
%  G(s) such that it achieves

%  1) Robust stability

%  2) Robust margins (via closed-loop magnitude peaks)
%         | P(jw)G(jw) |
%         |------------| < 1.1
%         |1+P(jw)G(jw)|

%  3) Robust tracking
%         |     R(jw)-Y(jw)|
%         |W(w) -----------| < 1
%         |        R(jw)   |

%  4) Robust disturbance rejection
%         |     R(jw(-Y(jw)|
%         |W(w) -----------| < 1
%         |        V(jw)   |


pause % Strike any key to continue
clc

%  5) Robust noise rejection
%         |     R(jw)-Y(jw)|          |     U(jw)|
%         |W(w) -----------| < 1,     |W(w) -----| < 1
%         |        N(jw)   |          |     N(jw)|

%  6) Robust control effort limitation
%         |U(jw)|        |U(jw)|
%         |-----| < 1,   |-----| < 1
%         |R(jw)|        |V(jw)|

%  7) <100hz control terms bandwidth (DSP board limitation)

%       where
%                |2*pi*wb|
%         W(w) = |-------|, wb=10hz (bandwidth),
%                |  jw+1 |

%         wn = 0.01,   wu = 0.33,    wv = 0.1

pause % Strike any key to continue
echo off
clc

% PROBLEM DATA

% computing plant templates
disp('Computing plant templates (22 plants at 9 frequencies)....')
drawnow
% system constants
j1 = 1.45e-6; j2 = j1; km = 34.5e-2; ds = 5e-6; dm1 = .45e-6; dm2 = dm1;
kgen = .05; kt = .1333/10; ka = 20;
% discretize uncertain parameter
cs = (1.11:.04:1.95)/100;
j=1;
for k=cs,
 nump = (kgen*kt*ka)*km*[ds k];
 denp = [j1*j1,j1*(ds+dm1)+j2*(ds+dm1),k*(j1+j2)+dm1*dm2,k*(dm1+dm2)];
 P(1,1,j) = tf(nump,denp); j=j+1;
end


w = [.1,1,10,50,140,155,164,180,500];  % working frequencies
nom = 22;

disp(' ')
disp('plottmpl(w,P,nom); %show templates')
drawnow
plottmpl(w,P,nom),title('Plant Templates')
qpause;close(gcf);

% BOUNDS
disp(' ')
disp('Computing bounds...')
disp(' ')
disp('bdb1=sisobnds(1,w,W1,P,[],nom,[],[],phase); %margins');
drawnow
W1=1.1; % robust margins weight
phase = 0:-5:-360;  % use a finer phase grid
bdb1 = sisobnds(1,w,W1,P,[],nom,[],[],phase);
disp('plotbnds(bdb1,1,phase); %show bounds')
drawnow
plotbnds(bdb1,1,phase),title('Robust Margins Bounds');
qpause;close(gcf);

disp(' ')
disp('bdb2=sisobnds(2,w,W2,P,nom,[],[],phase); %tracking');
drawnow
W2 = tf([1,1],2*pi*10);  % Performance weight
bdb2 = sisobnds(2,w,W2,P,[],nom,[],[],phase);
disp('plotbnds(bdb2,2,phase); %show bounds')
drawnow
plotbnds(bdb2,2,phase),title('Robust Tracking Bounds');
qpause;close(gcf);

disp(' ')
disp('bdb4=sisobnds(4,w,W4,P,[],nom,[],[],phase); %control effort');
drawnow
W4 = 3.3; % controller effort weight
bdb4 = sisobnds(4,w,W4,P,[],nom,[],[],phase);
disp('plotbnds(bdb4,4,phase); %show bounds')
drawnow
plotbnds(bdb4,4,phase),title('Robust Controller Effort Bounds');
qpause;close(gcf);

disp(' ')
disp('bdb=grpbnds(bdb1,bdb2,bdb4); %group bounds')
drawnow
bdb = grpbnds(bdb1,bdb2,bdb4);
disp('plotbnds(bdb,[],phase); %show all bounds')
drawnow
plotbnds(bdb,[],phase),title('All Bounds');
qpause;close(gcf);

disp(' ')
disp('ubdb=sectbnds(bdb); %intersect bounds')
drawnow
ubdb = sectbnds(bdb);
disp('plotbnds(ubdb,[],phase); %show bounds')
drawnow
plotbnds(ubdb,[],phase),title('Intersection of Bounds');
qpause;close(gcf);

% DESIGN
disp(' ')
disp('Design')
disp('lpshape(wl,ubdb,P0,C0,phase); %loop shaping')
drawnow
P0 = P(1,1,nom);  % nominal plant
wl = [logspace(-2,2,100),logspace(log10(110),log10(170),100),logspace(log10(171),3.3)];
C0 = tf(.0416*[1/80^2,2*.3/80,1],[1/750^2,2*.5/750,1]);
lpshape(wl,ubdb,P0,C0,phase); qpause;

disp(' ')
disp('pfshape(4,wl,W4,P,[],C0,[],F0); %pre-filter shaping')
drawnow
F0 = tf(1,[1/477,1]);
pfshape(4,wl,W4,P,[],C0,[],F0);  qpause

% ANALYSIS
disp(' ')
disp('Analysis....')

disp(' ')
disp('chksiso(1,wl,W1,P,[],C0); %margins spec')
drawnow
chksiso(1,wl,W1,P,[],C0);
qpause;close(gcf);

disp(' ')
disp('chksiso(2,wl(ind),W2/F0,P,[],C0); %tracking spec')
drawnow
ind = find(wl<=500);
chksiso(2,wl(ind),W2/F0,P,[],C0);
qpause;close(gcf);

disp(' ')
disp('chksiso(4,wl(ind),W4,P,[],C0,[],F0); %control effort spec')
drawnow
chksiso(4,wl(ind),W4,P,[],C0,[],F0);
qpause;close(gcf);
