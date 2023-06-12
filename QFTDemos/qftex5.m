% QFTEX5 ACC benchmark.

% Author: Y. Chait
% 2/15/93
% 4.18.03 revised for V2 (YC)
% Copyright (c) 2002 by Terasoft, Inc.


clc
clear
echo on

% Example #5 (ACC benchmark) describes the application of QFT to a
% feedback design problem with an ideal parametric uncertain flexible
% mechanical system.

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

%       The plant model P(s) a parametric uncertainty model:
%	        |                  k                                      |
%        P(s) = { -------------------------- :  m1=m2=1; k in [0.5,2] }
%	        |      m1*s^2 (m2*s^2+(1+m2/m1)k)                         |

pause % Strike any key to continue
clc

%	The performance specifications are: design a controller
%  G(s) such that it achieves

%  1) Robust stability

pause % Strike any key to continue
echo off
clc

% PROBLEM DATA

% compute plant templates
disp('Computing plant templates (100 plants)....')
drawnow
c=1;
for k = linspace(0.5,2,100),
% nump(c,:)=j;denp(c,:)=conv([1,0],conv([1,0],[1,0,2*j]));  c=c+1;
 P(1,1,c) = tf(k,conv([1,0,0],[1,0,2*k]));  c = c + 1;
end
w = [0.01,0.1,0.9,0.99,1.01,1.5,2.01,20];  % working frequencies
%
disp(' ')
disp('plottmpl(w,P); %show templates')
drawnow
plottmpl(w,P), title('Plant Templates');
qpause;close(gcf);

% BOUNDS
disp(' ')
disp('Computing bounds...')
disp(' ')
disp('bdb1=sisobnds(1,w,W1,P,[],[],[],[],phs); %margins')
drawnow
W1 = 2.25; % robust margins weight
phs = [0:-5:-30,-150:-5:-210,-320:-5:-360];  % a finer phase grid
bdb1 = sisobnds(1,w,W1,P,[],[],[],[],phs);
disp('plotbnds(bdb1,[],phs); %show bounds')
drawnow
plotbnds(bdb1,[],phs),title('Robust Margins Bounds');
qpause;close(gcf);

% DESIGN
disp(' ')
disp('Design')
disp('lpshape(wl,ubdb,P0,C0,phs); %loop shaping')
drawnow
wl=[logspace(-1.5,-.1,50),logspace(-.09,log10(2.09),50),logspace(log10(2.1),1.5,50)];
P0 = P(1,1,1);  % nominal plant
C0 = tf(.032*[10,1],[1/0.5/0.5 1/0.5 1]);
lpshape(wl,bdb1,P0,C0,phs); qpause;

% ANALYSIS
disp(' ')
disp('Analysis....')

disp(' ')
disp('chksiso(1,wl,W1,P,0,C0); %margins spec')
drawnow
chksiso(1,wl,W1,P,0,C0);
qpause;close(gcf);
