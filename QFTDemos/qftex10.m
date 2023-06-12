% QFTEX10 Inverted pendulum.

% Author: O. Yaniv, Y. Chait
% 2/15/94
% 4.18.03 revised for V2 (YC)
% Copyright (c) 2002 by Terasoft, Inc.

clc
clear
echo on

% Example #10 (inverted pendulum) describes the application of QFT to a
% feedback design problem with a parametrically uncertain flexible
% mechanical system with both open-loop instability and pure time delay.

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

%       The plant P(s) (motor current to pendulum angle)
%       has a parametric uncertainty model:
%	        |    K*a             wn^2       s^2                   |
%	P(s)  = {------------ ----------------- -------  exp^(-s*tau) }
%	        |s(s+a)-K*k*a s^2+2*z*wn*s+wn^2 L*s^2-g               |

%       L in [0.35,0.45] m, K in [1.5,1.7] m/v/s, a in [15,17]
%       wn in [50,70] r/s, z in [0.01,0.02], tau=0.014 s
%       g=9.81 m/s^2, k=11 v/m

pause % Strike any key to continue
clc

%	The performance specifications are: design a controller
%  G(s) such that it achieves

%  1) Robust stability

%  2) Robust margins (via closed-loop magnitude peaks)
%         | P(jw)G(jw) |
%         |------------| < 2.1,  w>0
%         |1+P(jw)G(jw)|

pause % Strike any key to continue
echo off
clc

% PROBLEM DATA

disp('Computing plant templates (32 cases at 16 frequencies)...')
drawnow
% Computing plant templates
ll = [.3 .45]; kk = [1.5 1.7]; aa = [15 17]; ww = [50 70]; xx = [.01 .02]; kpf = .1; g = 9.8;
j=1;
for l=ll,
  for k=kk,
    for alfa=aa,
      for wn=ww,
        for xi=xx,
           nump = k*alfa*wn^2*(1/l)*[1 0 0];
           dddd = conv([1 alfa -kpf*k*alfa],[1 2*xi*wn wn^2]);
           denp = conv(dddd,[1 0 -g/l]);
           P(1,1,j) = tf(nump,denp);  j=j+1;
        end
      end
    end
  end
end
P.iodelay = 0.014;  % same delay for all plant cases

disp(' ')
disp('plottmpl(wb,P); %show templates')
drawnow
wb = [1 10 42 47 49 50 60 70];
plottmpl(wb,P),title('Plant Templates')
qpause;close(gcf);

% BOUNDS
disp(' ')
disp('Computing bounds...')
disp(' ')
disp('bdb1=sisobnds(1,w,W1,P); %margins');
drawnow
W1 = 2.1;
w = [0.1 0.2 0.5 1 2 5 10 30 35 40 42 47 49 50 60 70];  % working frequency range
bdb1 = sisobnds(1,w,W1,P);
disp('plotbnds(bdb1); %show bounds')
drawnow
plotbnds(bdb1),title('Robust Margins Bounds');
qpause,close(gcf);

% DESIGN
disp(' ')
disp('Design')
disp('lpshape(wl,bdb1,P0,C0); %loop shaping')
drawnow
P0 = P(1,1,1);  % nominal plant
wl = sort([linspace(45,55,100),logspace(-2,log10(150),100)]);
nc0 = [4.4971e+2,1.0312e+4,4.3164e+4,5.4188e+3,1.3846e+3];
dc0 = [1,8.6257e+1,2.9683e+3,4.8682e+2,1.0848e+2,4.5962];
C0 = tf(nc0,dc0);
lpshape(wl,bdb1,P0,C0);
qpause;

% ANALYSIS
disp(' ')
disp('Analysis....')
disp(' ')

disp('chksiso(1,wl,W1,P,[],C0); %margins spec')
drawnow
chksiso(1,wl,W1,P,[],C0);
qpause; close(gcf);
