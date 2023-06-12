% QFTEX13 QFT tracking - discrete-time.

% Author: Y. Chait
% 2/15/93
% 4.21.03 revised for V2 (YC)
% Copyright (c) 2002 by Terasoft, Inc.

clc
clear
echo on

% Example #13 (2 DOF) describes the application of QFT to a 2
% degree-of-freedom, discrete-time, feedback design problem with a
% parametrically uncertain plant and QFT tracking specification.

% Please refer to manual for more details....

% Strike any key to advance from one plot to another....
pause % Strike any key to continue
clc

%	Consider a discrete-time, siso, negative unity feedback system

%                                       |V(s)            |D(s)
%            ----             ----      |      ----      |
%      ----->|F(s)|--->x---->|G(s)|---->x---->|P(s)|---->x---->
%      R(s)  ----      |      ----             ----      | Y(s)
%                      |             --                  |
%                      -------------|-1|------------------
%                                    --

%       The plant P(z) has parametric a uncertainty model:
%          P(z)= Z(zoh(s)*P(s))
%          Z(.) = Z-transform
%          zoh(s) = zero-order hold
%	           |    k                                |
%	   P(s)  = { ------- :  k in [1,10], a in [1,10] }
%	           | s (s+a)                             |

pause % Strike any key to continue
clc

%	The performance specifications are: given a sampling time
%       of T=0.001 sec design a controller
%       G(z) and a pre-filter F(z) such that they achieve

%  1) Robust stability

%  2) Robust margins (via closed-loop magnitude peaks)
%         | P(jw)G(jw) |
%         |------------| < 1.2,  z=exp^(i*w*T), w>0
%         |1+P(jw)G(jw)|

%  3) Robust tracking (related to tracking step responses)
%                |       P(jw)G(jw) |
%         a(w) < |F(jw) ------------| < b(w),  z=exp^(i*w*T), w<10
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

% define sampling time
Ts = 0.001;

% computing the boundary of plant templates
disp('Computing plant templates (16 plants)....')
disp('Computing zero-order hold equivalence....')
drawnow
nump =                  [0    4.998333749117734e-007    4.996667917200881e-007
                         0    1.249583437390456e-006    1.249166979189198e-006
                         0    2.499166875002956e-006    2.498333958156352e-006
                         0    3.748750312393412e-006    3.747500937345549e-006
                         0    4.998333749783868e-006    4.996667916423725e-006
                         0    4.983374916944783e-006    4.966791333993470e-006
                         0    1.245843729202889e-005    1.241697833509470e-005
                         0    2.491687458405778e-005    2.483395667007837e-005
                         0    3.737531187608667e-005    3.725093500517307e-005
                         0    4.983374916789352e-005    4.966791334037879e-005
                         0    9.993336664848584e-007    9.986676662299132e-007
                         0    2.246628793445282e-006    2.243261376988492e-006
                         0    4.486530320368942e-006    4.473090906786936e-006
                         0    9.993336665292674e-006    9.986676661410954e-006
                         0    2.246628793467487e-005    2.243261376966288e-005
                         0    4.486530320413351e-005    4.473090906698119e-005];

denp =[1.000000000000000e+000   -1.999000499833375e+000    9.990004998333750e-001
    1.000000000000000e+000   -1.999000499833375e+000    9.990004998333750e-001
    1.000000000000000e+000   -1.999000499833375e+000    9.990004998333750e-001
    1.000000000000000e+000   -1.999000499833375e+000    9.990004998333750e-001
    1.000000000000000e+000   -1.999000499833375e+000    9.990004998333750e-001
    1.000000000000000e+000   -1.990049833749168e+000    9.900498337491682e-001
    1.000000000000000e+000   -1.990049833749168e+000    9.900498337491682e-001
    1.000000000000000e+000   -1.990049833749168e+000    9.900498337491682e-001
    1.000000000000000e+000   -1.990049833749168e+000    9.900498337491682e-001
    1.000000000000000e+000   -1.990049833749168e+000    9.900498337491682e-001
    1.000000000000000e+000   -1.998001998667333e+000    9.980019986673330e-001
    1.000000000000000e+000   -1.995510109829571e+000    9.955101098295704e-001
    1.000000000000000e+000   -1.991040378772884e+000    9.910403787728836e-001
    1.000000000000000e+000   -1.998001998667333e+000    9.980019986673330e-001
    1.000000000000000e+000   -1.995510109829571e+000    9.955101098295704e-001
    1.000000000000000e+000   -1.991040378772884e+000    9.910403787728836e-001];
for l=1:16,
    P(1,1,l) = tf(nump(l,:),denp(l,:));
end
P.Ts = Ts;

disp(' ')
disp('plottmpl(w,P); %show templates')
drawnow
w = [.1,.5,1,2,15,100,1000,pi/Ts];  % working frequencies
plottmpl(w,P), title('Plant Templates')
qpause;close(gcf);

% BOUNDS
disp(' ')
disp('Computing bounds...')
disp(' ')
disp('bdb1=sisobnds(1,w,W1,P); %margins');
drawnow
W1 = 1.2; % margins weight
bdb1 = sisobnds(1,w,W1,P);
disp('plotbnds(bdb1); %show bounds');
plotbnds(bdb1),title('Robust Margins Bounds');
qpause;close(gcf);

disp(' ')
disp('bdb7=sisobnds(7,wbd7,W7,P); %tracking');
drawnow
wbd7 = [.1,.5,1,15];
mu = tf(0.6584*[1,30],[1,4,19.752]);
ml = tf(120,[1,17,82,120]);
W7 = [mu;ml]; % tracking weight
bdb7 = sisobnds(7,wbd7,W7,P);
disp('plotbnds(bdb7); %show bounds');
drawnow
plotbnds(bdb7),title('Robust Tracking Bounds');
qpause;close(gcf);

disp(' ')
disp('bdb=grpbnds(bdb1,bdb7); %group bounds')
drawnow
bdb = grpbnds(bdb1,bdb7);
disp('plotbnds(bdb); %show all bounds')
drawnow
plotbnds(bdb),title('Margins and Tracking Bounds');
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
P0 = P(1,1,1); % nominal plant
wl = [logspace(-2,log10((pi-.05)/Ts),200),logspace(log10((pi-.049)/Ts),log10(pi/Ts),100)];
C0 = tf(2235*[1,-.998],[1,-.5]);  C0.Ts = Ts;
lpshape(wl,ubdb,P0,C0);  qpause

disp(' ')
disp('pfshape(7,wl,W7,P,[],C0,[],F0); %pre-filter shaping')
drawnow
F0 = tf(0.2489027e-4*[1,0],[1,-1.9912137,0.9912386]);
F0.Ts = Ts;
pfshape(7,wl,W7,P,[],C0,[],F0); qpause

% ANALYSIS
disp(' ')
disp('Analysis....')

disp(' ')
disp('chksiso(1,wl,W1,P,R,G); %margins spec')
drawnow
chksiso(1,wl,W1,P,[],C0);
qpause;close(gcf);

disp(' ');
disp('chksiso(7,wl,W7,P,R,G,[],F); %tracking spec')
drawnow
ind = find(wl<=15);
chksiso(7,wl(ind),W7,P,[],C0,[],F0);
qpause;close(gcf);
