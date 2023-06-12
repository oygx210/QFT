% QFTEX7 Cascaded: Inner-outer.

% Author: Y. Chait
% 7/15/93
% 4.19.03 revised for V2 (YC)
% Copyright (c) 2002 by Terasoft, Inc.

clc
clear
echo on

% Example #7 (inner-outer cascaded loop) describes the application of QFT
% to a cascaded-loop robust feedback design problem using a sequential
% design of the inner-loop followed by the outer-loop.  When compared to
% Example 1, it illustrates the advantage of cascaded-loop design over
% single-loop design in terms of reduced control bandwidth.

% Please refer to manual for more details....

% Strike any key to advance from one plot to another....
pause % Strike any key to continue
clc

%	Consider a continuous-time, siso, negative unity feedback system

%                                       |V(s)                     |D(s)
%                 -----        ------   |   -----        -----    |
%       ---->x-->|G1(s)|-->x--|G2(s)|-->x--|P2(s)|-->x--|P1(s)|-->x--->
%       R(s) |    -----    |   -----        -----    |   -----    |Y(s)
%            |             |       --                |            |
%            |             -------|-1|---------------x<--         |
%            |         --          --                  N2(s)      | N1(s)
%            ---------|-1|----------------------------------------x<--
%                      --

%       The plant P1(s) and P2(s) have parametric uncertainty models:
%	           |      1                                   |
%	   P1(s) = { ----------- :   a in [1,5], b in [20,30] }
%	           | (s+a) (s+b)                              |

%          P2(s) = {k: k in [1,10]}

pause % Strike any key to continue
clc

%	The performance specifications are: design the controllers
%  G1(s) and G2(s) such that they achieve

%  1) Robust stability

%  2) Robust margins: 50 degrees phase margin in each loop


%  3) Robust output disturbance rejection
%         |Y(jw)|       |(jw)^3+64(jw)^2+748(jw)+2400|
%         |-----| < 0.02|----------------------------|, w<10
%         |D(jw)|       |    (jw)^2+14.4(jw)+169     |

%  4) Robust input disturbance rejection
%         |Y(jw)|
%         |-----| < 0.01, w<50
%         |V(jw)|
%

pause % Strike any key to continue
echo off
clc

% PROBLEM DATA
% inner loop
disp('Closing inner loop first....')
drawnow

% compute the boundary of the plant templates
disp('Computing inner plant templates (5 points)....')
c=1;
for k = linspace(1,10,5),
  P2(1,1,c) = tf(k,[1]);  c = c + 1;
end

disp(' ')
disp('plottmpl(w,P2); %show templates')
drawnow
w = [.1,1,5,10,50,500];
plottmpl(w,P2), title('Plant Templates')
qpause;close(gcf);

% BOUNDS
disp(' ')
disp('Computing bounds...')
disp(' ')
disp('bdb1=sisobnds(1,w,W1,P2); %margins')
drawnow
W1 = 1.2;  % margin weight
bdb1 = sisobnds(1,w,W1,P2);
disp('plotbnds(bdb1); %show bounds')
drawnow
plotbnds(bdb1),title('Robust Stability Bounds');
qpause;close(gcf);

% DESIGN
disp(' ')
disp('Design')
disp('lpshape(wl,bdb1,P20,C20); %loop shaping')
drawnow
wl = logspace(-1,6,100);  % define a frequency array for loop shaping
P20 = P2(1,1,1); % nominal loop
C20 = tf(5,conv([1/500,1],[1/22000,1]));
lpshape(wl,bdb1,P20,C20); qpause;

% ANALYSIS
disp(' ')
disp('Analysis....')

disp(' ')
disp('chksiso(1,wl,W1,P2,[],C20); %margins spec')
drawnow
chksiso(1,wl,W1,P2,[],C20);
qpause;close(gcf);

% OUTER LOOP
disp(' ')
disp('Closing outer loop....')
disp(' ')
drawnow
% compute the boundary of the plant templates
disp('Computing equivalent plant templates (5*14 points)....')
aa = linspace(1,5,5);
bb = linspace(20,30,3);
c = 1; b = bb(1);
for a = aa,
   P1(1,1,c) = tf(1,[1,a+b,a*b]);  c = c + 1;
end
b = bb(3);
for a = aa,
   P1(1,1,c) = tf(1,[1,a+b,a*b]);  c = c + 1;
end
a = aa(1);
for b = bb(2:3),
   P1(1,1,c) = tf(1,[1,a+b,a*b]);  c = c + 1;
end
a = aa(5);
for b = bb(2:3),
   P1(1,1,c) = tf(1,[1,a+b,a*b]);  c = c + 1;
end

% compute closed-loop transfer function of inner loop
T2 = P2*C20/(1+P2*C20);

% compute product of 2nd plant and closed-loop transfer functions
P12 = multmpl(P1,T2,2);  % uncorrelated uncertainties

disp(' ')
disp('plottmpl(w,P12); %show templates')
drawnow
plottmpl(w,P12), title('Plant Templates')
qpause;close(gcf);

% BOUNDS
disp(' ')
disp('Computing bounds...')
disp(' ')
disp('bdb1=sisobnds(1,w,W1,P12); %margins')
drawnow
W1 = 1.2;  % define weight
bdb1 = sisobnds(1,w,W1,P12);
disp('plotbnds(bdb1); %show bounds')
drawnow
plotbnds(bdb1),title('Robust Stability Bounds');
qpause;close(gcf);

disp(' ')
disp('bdb2=sisobnds(2,wbd2,W2,P12); %output disturbance rejection')
drawnow
wbd2 = [.1,1,5,10]; % performance frequency range
W2 = tf(0.02*[1,64,748,2400],[1,14.4,169]); % performance weight
bdb2 = sisobnds(2,wbd2,W2,P12);
disp('plotbnds(bdb2); %show bounds')
drawnow
plotbnds(bdb2),title('Robust Output Disturbance Rejection Bounds');
qpause;close(gcf);

disp('bdb3=genbnds(10,wbd3,W3,A,B,C,D,P12(1,1,1)); %input disturbance rejection')
drawnow
wbd3 = [.1,1,5,50]; % performance frequency rangeay
W3 = 0.01; % performance weight
A = P12/C20; B = 0;  C = 1;  D = P12;
bdb3 = genbnds(10,wbd3,W3,A,B,C,D,P12(1,1,1));
disp('plotbnds(bdb3); %show bounds');
drawnow
plotbnds(bdb3),title('Robust Input Disturbance Rejection Bounds');
qpause;close(gcf);

disp(' ')
disp('bdb=grpbnds(bdb1,bdb2,bdb3); %group bounds')
drawnow
bdb = grpbnds(bdb1,bdb2,bdb3);
disp('plotbnds(bdb)')
drawnow
plotbnds(bdb),title('All Bounds');
qpause;close(gcf);

disp(' ')
disp('ubdb=sectbnds(bdb); %intersect bounds')
drawnow
ubdb = sectbnds(bdb);
disp('plotbnds(ubdb)')
drawnow
plotbnds(ubdb),title('Intersection of Bounds');
qpause;close(gcf);

% DESIGN
disp(' ')
disp('Design')
disp('lpshape(wl,ubdb,P120,C10); %loop shaping')
drawnow
wl = logspace(-2,log10(100),200);  % define a frequency array for loop shaping
P120 = P12(1,1,1); % nominal loop
C10 = tf(460,[1/300,1]);
lpshape(wl,ubdb,P120,C10); qpause;

% ANALYSIS
disp(' ')
disp('Analysis....')

% redefine wl to reflect bandwidth of interest
wl = logspace(-1,3,150);

disp(' ')
disp('chksiso(1,wl,W1,P12,[],C10); %margins spec')
drawnow
chksiso(1,wl,W1,P12,[],C10);
qpause;close(gcf);

disp(' ')
disp('chksiso(2,wl(ind),W2,P12,[],C10); %output disturbance rejection spec')
drawnow
ind = find(wl<=10);
chksiso(2,wl(ind),W2,P12,[],C10);
qpause;close(gcf);

disp(' ')
disp('hkgen(10,wl(ind),W3,A,B,C,D,C10); %input disturbance rejection')
drawnow
ind = find(wl<=100);
A = P12/C20; B = 0;  C = 1;  D = P12;
chkgen(10,wl(ind),W3,A,B,C,D,C10);
qpause;close(gcf);

disp(' ')
disp('Comparing controller bandwidths between single and inner-outer cascased designs....')
drawnow
Cs = tf(379*[1/42,1],[1/247^2,1/247,1]);  % single loop design from QFTEX1
wl = logspace(-1,6,100); % show high-frequency dynamics
mCs = 20*log10(abs(squeeze(freqresp(Cs,wl))));
mC1 = 20*log10(abs(squeeze(freqresp(C10,wl))));
mC2 = 20*log10(abs(squeeze(freqresp(C20,wl))));
%semilogx(wl,mCs,'y',wl,mC1,'g',wl,mC2,'r')
%title('Comparing controller bandwidths')
%text(wl(55),mCs(55),'single loop: qftex1')
%text(wl(1),mC1(1),'cascaded-outer loop')
%text(wl(1),mC2(1),'cascaded-inner loop')
h = semilogx(wl,mCs,'y',wl,mC1,'g',wl,mC2,'r');
legend(h, 'single loop: qftex1', 'cascaded-outer loop', 'cascaded-innerloop')
qpause;close(gcf);
