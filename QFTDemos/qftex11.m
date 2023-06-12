% QFTEX11 Active vibration isolation.

% Author: Y. Chait
% 11/19/93
% 4.20.03 revised for V2 (YC)
% Copyright (c) 2002 by Terasoft, Inc.


clc
clear
echo on

% Example #11 (active vibration isolation) describes the application of
% QFT to a feedback design problem with an experimental plant model.

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

%       The plant P(s) model (mount input to measured engin acceleration)
%       is not known, however, an experimental frequency response is
%       available

pause % Strike any key to continue
echo off
clc

disp('Loading experimental data from mat file');
drawnow
load lorddata

frq1 = frf1(:,1);
mag1 = frf1(:,2);
phs1 = frf1(:,3);
frq2 = frf2(:,1);
mag2 = frf2(:,2);
phs2 = frf2(:,3);

fmin = 30;
fmid = 1000;
fmax = 6000;

s1 = find(frq1==fmid);
if isempty(s1), s1=1; end
s1 = s1 + 1;
e1 = find(frq1==fmax);
if isempty(e1), e1=max(size(frq1)); end

s2 = find(frq2==fmin);
if isempty(s2), s2=1; end
e2 = find(frq2==fmid);
if isempty(e2), e2=max(size(frq2)); end

lfrq = [frq2(s2:e2);frq1(s1:e1)];
mag = [mag2(s2:e2);mag1(s1:e1)];
phso = [phs2(s2:e2);phs1(s1:e1)];
phs = (180/pi)*unwrap((pi/180)*phso);

% plot the experimental transfer function
disp('Plotting experimental transfer function');
drawnow
subplot(211)
semilogx(lfrq,mag), grid
title('Measured Frequency Response')
subplot(212)
semilogx(lfrq,phs), grid
qpause;close(gcf);

% BOUNDS
disp(' ')
disp('Computing bounds...')
disp(' ')
disp('bdb1=sisobnds(1,wbd1,W1,Pfr); %margins')
drawnow
wl = lfrq*2*pi;
wbd1 = wl(150);
P = 10 .^(mag(:)'/20).*exp(i*pi/180*phs(:)');
Pfr = frd(P.',wl);
W1 = 2;
bdb1=sisobnds(1,wbd1,W1,Pfr);
disp('plotbnds(bdb1); %show bounds')
drawnow
plotbnds(bdb1),title('Robust Stability Bounds');
qpause;close(gcf);

disp(' ')
disp('bdb2=sisobnds(2,wbd2,W2,Pfr); %output dist. rejection')
drawnow
wbd2 = wbd1;
W2 = 0.1;
bdb2 = sisobnds(2,wbd2,W2,Pfr);
disp('plotbnds(bdb2); %show bounds')
plotbnds(bdb2),title('Robust Output Disturbance Rejection Bounds');
qpause;close(gcf);

disp(' ')
disp('bdb=grpbnds(bdb1,bdb2,bdb3); %group bounds')
drawnow
bdb = grpbnds(bdb1,bdb2);
disp('plotbnds(bdb); %show all bounds')
plotbnds(bdb),title('All Bounds');
qpause;close(gcf);

% DESIGN
disp(' ')
disp('Design')
disp('lpshape(wl,bdb,Pfr,C0); %loop shaping')
drawnow
nc0=0.77*conv([1/.75,1],conv([1/603,1],[1/2488^2,2*.36/2488,1]));
dc0=conv([1/367,1],conv([1/60^2,2*.7/60,1],[1/2027^2,2*.73/2027,1]));
C0 = tf(nc0,dc0);
lpshape(wl,bdb,Pfr,C0);
qpause;
numC=nc0; denC=dc0;

% ANALYSIS
disp(' ')
disp('Analysis....')

G=freqcp(numC,denC,wl);

disp(' ')
disp('chksiso(1,wl,W1,Pfr,[],C0); %margins spec')
drawnow
chksiso(1,wl,W1,Pfr,[],C0);
qpause;close(gcf);

disp(' ')
disp('chksiso(2,wl,W2,Pfr,[],C0); %output dist. rejection spec')
drawnow
w1 = 75*2*pi; w2 = 225*2*pi;
ind = find(wl<=w2 & wl>=w1);
chksiso(2,wl(ind),W2,Pfr,[],C0);
qpause;close(gcf);
