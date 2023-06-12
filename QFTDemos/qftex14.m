% QFTEX14 CD mechanism - sampled-data.

% Author: Y. Chait
% 1/17/94
% 4.26.03 revised for V2 (YC)
% Copyright (c) 2002 by Terasoft, Inc.


clc
clear
echo on

% Example #14 (compact disk) describes the application of QFT to a
% high performance feedback design problem with a nonminimum phase
% zero and an experimental plant model including parametric uncertainty.

% Please refer to manual for more details....

% Strike any key to advance from one plot to another....
pause % Strike any key to continue
clc

%	Consider a continuous-time (plant, zero-order-hold, sampling and
%       computational delay all lumped together), siso,
%       negative unity feedback system

%                   E(s)  ----             ----
%             ---->x---->|G(s)|---------->|P(s)|---------->
%             R(s) |      ----             ----        | Y(s)
%                  |               --                  |
%                  ---------------|-1|------------------
%                                  --

%       The plant P(s) model is not known, however, a nominal frequency
%       response has been obtained experimentally.  The three significant
%       natural frequencies are allowed to vary by 5% from nominal values.
%       The computed (w/o identification) uncertain frequency response
%       model is

pause % Strike any key to continue
echo off
clc

% PROBLEM DATA
load sisocd
Pw=frdata(P);

% define frequency array for displaying plant templates
in = [25,81,113,132,154,187,241];
nompt = 63;

disp(' ')
disp('Plotting 125 cases....')
drawnow
subplot(211),loglog(wp,abs(squeeze(Pw))), title('Magnitude');
subplot(212),semilogx(wp,180/pi*qatan4(squeeze(Pw))), title('Phase');
qpause;close(gcf);
disp(' ')
disp('Plotting plant templates (pre-computed)....')
disp('plottmpl(wp(in),P,nompt)')
drawnow
plottmpl(wp(in),P,nompt), title('Plant Templates')
qpause;close(gcf);

% BOUNDS
disp(' ')
disp('Computing bounds (pre-computed)...')
disp(' ')

disp('bdb1=sisobnds(1,wp,wbd1,W1,P,[],nompt,[],[],phs); %margins')
drawnow
W1 = 3;  % define weight
%bdb1 = sisobnds(1,wp(in),W1,P,[],nompt,[],[],phs);
disp('plotbnds(bdb1,[],phs); %show bounds')
plotbnds(bdb1,[],phs),title('Robust Margins Bounds');
qpause;close(gcf);

disp(' ')
disp('bdb2=sisobnds(2,wp(in),W2,P,[],nompt,[],[],phs); %sensitivity')
drawnow
wbd2 = wp(in); % the frequency array
%bdb2 = sisobnds(2,wp(in),W2,P,[],nompt,[],[],phs);
disp('plotbnds(bdb2,[],phs); %show bounds')
plotbnds(bdb2,[],phs),title('Robust Sensitivity Bounds');
qpause;close(gcf);

disp(' ')
disp('bdb=grpbnds(bdb1,bdb2); %group bounds')
drawnow
bdb = grpbnds(bdb1,bdb2);
disp('plotbnds(bdb,[],phs); %show all bounds')
plotbnds(bdb,[],phs),title('All Bounds');
qpause;close(gcf);

disp(' ')
disp('ubdb=sectbnds(bdb); %intersecting bounds')
drawnow
ubdb = sectbnds(bdb);
disp('plotbnds(ubdb,[],phs); %show bounds')
drawnow
plotbnds(ubdb,[],phs),title('Intersection of Bounds');
qpause;close(gcf);

% DESIGN
disp(' ')
disp('Design')
disp('lpshape(wl,ubdb,P(:,:,nompt),C0,phs); %loop shaping')
drawnow
wl = wp;  % define a frequency array for loop shaping
C0 = getqft('sisocdc.shp');
lpshape(wl,ubdb,P(:,:,nompt),C0,phs);
qpause

% ANALYSIS
disp(' ')
disp('Analysis....')
disp(' ')

disp(' ')
disp('chksiso(1,wl,W1,P,[],C0); %sensitivity spec')
drawnow
chksiso(1,wl,W1,P,[],C0);
qpause;close(gcf);

disp(' ')
disp('chksiso(2,wl,W2,P,[],C0); %sensitivity spec')
drawnow
chksiso(2,wl,W2,P,[],C0);
qpause;close(gcf);
