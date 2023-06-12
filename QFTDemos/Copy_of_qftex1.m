% QFTEX1 Main Example.

% Author: Y. Chait
% 10/10/94
% 4.18.03 revised for V2 (YC)
% Copyright (c) 2002 by Terasoft, Inc.

clc
clear all;
echo on

% Get current path to working directory and split
pathParts = strsplit( pwd, filesep );
% Go up one level and generate new path
src = fullfile( pathParts{1:end-1} );
% Add QFT2 to path
addpath( genpath(src) );

% Example #1 (main example) is used in the manual to describe in great
% detail the application of QFT (including use of relevant toolbox functions)
% to a feedback design problem with a parametrically uncertain plant and
% several robust performance specifications.

% Please refer to manual for more details....

% Strike any key to advance from one plot to another....
% pause % Strike any key to continue
clc

%	Consider a continuous-time, siso, negative unity feedback system

%                                   |V(s)            |D(s)
%                         ----      |      ----      |
%             ---->x---->|G(s)|---->x---->|P(s)|---->x---->
%             R(s) |      ----             ----        | Y(s)
%                  |               --                  |
%                  ---------------|-1|------------------
%                                  --

%       The plant P(s) has parametric a uncertainty model:
%	           |      k                                   |
%	   P(s)  = { ----------- :    k in [1,10], a in [1,5] }
%	           | (s+a) (s+b)      b in [20,30]            |

% pause % Strike any key to continue
clc

%	The performance specifications are: design a controller
%       G(s) such that it achieves

%  1) Robust stability

%  2) Robust margins (via closed-loop magnitude peaks)
%         | P(jw)G(jw) |
%         |------------| < 1.2,  w>0
%         |1+P(jw)G(jw)|

%  3) Robust output disturbance rejection
%         |Y(jw)|       |(jw)^3+64(jw)^2+748(jw)+2400|
%         |-----| < 0.02|----------------------------|, w<10
%         |D(jw)|       |    (jw)^2+14.4(jw)+169     |

%  4) Robust input disturbance rejection
%         |Y(jw)|
%         |-----| < 0.01, w<50
%         |V(jw)|
%

% pause % Strike any key to continue
echo off
clc

% PROBLEM DATA

% compute the boundary of the plant templates
disp('Computing plant templates (40 points)....')
drawnow

c = 1;
for k = linspace( 1, 3, 3 )
    P( 1, 1, c ) = tf( k, [ 1 0 ]);  c = c + 1;
end
nompt = 1;  % define nominal plant case

w = [ 1 2 4 8 ]; % working frequencies

disp(' ')
disp('plottmpl(w,P,nompt); %show templates')
drawnow
plottmpl(w,P,nompt), title('Plant Templates')
% qpause;close(gcf);

%% BOUNDS
% disp(' ')
% disp('Computing bounds...')
% disp(' ')
% disp('bdb1=sisobnds(1,w,W1,P,R,nompt); %margins')
% drawnow
% R=0;
% wbd1 = w; % compute bounds at all frequencies in w
% W1 = min( wbd1./5, 1.4 );  % margin weight
% bdb1 = sisobnds(1,w,W1,P,R,nompt);
% disp('plotbnds(bdb1); %show bounds')
% drawnow
% plotbnds(bdb1),title('Robust Margins Bounds');
% % qpause;close(gcf);

close all;
disp(' ')
disp('bdb2=sisobnds(2,wbd2,W2,P,R,nompt); %output disturbance rejection')
drawnow
R = 0;
wbd2 = w; % performance frequency range
W2 = min( w./5, 1.4 );% upper and lower weights
bdb2 = sisobnds(2,wbd2,W2,P,R,nompt);
disp('plotbnds(bdb2); %show bounds')
drawnow
plotbnds(bdb2),title('Robust Output Disturbance Rejection Bounds');
% qpause;close(gcf);



%% DESIGN
disp(' ')
disp('Design')
disp('lpshape(wl,ubdb,L0,G); %loop shaping')
disp(' ')
disp('There are two possible controllers here:')
disp(' ')
disp('            379*(s/42 + 1) ')
disp('  1)     --------------------    strictly proper ')
disp('         s^2/247^2+ s/247 + 1 ')
disp(' ')
disp('           379*(s/42 + 1) ')
disp('  2)       -----------------     proper ')
disp('              s/165 + 1 ')
disp(' ');
disp('  3)     Following Example 1 (Chapter 4) in manual');
disp(' ');
ans=input('Enter your choice (1,2, or 3) ==>  ');
G = tf(379*[1/42,1], [1/247^2,1/247,1]);
if ans==2,
   G = tf(379*[1/42,1], [1/165,1]);
elseif ans==3,
   G = tf(1, 1);
end
wl = logspace(-2,3,100);  % define a frequency array for loop shaping
L0=P(1,1,nompt);
L0.ioDelay = 0; % no delay
lpshape(wl,ubdb,L0,G);
qpause;

% ANALYSIS
disp(' ')
disp('Analysis....')
disp('Re-define a more dense plant template (100 points)....')
drawnow

% redfine a more dense plant template boundary - 100 points
c = 1; k = 10; b = 20;
for a = linspace(1,5,25),
 P(1,1,c) = tf(k,[1,a+b,a*b]);  c = c + 1;
end
k = 1; b = 30;
for a = linspace(1,5,25),
 P(1,1,c) = tf(k,[1,a+b,a*b]);  c = c + 1;
end
b = 30; a = 5;
for k = linspace(1,10,25),
 P(1,1,c) = tf(k,[1,a+b,a*b]);  c = c + 1;
end
b = 20; a = 1;
for k = linspace(1,10,25),
 P(1,1,c) = tf(k,[1,a+b,a*b]);  c = c + 1;
end

disp(' ')
disp('chksiso(1,wl,W1,P,R,G); %margins spec')
drawnow
chksiso(1,wl,W1,P,R,G);
qpause;close(gcf);

disp(' ')
disp('chksiso(2,wl,W2,P,R,G); %output disturbance rejection spec')
drawnow
ind=find(wl<=10);
chksiso(2,wl(ind),W2,P,R,G);
qpause;close(gcf);

disp(' ')
disp('chksiso(3,wl,W3,P,R,G); %input disturbance rejection spec')
drawnow
chksiso(3,wl(ind),W3,P,R,G);
qpause;close(gcf);
