% QFTEX8 Cascaded: Outer-inner.

% Author: Y. Chait
% 5/4/94
% 4.20.03 revised for V2 (YC)
% Copyright (c) 2002 by Terasoft, Inc.

clc
clear
echo on

% Example #8 (outer-inner cascaded loop) describes the application
% of QFT to a cascaded-loop robust feedback design problem using a
% a sequential design of the outer-loop followed by the inner-loop.
% When compared to inner-outer design (Example 7), it illustrates
% the possible advantage of an improved control bandwidth distribution
% between the two loops (via the concept of "free uncertainty").

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

%     P2(s) = {k: k in [1,10]}

pause % Strike any key to continue
clc

%	The performance specifications are: design the controllers
%  G1(s) and G2(s) such that they achieve

%   1) Robust stability

%   2) Robust margins: 50 degrees phase margin in each loop


%   3) Robust output disturbance rejection
%         |Y(jw)|       |(jw)^3+64(jw)^2+748(jw)+2400|
%         |-----| < 0.02|----------------------------|, w<10
%         |D(jw)|       |    (jw)^2+14.4(jw)+169     |

%   4) Robust input disturbance rejection
%         |Y(jw)|
%         |-----| < 0.01, w<50
%         |V(jw)|
%

pause % Strike any key to continue
echo off
clc

% PROBLEM DATA
% outer loop
disp('Closing outer loop first....')
drawnow

% compute the boundary of the plant templates
disp('Computing outer plant (14 points)....')
% outer loop tf
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
%
w = [.1,1,5,10,15,25,50,500]; % working frequencies
plottmpl(w,P1), title('Outer Plant Templates')
qpause;close(gcf);

% inner loop tf
% compute the boundary of the plant templates
disp('Computing inner plant templates (3 points)....')
drawnow
c=1;
for k = linspace(1,10,5),
  P2(1,1,c) = tf(k,[1]);  c = c + 1;
end

plottmpl(w,P2), title('Inner Plant Templates')
qpause;close(gcf);

% BOUNDS
disp(' ')
disp('Computing bounds...')
disp(' ')
disp('bdb1=sisobnds(1,w,W1,P1); %margins')
drawnow
W1 = 1.2;  % margin weight
bdb1 = sisobnds(1,w,W1,P1);
disp('plotbnds(bdb1); %show bounds')
drawnow
plotbnds(bdb1),title('Robust Robust Stability Bounds');
qpause;close(gcf);

disp(' ')
disp('bdb2=genbnds(10,wbd2,W2,A,B,C,D,P1(1,1,1)); %output disturbance rejection bounds')
drawnow
wbd2 = [.1,1,5,10]; % performance frequency range
W2 = tf(0.02*[1,64,748,2400],[1,14.4,169]); % performance weight
A = 1; B = 0; C = 1; D = multmpl(P1,P2,2);
bdb2 = genbnds(10,wbd2,W2,A,B,C,D,P1(1,1,1));
disp('plotbnds(bdb2); %show bounds')
plotbnds(bdb2),title('Robust Output Disturbance Rejection Bounds');
qpause;close(gcf);

disp(' ')
disp('bdb=grpbnds(bdb1,bdb2); %group bounds')
drawnow
bdb = grpbnds(bdb1,bdb2);
disp('plotbnds(bdb); %show all bounds')
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
disp('lpshape(wl,bdb1,P10,C10); %loop shaping')
drawnow
wl = logspace(-2,2,200);  % define a frequency array for loop shaping
P10 = P1(1,1,1);
C10 = tf(2.0488e+14,[1,2.0053e+5,4.106052e+9,5.304e+11]);
lpshape(wl,bdb1,P10,C10);  qpause;

% CLOSING INNER LOOP
disp('Closing inner loop....')
drawnow
% BOUNDS
disp(' ')
disp('Computing bounds...')
disp(' ')

disp('bdb1=genbnds(10,w,W1,A,B,C,D,P2(1,1,1); %inner loop margins')
drawnow
[i1,i2] = size(P2);
W1 = 1.2;  % margin weight
A = 0; B = P2; C = 1; D = P2;
bdb1 = genbnds(10,w,W1,A,B,C,D,P2(1,1,1));
disp('plotbnds(bdb1); %show bounds')
drawnow
plotbnds(bdb1),title('Robust inner-loop margin Bounds');
qpause;close(gcf);

disp(' ')
disp('bdb2=genbnds(10,w,W1,A,B,C,D,P2(1,1,1)); %outer loop margins')
drawnow
P1P2 = multmpl(P1,P2,2);
L10 = C10*P1(1,1,1);
A = L10*P1P2; B = L10*P2(1,1,1)*P1P2; C = P1(1,1,1)*P2(1,1,1)+B;
x = P1(1,1,1)*P2(1,1,1)+C10*P1(1,1,1)*P1; D = multmpl(x,P2,2);
W2 = 1.2;
bdb2 = genbnds(10,w,W1,A,B,C,D,P2(1,1,1));
disp('plotbnds(bdb2); %show bounds')
drawnow
plotbnds(bdb2),title('Robust Outer-Loop Margin Bounds');
qpause;close(gcf);

disp(' ')
disp('bdb3=genbnds(10,wbd3,W3,A,B,C,D,P2(1,1,1)); %output dist. rej.')
W3 = tf(0.02*[1,64,748,2400],[1,14.4,169]); % performance weight
A = P1(1,1,1)*P2(1,1,1); B = A*multmpl(P2,(P1/P1),2);  % *(P1/P1) needed to have same array sizes for B and D
wbd3 = w(1:4);
bdb3 = genbnds(10,wbd3,W3,A,B,C,D,P2(1,1,1));
disp('plotbnds(bdb3); %show bounds')
drawnow
plotbnds(bdb3),title('Robust output disturbance bounds');
qpause;close(gcf);

disp(' ')
disp('bdb4=genbnds(10,wbd4,W4,A,B,C,D,P2(1,1,1)); %input dist. rej.')
drawnow
W4 = 0.01;% weight computed at w
A = P1(1,1,1)*P2(1,1,1)*P1P2;  B = 0;
wbd4 = w(1:7);
bdb4 = genbnds(10,wbd4,W4,A,B,C,D,P2(1,1,1));
disp('plotbnds(bdb4); %show bounds')
drawnow
plotbnds(bdb4),title('Robust input disturbance bounds');
qpause;close(gcf);

disp(' ')
disp('bdb=grpbnds(bdb1,bdb2,bdb3,bdb4); %group bounds')
drawnow
bdb = grpbnds(bdb1,bdb2,bdb3,bdb4);
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
disp('lpshape(wl,ubdb,P20,C20); %loop shaping')
drawnow
wl = logspace(-1,6,200);  % define a frequency array for loop shaping
P20 = P2(1,1,1); % nominal loop
C20 = tf(4e+9,[1,200400,80000000]);
lpshape(wl,ubdb,P20,C20); qpause;

% extract actual outer controller
% C1a = C10*(1+C20*P2(1,1,1))/(C20*P2(1,1,1)); % should be similar to nc1=51220; dc1=[1,130];

% ANALYSIS
disp(' ')
disp('Analysis....')

% redefine wl with less points
wl = logspace(-2,3,100);

disp(' ')
disp('chksiso(1,wl,W1,P2,[],C10); %inner margins spec')
drawnow
W1 = 1.2;
chksiso(1,wl,W1,P2,[],C10);
qpause;close(gcf);

disp(' ')
disp('chksiso(1,wl,W1,P12,[],C10); %outer margins spec')
drawnow
T2 = P2*C20/(1+P2*C20);
P12 = multmpl(P1,T2,2);
chksiso(1,wl,W1,P12,[],C10);
qpause;close(gcf);

disp(' ')
disp('chksiso(2,wl,W2,P12,[],C10); %output dist. rej. spec')
drawnow
ind = find(wl<=10);
W2 = tf(0.02*[1,64,748,2400],[1,14.4,169]); % performance weight
chksiso(2,wl(ind),W2,P12,[],C10);
qpause;close(gcf);

disp(' ')
disp('chkgen(10,wl,W3,A,B,C,D,C10); %input dist. rej. spec')
drawnow
ind = find(wl<=100);
A = P12/C20; B = 0; C = 1;  D = P12;
W3 = 0.01;
chkgen(10,wl(ind),W3,A,B,C,D,C10);
qpause;close(gcf);

% compare all 3 designs
disp(' ')
disp('Comparing controller bandwidths between single and both cascaded designs....')
drawnow
wl = logspace(1,6,50); % include higher frequencies
% single loop
Cs = tf(379*[1/42,1],[1/247^2,1/247,1]);  % single loop design from QFTEX1
mCs = 20*log10(abs(squeeze(freqresp(Cs,wl))));
% inner outer
C10 = tf(460,[1/300,1]);
C20 = tf(5,conv([1/500,1],[1/22000,1]));
mC1io = 20*log10(abs(squeeze(freqresp(C10,wl))));
mC2io = 20*log10(abs(squeeze(freqresp(C20,wl))));
% outer-inner
C10 = tf(2.0488e+14,[1,2.0053e+5,4.106052e+9,5.304e+11]);
C20 = tf(4e+9,[1,200400,80000000]);
mC1oi = 20*log10(abs(squeeze(freqresp(C10,wl))));
mC2oi = 20*log10(abs(squeeze(freqresp(C20,wl))));
h = semilogx(wl,mCs,'y',wl,mC1io,'r--',wl,mC2io,'r:',wl,mC1oi,'g--',wl,mC2oi,'g:');
legend(h, 'single loop: qftex1','Inner-Outer C1','Inner-Outer C2',...
'Outer-Inner C1','Outer-Inner C2')
qpause;close(gcf)
