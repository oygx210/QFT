function [num,den] = c2dNotch(d1,d2,w1,w2,Ts);
%[num,den] = c2dNotch(d1,d2,w1,w2,Ts);
%c2d of the notch filter
%         s^2+2*d1*w1*s+w1^2   w2^2
%   P =  -------------------  ----
%         s^2+2*d2*w2*s+w2^2   w1^2
%First,
%translated to z-domain using bilinear transformation
%num = [4*x1*d1*w1/Ts+4*x1^2/Ts^2+w1^2 2*w1^2-8*x1^2/Ts^2 4*x1^2/Ts^2+w1^2-4*x1*d1*w1/Ts];
%      x1=w1*Ts/2; x1 = x1/tan(x1)
%den = [4*x2*d2*w2/Ts+4*x2^2/Ts^2+w2^2 2*w2^2-8*x2^2/Ts^2 4*x2^2/Ts^2+w2^2-4*x2*d2*w/2Ts];
%      x2=w2*Ts/2; x2 = x2/tan(x2)
%Input:
%d1,w1 - the complex zero damping factor and frequency
%d2,w2 - the complex pole damping factor and frequency
%Ts    - Sampling time
%z     - exp(j*wTs)
%Output:
%The transfer function num(z)/den(z)
% Author: Oded yaniv
% 12/06/03

%Note
%z+a+b/z = cos(wT)+a+b*cos(wT) + j(sin(wT)-b*sin(wT)), now 0<wT<pi and |b|<1 for stability,
%hence the phase of our expression runs from 0-deg to +pi-deg.
Tsd = Ts/2;
x1  = w1*Tsd; 
x2  = w2*Tsd; 
x1  = x1/tan(x1);%to ensure peak notch at resonance
x2  = x2/tan(x2);%to ensure peak notch at resonance
%now changes Tsd to Tsd^2
Tsd = Tsd^2;  %Ts^2/4
w12 = w1*w1;
w22 = w2*w2;
%num = [4*d1/w1*x1/Ts+4*x1^2/Ts^2/w12+1 2-8*x1^2/Ts^2/w12 4*x1^2/Ts^2/w12+1-4*d1/w1*x1/Ts];
%den = [4*d2/w2*x2/Ts+4*x2^2/Ts^2/w22+1 2-8*x2^2/Ts^2/w22 4*x2^2/Ts^2/w22+1-4*d2/w2*x2/Ts];
num = [4*d1/w1*x1/Ts+x1^2/Tsd/w12+1 2-2*x1^2/Tsd/w12 x1^2/Tsd/w12+1-4*d1/w1*x1/Ts];
den = [4*d2/w2*x2/Ts+x2^2/Tsd/w22+1 2-2*x2^2/Tsd/w22 x2^2/Tsd/w22+1-4*d2/w2*x2/Ts];
num = num/(sum(num)/sum(den));

