% Inverted Pendulum p.471
clear;
clc;

h = .5; % m
g = 9.8; 
m = 0.05; 
M = 1.5; 

 A = [0 1 0 0 ;...
    0 0 -m*g/M 0;...
    0 0 0 1;...
    0 0 (M+m)*g/M/h 0];
B = [0 1/M 0 -1/M/h]';
C = [1 0 0 0 ;
     0 0 1 0];

sys1 = ss(A,B,C,0);

p = tf(sys1); 
p11 = p(1); % transfer from the force to the cart position
p21 = p(2); % transfer from the force to the angle of the pendulum

sys2 = ss(p21); 

% LQR control design
Q = diag([1 0 1 0]); 
R = 1; 
K = lqr(A,B,Q,R);