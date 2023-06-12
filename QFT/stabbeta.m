function data = stabbeta(b, w, L)
% yc 10.20.98
% compute stabilzing beta for reset control

% known data
% w = frequency vector
% L = frequency response of nominal loop
% b = reset FORE pole (must be defined in workspace)

C=freqcp(1,[1,b],w);  % compute freq response of FORE

l=abs(L);
alpha=qatan4(L);
cr=real(C);ci=imag(C);

A=(cr+cr.*l.*cos(alpha)+l.*ci.*sin(alpha));
B=l.*cos(alpha)+l.*l;

u1=find(B>0);
u2=find(B<=0);
data=-[min(A(u1)./B(u1)), max(A(u2)./B(u2))];
