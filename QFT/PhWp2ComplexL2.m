function [w1,w2] = PhWp2ComplexL2(ph,wp,d);
%[w1,w2] = PhWp2ComplexL2(ph,wp,d);
%Do not allow w2/w1 >= sqrt(32) due to hardware limitations
%Finds the element
%               s^2+2*d*w1*s+w1^2
%               -------------------
%               s^2+2*d*w2*s+w2^2
%  whose maximum phase is ph(deg) and occures at wp
%
% The solution:
%               2*d1*w1*w           2*d2*w2*w
%  phi = atan( -----------) - atan( ----------)
%                w1^2-w^2            w2^2-w^2
%
% whose extremum is when its derivative is zero, that is
%
%           1                                 1
%  ---------------------------- = ----------------------------
%  1 + [2*d1*w1*w/(w1^2-w^2)]^2   1 + [2*d2*w2*w/(w2^2-w^2)]^2
%Thus
%       2*d1*w1*w/(w1^2-w^2)  = -2*d2*w2*w/(w2^2-w^2) ===>
%       d1*w1/(w1^2-w^2) = -d2*w2/(w2^2-w^2)                (*)
%                       2*d1*w1*w           - 2*d1*w1*w
%and   phi_ext = atan( -----------) - atan(-------------)
%                       w1^2-w^2              w1^2-w^2
%              = -180+2*atan(2*d1*w1*w-(w1^2-w^2))          (**)
%       ===> (for w=1) (*) gives
%       w1*w2*[d1*w2+d2*w1] = [d2*w2+d1*w1]
%       and for d1=d2=d we get w^2 = w1*w2
%
%       The solution of (**) is the solution of the quadratic equation
%       w1^2*tan(X) - 2*d*w1 - tan(X) = 0, where X=(180-ph)/2
%
%See also ZP2PhWp, PhWp2ZP

%phi = ph*pi/360; %0.5*phi
phi = (180+ph)*pi/360;
C   = -tan(phi);
A   = -C;
B   = -d;
det = B^2-A*C;
if(det<0)
   error('No solution')
else
   if(A>0) %Lead phase ph>0
      w1 = (-B+sqrt(det))/A;
   else
      w1 = (-B-sqrt(det))/A;
   end
   if(w1>5.5) w1=5.5; end   %w1/w2 bunded by about sqrt(32)
   if(w1<1/5.5) w1=1/5.5; end
   w2 = 1/w1;
   w1 = w1*wp;
   w2 = w2*wp;
end
%ph = atan(2*d*w1*wp/(w1^2-wp^2))-atan(2*d*w2*wp/(w2^2-wp^2));ph=ph*180/pi
%sys = tf([1 2*d*w1 w1^2],[1 2*d*w2 w2^2]);
%nichols(sys)
