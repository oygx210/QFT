function h=dintegtr(int,w,T,qflag)
% DINTEGTR Discrete-time integrator frequency response. (Utility Function)
%          DINTEGTR computes the frequency response of discrete integrators.

% Author: Craig Borghesani
% Date: 7/1/03 1:41PM
% Copyright (c) 2003, Terasoft, Inc.

h=ones(1,length(w));
z=exp(sqrt(-1)*w(:)'*T);
if int==1,
 h=(T*z./(z-1)).^(-qflag);
elseif int==2,
 h=(T^2*z./(z-1).^2).^(-qflag);
elseif int==3,
 h=((T^3*z.*(z+1))./(2*(z-1).^3)).^(-qflag);
end
