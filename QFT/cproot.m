function h = cproot(zta,wn,w,flag)
% CPROOT Frequency response of complex root. (Utility Function)
%        CPROOT computes the magnitude response of a complex pole or
%        zero.

% Author: Craig Borghesani
% 9/3/93
% Copyright (c) 2003, Terasoft, Inc.


if flag(2)==0,
 s=sqrt(-1)*w(:)';
 ht=s.^2/(wn+(wn==0))^2+2*zta/(wn+(wn==0))*(wn~=0)*s+(wn~=0);
 zero=find(abs(ht)==0);
 if length(zero), ht(zero)=ones(1,length(zero))*eps; end
 h=ht.^flag(1);
else
 z=exp(sqrt(-1)*w(:)'*flag(2));
 T=flag(2);
 a=zta*wn; b=wn*sqrt(1-zta^2);
% ht=z-2*exp(-a*T)*cos(b*T)+exp(-2*a*T)./z;
 ht = 1 - 2*exp(-a*T)*cos(b*T)./z + exp(-2*a*T)./(z.^2);
 zero=find(abs(ht)==0);
 if length(zero), ht(zero)=ones(1,length(zero))*eps; end
 h=ht.^flag(1)*(abs(ht(1))^(-flag(1)));
end
