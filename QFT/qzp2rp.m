function [coeffs,poles,k] = qzp2rp(z,p,k1)
%RESIDUE Partial-fraction expansion or residue computation.
%	[R,P,K] = RESIDUE(Z,P,K1,Q) finds the residues, poles and direct term of
%	a partial fraction expansion of the zero, pole and gain.
%       If there are no multiple roots,
%	   B(s)       R(1)       R(2)             R(n)
%	   ----  =  -------- + -------- + ... + -------- + K(s)
%	   A(s)     s - P(1)   s - P(2)         s - P(n)
%	Vectors Z and P specify the zeros and poles.
%       The residues are returned in the column
%	vector R, the pole locations in column vector P, and the direct
%	terms in row vector K.
%       Vector Q has the format as in RESIDUE.M
%	The direct term coefficient vector K is empty if length(z) < length(p);
%	otherwise, for proper systems
%	   length(K) = length(B)-length(A)+1 = 1
%
%	If P(j) = ... = P(j+m-1) is a pole of multplicity m, then the
%	expansion includes terms of the form
%	             R(j)        R(j+1)                R(j+m-1)
%	           -------- + ------------   + ... + ------------
%	           s - P(j)   (s - P(j))^2           (s - P(j))^m
%
% ******** DOES NOT HOLD FOR REPEATED COMPLEX POLES *************
%
%	Reference:  see RESIDUE.M, RESI2.M

% Author: Yossi Chait
% 8/31/93
% Copyright (c) 2003, Terasoft, Inc.


k = [];
if length(z)==length(p), k = 1; end

% create the q matrix.  first row is the sorted poles and the second row
% designates repeated roots
[jk,i1] = sort(-abs(p));
p = p(i1);
pl = length(p);
q = zeros(1,pl);
ct = 1;

while ct<=pl,
 rt = sort(find(p(ct)==p));
 rtl = length(rt);
 q(rt) = (rtl-1):-1:0;
 if all(diff(rt)==0),
  ct = ct + rtl;
 else
  ct = ct + 1;
 end
end
q=[conj(p(:)');q];

coeffs = zeros(length(p),1);

for j1 = 1:length(p)
  p1 = q(1,j1); n = q(2,j1);
  temp = p(find(p~=p1));
  if length(k) tmp = k*prod(-p+p1); else, tmp = 0; end
  if n==0,
     coeffs(j1) = (prod(-z+p1)-tmp)/prod(-temp+p1);
  else
     v = poly(temp);
     u = poly(z) - tmp;
     for j = 1:n, [u,v] = polyder(u,v); end
       c = 1; if k < n, c = prod(1:n-k); end
       coeffs(j1) = polyval(u,p1) / prod(-temp+p1)^(2^n) / c;
  end
end

poles = p(:);

% scale by input gain
k = k1*k;
coeffs = k1*coeffs;

% make sure we have true complex conjugate pairs
indi1 = find(imag(p)>0);
indi2 = find(imag(p)<0);
if length(indi1),
  for j=1:length(indi1),
     coeffs(indi2(j)) = coeffs(indi1(j))';
  end
end

% make sure residues of real poles are zero
ind = find(imag(p)==0);
if length(ind),
     coeffs(ind) = real(coeffs(ind));
end
