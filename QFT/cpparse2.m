function [out_cp,cp] = cpparse(in_bnd,phs_deg)
% CPPARSE Complex parsing of QFT bound information
%
%         OUT_CP = CPPARSE(IN_BND) parses the single column
%         QFT bound in IN_BND and uses the default phase vector.
%
%         OUT_CP = CPPARSE(IN_BND,PHS) parse the single column
%         QFT bound in IN_BND and uses the phase vector in
%         PHS (degrees) (a row vector).
%
%         This function is entirely experimental and to be used
%         by Dr. Yossi Chait only.

% Author: Craig Borghesani <cborg@terasoft.com>
% Date: 7/23/96 1:11AM
% Revision:
% Copyright (c) 1996 by Terasoft, Inc.

% obtain default phase vector
if nargin == 1,
   defs = qftdefs;
   phs_deg = defs(1,1):defs(1,2):defs(1,3);
end

phs_deg = phs_deg';
ld = length(phs_deg);

% determine location of bounds in bound vector
[loc_abv, loc_blw] = wherebnd(in_bnd);

abv_vec = [];
for k = 1:size(loc_abv),
   abv_tmp = loc_abv(k,2):loc_abv(k,3);
   abv_vec = [abv_tmp,abv_vec];
end

% if no above bound
if ~sum(abv_vec), abv_vec = []; end

blw_vec = [];
for k = 1:size(loc_blw),
   blw_tmp = loc_blw(k,2):loc_blw(k,3);
   blw_vec = [blw_tmp,blw_vec];
end

% if no below bound
if ~sum(blw_vec), blw_vec = []; end

% compute complex values for above bound
abv_cplx = mp2cp(10 .^(in_bnd(abv_vec)/20),phs_deg(abv_vec));

% compute complex values for below bound
blw_cplx = mp2cp(10 .^(in_bnd(blw_vec+ld)/20),phs_deg(blw_vec));

%%% combine complex vectors
out_cp = abv_cplx;

% if there is an above bound and there is no below bound, the
% next line will not change out_cp.  however, if there is a
% below bound (and since we have an above bound), the below
% bound is flipped such as to create the "circle" created
% in the nichols chart.  this can be verified by plotting the
% resulting complex matrix with plot(real(out_cp),imag(out_cp))
if length(abv_cplx),
   if length(abv_tmp) < 72,
      i1=length(abv_vec);
      mr2=in_bnd(blw_vec(1)); 
      mr1=in_bnd(blw_vec(1)+ld);
      ml1=in_bnd(blw_vec(i1)+ld);
      ml2=in_bnd(blw_vec(i1));
      dl=(ml2-ml1)/10;
      dr=(mr2-mr1)/10;
      mr=[mr2:-dr:mr1]';
      ml=[ml2:-dl:ml1]';
      abv_cplxl = mp2cp(10 .^(ml/20),ones(11,1)*phs_deg(abv_vec(i1)));
      abv_cplxr = mp2cp(10 .^(mr/20),ones(11,1)*phs_deg(abv_vec(1)));
      out_cp = abv_cplx;
      out_cp = [abv_cplxr;out_cp;abv_cplxl;flipud(blw_cplx)];
   else
      out_cp = [out_cp;flipud(blw_cplx)];
   end
else

% if there is no above bound, we have a below bound which spans
% the entire phase range.  no flipping necessary.
   out_cp = blw_cplx;
end

temp = out_cp ./ (1 + out_cp);

cp = [real(out_cp),imag(out_cp)];
out_cp = [real(temp),imag(temp)];
