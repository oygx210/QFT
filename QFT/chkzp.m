function [do_it,repeat,msg]=chkzp(z,p,T)
% CHKZP Check zeros and poles. (Utility Function)
%       CHKZP determines if the poles and zeros selected by the user
%       are proper and stable before the call to QFWBAL.
%       It also checks for repeated real or complex poles.

% Author: Craig Borghesani
% Date: 10/29/93
% Revised: 2/16/96 1:03 PM V1.1 updates
% Copyright (c) 2003, Terasoft, Inc.
%       $Revision: 1.5 $

%%%%%% V5 change to accomodate nargin change
nargval = nargin;

if nargval==3,
 if T == 0,
  nargval=2;
 end
end

str1=[]; str2=[]; str3=[]; msg = [];

repeat=0;
do_it=1;

% check for properness
lz=length(z); lp=length(p);
if lz>lp,
 do_it=0;
 str1 = '[Improper]';
end

if lz > 1 | lp > 1,
 % check for stability
 if nargval==2,
  axis_poles=find(real(p)==0);
  unstable_poles=find(real(p)>0);
 else
  axis_poles=find(real(p)==1);
  unstable_poles=find(real(p)>1);
 end

 if length(axis_poles),
  do_it=0;
  str2 = '[Imag axis poles]';
 end

 if length(unstable_poles),
  do_it=0;
  str3 = '[Unstable poles]';
 end

 if ~do_it,
  if nargout==2,
   errordlg(['Cannot be balanced due to ',str1,str2,str3],'Message','on');
  else
   disp(['Cannot be balanced due to ',str1,str2,str3]);
  end
 end

% check for repeated poles
 if lp & (~length(unstable_poles)) & (~length(axis_poles)),
  sp_real = sort(p(find(imag(p)==0)));
  sp_imag = sort(p(find(imag(p)~=0)));

  if any(diff([0;sp_real])==0),
   repeat = 1;
  end

  if any([diff([0;real(sp_imag)])==0 & diff([0;imag(sp_imag)])==0]),
   repeat = 1;
  end
 end
else
 do_it = 0;
 if nargout==2,
  errordlg('No further reduction possible','Message','on');
 elseif length(lp)==1,
  if p < 0, do_it = 1; end
 end
end
