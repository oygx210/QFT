function [a,b,c,d,t,Terms]=getqft(fname)
% GETQFT Interactively retrieve QFT element matrix from mat-file.
%        C = GETQFT places the contents of the chosen file into the
%        command line LTI model named C.
%
%        C = GETQFT(FNAME) reads directly from the filename specified
%        by FNAME.  FNAME must be in quotes, i.e., 'filename.shp'.
%
%        FNAME can also be one of the file extensions associated with
%        a shaping environment. i.e. '*.shp'.
%
%        The shaping environments create files with the following
%        extensions:
%              LPSHAPE - *.shp (continuous-time)
%              LPSHAPE - *.dsh  (discrete-time)
%              PFSHAPE - *.fsh (continuous-time)
%              PFSHAPE - *.dfs  (discrete-time)
%
%        See also PUTQFT, LPSHAPE, PFSHAPE.

%        [NUM,DEN]=GETQFT places the contents of the chosen file into the
%        numerator/denominator format.  Elements in NUM and DEN are in
%        descending powers of s.
%
%        [Z,P,K]=GETQFT places the contents of the chosen file into the
%        zero/pole/gain format.  Z contains the zeros, P contains the poles,
%        and K contains the D.C. gain.
%
%        [A,B,C,D]=GETQFT places the contents into a balanced form.
%
%        [...,Ts]=GETQFT returns the desired discrete-time format and the
%        sampling-time, Ts.
%


% Author: Craig Borghesani
% Date: 9/14/93
% Revised: 2/16/96 1:18 PM V1.1 updates, 7/3/03 9:59AM v2.5 updates.
% Copyright (c) 2003, Terasoft, Inc.

% V5 initialization
a = []; b = []; c = []; d = []; t = [];

if nargin==0,
 fname = '*.*';
 [filename,pathname]=uigetfile(fname,'Get Element Matrix');
 if isstr(filename),
   eval(['load ''',pathname,filename,''' -mat']);
 else
   return;
 end

elseif strcmp(fname,'*.shp') | strcmp(fname,'*.dsh') | ...
       strcmp(fname,'*.fsh') | strcmp(fname,'*.dfs'),
 [filename,pathname]=uigetfile(fname,'Get Element Matrix');
 if isstr(filename),
   eval(['load ''',pathname,filename,''' -mat']);
 else
   return;
 end

else
 eval(['load ''',fname,''' -mat']);

end

Ts = T;
if isempty(Ts), Ts = 0; end

des_form = '';

if nargout == 1,
 des_form = 'lti';
end

if nargout == 2,
 des_form = 'num';
 if Ts > 0,
  disp('Retrieving discrete-time numerator/denominator');
 else
  disp('Retrieving continuous-time numerator/denominator');
 end
end

if nargout == 3 & Ts > 0,
 des_form = 'num';
 disp('Retrieving discrete-time numerator/denominator');
elseif nargout == 3,
 des_form = 'zpk';
 disp('Retrieving continuous-time zero/pole/gain');
end

if nargout == 4 & Ts > 0,
 des_form = 'zpk';
 disp('Retrieving discrete-time zero/pole/gain');

elseif nargout == 4,
 des_form = 'abcd';
 disp('Retrieving continuous-time state-space form');
end

if nargout == 5,
 des_form = 'abcd';
 disp('Retrieving discrete-time state-space form');
end

% complex lead: 7
% notch: 6
% lead/lag: 5
% second-order zero: 4
% second-order pole: 3
% real zero: 2
% real pole: 1
% continuous integrator: 0.7
% discrete delay/predictor: 0.6
% discrete integrator: 0.5
% gain: 0

for k = 1:size(cont_r, 1),
   switch cont_r(k,4)
      case 0, % gain
         Terms.Gain = cont_r(k, 1);

      case 0.5, % discrete integrator
         Terms.DiscreteIntegrator = cont_r(k, 1:2);

      case 0.6, % discrete delay/predictor
         Terms.DiscreteDelayPredictor = cont_r(k, 1);

      case 0.7, % continuous integrator/differentiator
         Terms.ContinuousIntegrator = cont_r(k, 1:2);

      case 1, % real pole
         if isfield(Terms, 'RealPole'),
            Terms.RealPole = [Terms.RealPole; cont_r(k, 1)];
         else
            Terms.RealPole = cont_r(k, 1);
         end

      case 2, % real zero
         if isfield(Terms, 'RealZero'),
            Terms.RealZero = [Terms.RealZero; cont_r(k, 1)];
         else
            Terms.RealZero = cont_r(k, 1);
         end

      case 3, % second-order pole
         if isfield(Terms, 'SecondOrderPole'),
            Terms.SecondOrderPole = [Terms.SecondOrderPole; cont_r(k, 1:2)];
         else
            Terms.SecondOrderPole = cont_r(k, 1:2);
         end

      case 4, % second-order zero
         if isfield(Terms, 'SecondOrderZero'),
            Terms.SecondOrderZero = [Terms.SecondOrderZero; cont_r(k, 1:2)];
         else
            Terms.SecondOrderZero = cont_r(k, 1:2);
         end

      case 5, % lead/lag
         if isfield(Terms, 'LeadLag'),
            Terms.LeadLag = [Terms.LeadLag; cont_r(k, 1:2)];
         else
            Terms.LeadLag = cont_r(k, 1:2);
         end

      case 6, % notch
         if isfield(Terms, 'Notch'),
            Terms.Notch = [Terms.Notch; cont_r(k, 1:3)];
         else
            Terms.Notch = cont_r(k, 1:3);
         end

      case 7, % complex lead
         if isfield(Terms, 'ComplexLead'),
            Terms.ComplexLead = [Terms.ComplexLead; cont_r(k, 1:3)];
         else
            Terms.ComplexLead = cont_r(k, 1:3);
         end

   end %switch

end %for

if strcmp(des_form,'num'),
 [a,b]=cntextr(cont_r,Ts);
 c = T;
elseif strcmp(des_form,'zpk'),
 [a,b,c]=cnt2zpk(cont_r,Ts);
 d = T;
elseif strcmp(des_form,'abcd'),
 [z,p,k]=cnt2zpk(cont_r,Ts);
 if T == 0,
  [do_it,repeat,jk] = chkzp(z,p,Ts);
  if do_it,
   if repeat,
    [sysb,hsv]=qfwbal(z,p,k,[],'z');
   else
    [r,p,k]=qzp2rp(z,p,k);
    [sysb,hsv]=qfwbal(r,p,k,[],'r');
   end
   [sr,sc] = size(sysb);
   a = sysb(1:sr-2,1:sc-2);
   b = sysb(1:sr-2,sc-1);
   c = sysb(sr-1,1:sc-2);
   d = sysb(sr-1,sc-1);
  else
   [a,b,c,d]=zp2ss(z,p,k);
  end
 else
  [a,b,c,d]=zp2ss(z,p,k,Ts);
 end
 t = T;
elseif strcmp(des_form, 'lti'),
 [z,p,k]=cnt2zpk(cont_r,Ts);
 a = zpk(z,p,k,Ts);
end
