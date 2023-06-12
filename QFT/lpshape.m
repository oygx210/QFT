function lpshape(w,bdb,P0,C0,phase)
% LPSHAPE IDE Controller design environment.
%         LPSHAPE produces a Nichols plot IDE controller design environment
%         using default settings.
%
%         LPSHAPE(P0) produces a Nichols plot IDE controller design environment
%         where P0 is an LTI/FRD model (continuous-time or discrete-time)
%         representing the nominal plant.
%
%         LPSHAPE(W,BDB,P0,C0,PHS) produces the Nichols plot IDE design
%         environment with the nominal plant P0 with a using the user-supplied
%         frequency vector W.  BDB contains the QFT bounds,
%         C0 is an initial controller (LTI model) and PHS is
%         the phase vector used to compute BDB.
%
%         See manual for a complete description of the IDE environment.
%
%         See also PFSHAPE, PUTQFT, GETQFT.

% Author: Craig Borghesani
% Date: 9/5/93
% Revised: 2/17/96 10:11PM V1.1 updates
%          4/16/03 10:10AM Only LTI support.
% Copyright (c) 2003, Terasoft, Inc.

if nargin == 0,
   P0 = zpk([], [], 1);
   C0 = zpk([], [], 1);
   lpshpdef([], [], P0, C0, []);

elseif nargin==1,
   if isa(w, 'lti'),
      C0 = zpk([], [], 1, w.Ts);
      P0 = w;
      w = [];
   else
      P0 = zpk([], [], 1);
      C0 = zpk([], [], 1);
   end
   lpshpdef(w, [], P0, C0, []);

elseif nargin==3,
   if ~isa(P0, 'lti'),
      error('P0 must be an LTI object.');
   end
   C0 = zpk([], [], 1, P0.Ts);
   lpshpdef(w, bdb, P0, C0, []);

elseif nargin==4,
   if ~isa(P0, 'lti') & ~isempty(P0),
      error('P0 must be an LTI object.');

   elseif isempty(P0) & ~isempty(C0) & isa(C0, 'lti'),
      P0 = zpk([], [], 1, C0.Ts);

   elseif isempty(P0) & isempty(C0),
      P0 = zpk([], [], 1);
      C0 = zpk([], [], 1);

   elseif isempty(C0) & ~isempty(P0),
      C0 = zpk([], [], 1, P0.Ts);

   elseif ~isa(C0, 'lti'),
      error('C0 must be an LTI object.');

   end

   if P0.Ts ~= C0.Ts,
      error('P0 and C0 must have identical sampling times.');
   end

   lpshpdef(w, bdb, P0, C0, []);

elseif nargin==5
   if ~isa(P0, 'lti') & ~isempty(P0),
      error('P0 must be an LTI object.');

   elseif isempty(P0) & ~isempty(C0) & isa(C0, 'lti'),
      P0 = zpk([], [], 1, C0.Ts);

   elseif isempty(P0) & isempty(C0),
      P0 = zpk([], [], 1);
      C0 = zpk([], [], 1);

   elseif isempty(C0) & ~isempty(P0),
      C0 = zpk([], [], 1, P0.Ts);

   elseif ~isa(C0, 'lti'),
      error('C0 must be an LTI object.');

   end

   if P0.Ts ~= C0.Ts,
      error('P0 and C0 must have identical sampling times.');
   end

   lpshpdef(w, bdb, P0, C0, phase);

else
   error('Incorrect number of inputs');

end
