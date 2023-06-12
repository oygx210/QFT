function pfshape(ptype,w,W,P,R,G,H,F0)
% PFSHAPE Pre-Filter IDE design environment linear systems.
%         PFSHAPE produces a Bode magnitude plot IDE pre-filter design environment
%         using default settings.
%
%         PFSHAPE(PTYPE,W,Ws,P,R,G,H,F0)  produces a Bode magnitude plot IDE
%         pre-filter design environment at the frequency set w, performance specifications Ws,
%         plant P, multiplicative uncertainty disk radius matrix R, controller complex matrices G and H,
%         and an initial filter design F0.  ALL are scalars or LTI/FRD arrays.
%
%          PTYPE specifies the closed-loop configuration.
%          PTYPE=1:  |FPGH/(1+PGH)|<Ws
%          PTYPE=2:  |F/(1+PGH)|<Ws
%          PTYPE=3:  |FP/(1+PGH)|<Ws
%          PTYPE=4:  |FG/(1+PGH)|<Ws
%          PTYPE=5:  |FGH/(1+PGH)|<Ws
%          PTYPE=6:  |FPG/(1+PGH)|<Ws
%          PTYPE=7:  Ws1<|FPG/(1+PGH)|<Ws2
%          PTYPE=8:  |FH/(1+PGH)|<Ws
%          PTYPE=9:  |FH/(1+PGH)|<Ws
%
%          In PTYPE=7, Ws is [Ws1,Ws2]
%
%         See also LPSHAPE, PUTQFT, GETQFT.

% Author: Craig Borghesani
% 9/6/93, 7/3/03 10:04AM : v2.5 updates.
% Copyright (c) 2003, Terasoft, Inc.

if nargin == 1,
 F0 = tf(1, 1);
 pfshpdef(ptype,[],[],[],[],[],[],F0);

elseif nargin==2,
 F0 = w;
 if ~isa(F0, 'lti') & ~isempty(F0),
   error('Pre-filter must be an LTI object.');
 end
 w = [];
 pfshpdef(ptype,w,[],[],[],[],[],F0);

elseif nargin==3,

 F0 = tf(1, 1);
 pfshpdef(ptype,w,W,[],[],[],[],F0);

elseif nargin==4,

 if ~isa(P, 'lti') & ~isempty(P),
   error('P must be an LTI object.');
 end

 F0 = tf(1, 1, P.Ts);

 pfshpdef(ptype,w,W,P,[],[],[],F0);

elseif nargin==5,

 if ~isa(P, 'lti') & ~isempty(P),
   error('P must be an LTI object.');
 end

 F0 = tf(1, 1, P.Ts);

 pfshpdef(ptype,w,W,P,R,[],[],F0);

elseif nargin==6,

 if ~isa(P, 'lti') & ~isempty(P),
   error('P must be an LTI object.');
 elseif ~isa(G, 'lti') & ~isempty(G),
   error('G must be an LTI object.');
 end

 F0 = tf(1, 1, P.Ts);

 pfshpdef(ptype,w,W,P,R,G,[],F0);

elseif nargin==7,

 if ~isa(P, 'lti') & ~isempty(P),
   error('P must be an LTI object.');
 elseif ~isa(G, 'lti') & ~isempty(G),
   error('G must be an LTI object.');
 elseif ~isa(H, 'lti') & ~isempty(H),
   error('H must be an LTI object.');
 end

 F0 = tf(1, 1, P.Ts);

 pfshpdef(ptype,w,W,P,R,G,H,F0);

elseif nargin==8,

 if ~isa(P, 'lti') & ~isempty(P),
   error('P must be an LTI object.');
 elseif ~isa(G, 'lti') & ~isempty(G),
   error('G must be an LTI object.');
 elseif ~isa(H, 'lti') & ~isempty(H),
   error('H must be an LTI object.');
 elseif ~isa(F0,'lti') & ~isempty(F0),
   error('F0 must be an LTI object.');

 elseif ~isempty(P) & isempty(F0),
    F0 = tf(1,1, P.Ts);

 elseif ~isempty(G) & isempty(F0),
    F0 = tf(1,1, G.Ts);

 elseif ~isempty(H) & isempty(F0),
    F0 = tf(1,1, H.Ts);

 end;%if ~isa(uF0,'lti')

 pfshpdef(ptype,w,W,P,R,G,H,F0);

else
 error('Improper number of inputs');

end
