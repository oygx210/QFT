function putqft(fname,Ts,a1,b1,c1,d1,Terms)
% PUTQFT Interactively set QFT element matrix to mat-file.
%        PUTQFT(C) creates a dialog box that interactively allows
%        for the naming of the desired mat-file in QFT Toolbox model format.
%        C is an LTI model.  This mat-file can then be opened by the
%        IDE shaping environments LPSHAPE and PFSHAPE.
%
%        PUTQFT(FNAME,C) directly writes the passed format into the
%        file specified by the string variable FNAME (i.e. 'myfile.mat').
%
%        PUTQFT(LTI) interactively places the LTI format into a mat-file.
%
%        The shaping environments look for files with the following
%        extensions:
%              LPSHAPE - *.shp (continuous-time)
%              LPSHAPE - *.dsh  (discrete-time)
%              PFSHAPE - *.fsh (continuous-time)
%              PFSHAPE - *.dfs  (discrete-time)
%
%        See also GETQFT, LPSHAPE, PFSHAPE.

%        PUTQFT(Ts,NUM,DEN) creates a dialog box that interactively allows
%        for the naming of the desired mat-file.  Ts is the sampling time
%        and NUM and DEN represent the transfer function with descending
%        powers of s or z.  This mat-file can then be opened by any of the
%        shaping environments.
%
%        PUTQFT(Ts,Z,P,K) interactively places the Zero/Pole/Gain format
%        into a mat-file that can be opened by the shaping environments.
%
%        PUTQFT(Ts,A,B,C,D) interactively places the State-Space format into
%        a mat-file that can be opened by the shaping environments.
%
%        PUTQFT(FNAME,Ts,...) directly writes the passed format into the
%        file specified by the string variable FNAME (i.e. 'myfile.mat').
%
%        Note: Ts MUST be specified in all cases.
%              For continuous-time, Ts equals [].
%              For discrete-time, Ts > 0.

% Author: Craig Borghesani
% Date: 9/14/93
% Revised: 2/21/96 12:30 PM V1.1 updates, 7/3/03 10:06AM v2.5 updates.
% Copyright (c) 2003, Terasoft, Inc.
%       $Revision: 1.6 $

nargval = nargin;
filepath = [];

% build controller matrix from Terms input
if nargval == 7,
   TermsFields = fieldnames(Terms);
   TermsData = struct2cell(Terms);

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

   cont_r = [];
   for k = 1:length(TermsFields),
      switch lower(TermsFields{k}),
         case 'gain', % gain
            cont_r(1,1:4) = [TermsData{k}, 0, 0, 0];

         case 'discreteintegrator', % discrete integrator
            cont_r(3,1:4) = [TermsData{k}, 0, 0.5];

         case 'discretedelaypredictor', % discrete delay/predictor
            cont_r(2,1:4) = [TermsData{k}, 0, 0, 0.6];

         case 'continuousintegrator', % continuous integrator
            cont_r(2,1:4) = [TermsData{k}, 0, 0.7];

         case 'realpole', % real pole
            rp = TermsData{k};
            for k2 = 1:size(rp, 1),
               cont_r = [cont_r; rp(k2), 0, 0, 1];
            end %for

         case 'realzero', % real zero
            rz = TermsData{k};
            for k2 = 1:size(rz, 1),
               cont_r = [cont_r; rz(k2), 0, 0, 2];
            end %for

         case 'secondorderpole', % second-order pole
            sop = TermsData{k};
            for k2 = 1:size(sop, 1),
               cont_r = [cont_r; sop(k2,:), 0, 3];
            end %for

         case 'secondorderzero', % second-order zero
            soz = TermsData{k};
            for k2 = 1:size(soz, 1),
               cont_r = [cont_r; soz(k2,:), 0, 4];
            end %for

         case 'leadlag', % lead/lag
            ll = TermsData{k};
            for k2 = 1:size(ll, 1),
               cont_r = [cont_r; ll(k2,:), 0, 5];
            end %for

         case 'notch', % notch
            ntch = TermsData{k};
            for k2 = 1:size(ntch, 1),
               cont_r = [cont_r; ntch(k2,:), 6];
            end %for

         case 'complexlead', % complex lead
            cpld = TermsData{k};
            for k2 = 1:size(cpld, 1),
               cont_r = [cont_r; cpld(k2,:), 7];
            end %for

         otherwise
            error(['Term not recognized: ', TermsField{:}]);

      end %switch

   end %for
end

if ~isstr(fname),
 if nargval == 3,
  den = a1; num = Ts; Ts = fname;
 elseif nargval == 4,
  k = b1; p = a1; z = Ts; Ts = fname;
 elseif nargval == 5,
  D = c1; C = b1; B = a1; A = Ts; Ts = fname;
 elseif nargval == 1,
  [z, p, k, Ts] = zpkdata(fname);
 end
 [filename, filepath] = uiputfile('*.*','Put Element Matrix');
 nargval = nargval + 1;
else
 filename = fname;
 if nargval == 4,
  num = a1; den = b1;
 elseif nargval == 5,
  z = a1; p = b1; k = c1;
 elseif nargval == 6,
  A = a1; B = b1; C = c1; D = d1;
 elseif nargval == 2,
  [z, p, k, Ts] = zpkdata(Ts, 'v');
  if Ts == 0,
    Ts = [];
  end
 end
end

if isstr(filename),
 if nargval == 2,
   cont_r = cntpars(zpk(z,p,k,Ts));

 elseif nargval == 4,
  cont_r = cntpars(tf(num,den,Ts));

 elseif nargval == 5,
   cont_r = cntpars(zpk(z,p,k,Ts));

 elseif nargval == 6,
  cont_r = cntpars(ss2tf(A,B,C,D,Ts));

 end

 T = Ts;

 eval(['save ''',filepath,filename,''' cont_r T -mat']);

else

 disp('Warning: File not created.');

end
