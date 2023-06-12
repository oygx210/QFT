function [cl] = cltmpl(prob,P,G,H,F,sgn,typ)
% CLTMPL Compute various closed-loop configurations for LTI arrays.
%         CLTMPL(ptype,P,G,H,F,sgn,typ) produces the closed-loop model
%         for a user specified input-output configuration.
%
%         PTYPE specifies the closed-loop configuration.
%         PTYPE=1 computes FPGH/(1-sgnPGH)
%         PTYPE=2 computes F/(1-sgnPGH)
%         PTYPE=3 computes FP/(1-sgnPGH)
%         PTYPE=4 computes FG/(1-sgnPGH)
%         PTYPE=5 computes FGH/(1-sgnPGH)
%         PTYPE=6 computes FPG/(1-sgnPGH)
%         PTYPE=7 computes FPG/(1-sgnPGH)
%         PTYPE=8 computes FH/(1-sgnPGH)
%         PTYPE=9 computes FH/(1-sgnPGH)
%
%         P, G and H can be scalers and/or LTI/FRD arrays
%         while F cannot be an array.
%
%         sgn is the feedback sign (-1 is the default)
%
%         typ = 1 for correlated; typ = 2 for uncorrelated.
%
%         Correlated implies that the i'th model in P is matched
%         with each i'th model in all other arrays.
%         Uncorrelated implies that the i'th model in P1 is matched
%         with each model in all other arrays.
%         If all models are of same array dimension, the result
%         is computed as a correlated case.
%
%         See also ADDTMPL and MULTMPL.

% Author: Craig Borghesani
% Date: 9/18/02 9:49AM, 7/3/03 9:58AM
% Copyright (c) 2002, Terasoft, Inc.

if ~isa(P,'lti'),
   error('CLTMPL only for LTI arrays.');
end

if nargin==2,
 [P,G,H,F,sgn,typ]=clcpdef(P,[],[],[],[],[]);
elseif nargin==3,
 [P,G,H,F,sgn,typ]=clcpdef(P,G,[],[],[],[]);
elseif nargin==4,
 [P,G,H,F,sgn,typ]=clcpdef(P,G,H,[],[],[]);
elseif nargin==5,
 [P,G,H,F,sgn,typ]=clcpdef(P,G,H,F,[],[]);
elseif nargin==6,
 [P,G,H,F,sgn,typ]=clcpdef(P,G,H,F,sgn,[]);
elseif nargin==7,
 [P,G,H,F,sgn,typ]=clcpdef(P,G,H,F,sgn,typ);
else
 error('Improper number of inputs');
end

rp=prod(size(P));
rg=prod(size(G));
rh=prod(size(H));

rm=max([rp rg rh]);
pgh=1;
if prob==1,
   if typ==1,
      for p=1:rm,
         x=(rp>1)*p+(rp==1); y=(rg>1)*p+(rg==1); z=(rh>1)*p+(rh==1);
         cl(1,1,p) = F * ((P(1,1,x)*G(1,1,y)*H(1,1,z)) / (1-(sgn)*(P(1,1,x)*G(1,1,y)*H(1,1,z))));
      end
   else
      for p=1:rp, for g=1:rg, for h=1:rh,
         upper = P(1,1,p)*G(1,1,g)*H(1,1,h);
         lower = 1 - (sgn) * upper;
         cl(1,1,pgh) = F * (upper/lower);
         pgh=pgh+1;
      end; end; end
   end
elseif prob==2,
   if typ==1,
      for p=1:rm,
         x=(rp>1)*p+(rp==1); y=(rg>1)*p+(rg==1); z=(rh>1)*p+(rh==1);
         cl(1,1,p) = F * (1 / (1-(sgn)*(P(1,1,x)*G(1,1,y)*H(1,1,z))));
      end
   else
      for p=1:rp, for g=1:rg, for h=1:rh,
         upper=1;
         lower=1-(sgn)*(P(1,1,p)*G(1,1,g)*H(1,1,h));
         cl(1,1,pgh)=F*(upper/lower);
         pgh=pgh+1;
      end; end; end
   end
elseif prob==3,
   if typ==1,
      for p=1:rm,
         x=(rp>1)*p+(rp==1); y=(rg>1)*p+(rg==1); z=(rh>1)*p+(rh==1);
         cl(1,1,p) = F * (P(1,1,x) / (1-(sgn)*(P(1,1,x)*G(1,1,y)*H(1,1,z))));
      end
   else
      for p=1:rp, for g=1:rg, for h=1:rh,
         upper=P(1,1,p);
         lower=1-(sgn)*(P(1,1,p)*G(1,1,g)*H(1,1,h));
         cl(1,1,pgh)=F*(upper/lower);
         pgh=pgh+1;
      end; end; end
   end
elseif prob==4,
   if typ==1,
      for p=1:rm,
         x=(rp>1)*p+(rp==1); y=(rg>1)*p+(rg==1); z=(rh>1)*p+(rh==1);
         cl(1,1,p)=F*(G(1,1,y)/(1-(sgn)*(P(1,1,x)*G(1,1,y)*H(1,1,z))));
      end
   else
      for p=1:rp, for g=1:rg, for h=1:rh,
         upper=G(1,1,g); lower=1-(sgn)*(P(1,1,p)*G(1,1,g)*H(1,1,h));
         cl(1,1,pgh)=F*(upper/lower); pgh=pgh+1;
      end; end; end
   end
elseif prob==5,
   if typ==1,
      for p=1:rm,
         x=(rp>1)*p+(rp==1); y=(rg>1)*p+(rg==1); z=(rh>1)*p+(rh==1);
         cl(1,1,p)=F*((G(1,1,y)*H(1,1,z))/(1-(sgn)*(P(1,1,x)*G(1,1,y)*H(1,1,z))));
      end
   else
      for p=1:rp, for g=1:rg, for h=1:rh,
         upper=G(1,1,g)*H(1,1,h);
         lower=1-(sgn)*(P(1,1,p)*G(1,1,g)*H(1,1,h));
         cl(1,1,pgh)=F*(upper/lower); pgh=pgh+1;
      end; end; end
   end
elseif prob==6,
   if typ==1,
      for p=1:rm,
         x=(rp>1)*p+(rp==1); y=(rg>1)*p+(rg==1); z=(rh>1)*p+(rh==1);
         cl(1,1,p)=F*((P(1,1,x)*G(1,1,y))/(1-(sgn)*(P(1,1,x)*G(1,1,y)*H(1,1,z))));
      end
   else
      for p=1:rp, for g=1:rg, for h=1:rh,
         upper=P(1,1,p)*G(1,1,g); lower=1-(sgn)*(P(1,1,p)*G(1,1,g)*H(1,1,h));
         cl(1,1,pgh)=F*(upper/lower); pgh=pgh+1;
      end; end; end
   end
elseif prob==7,
   if typ==1,
      for p=1:rm,
         x=(rp>1)*p+(rp==1); y=(rg>1)*p+(rg==1); z=(rh>1)*p+(rh==1);
         cl(1,1,p)=F*((P(1,1,x)*G(1,1,y))/ ...
               (1-(sgn)*(P(1,1,x)*G(1,1,y)*H(1,1,z))));
      end
   else
      for p=1:rp, for g=1:rg, for h=1:rh,
         upper=P(1,1,p)*G(1,1,g); lower=1-(sgn)*(P(1,1,p)*G(1,1,g)*H(1,1,h));
         cl(1,1,pgh)=F*(upper/lower); pgh=pgh+1;
      end; end; end;
   end
elseif prob==8,
   if typ==1,
      for p=1:rm,
         x=(rp>1)*p+(rp==1); y=(rg>1)*p+(rg==1); z=(rh>1)*p+(rh==1);
         cl(1,1,p)=F*(H(1,1,z)/(1-(sgn)*(P(1,1,x)*G(1,1,y)*H(1,1,z))));
      end
   else
      for p=1:rp, for g=1:rg, for h=1:rh,
         upper=H(1,1,h); lower=1-(sgn)*(P(1,1,p)*G(1,1,g)*H(1,1,h));
         cl(1,1,pgh)=F*(upper/lower); pgh=pgh+1;
      end; end; end
   end
elseif prob==9,
   if typ==1,
      for p=1:rm,
         x=(rp>1)*p+(rp==1); y=(rg>1)*p+(rg==1); z=(rh>1)*p+(rh==1);
         cl(1,1,p)=F*((P(1,1,x)*H(1,1,z))/(1-(sgn)*(P(1,1,x)*G(1,1,y)*H(1,1,z))));
      end
   else
      for p=1:rp, for g=1:rg, for h=1:rh,
         upper=P(1,1,p)*H(1,1,h); lower=1-(sgn)*(P(1,1,p)*G(1,1,g)*H(1,1,h));
         cl(1,1,pgh)=F*(upper/lower); pgh=pgh+1;
      end; end; end
   end
end
