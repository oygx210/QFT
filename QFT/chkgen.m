function [err] = chkgen(prob,w,W,A,B,C,D,G,pos)
% CHKGEN Analyze closed-loop response with respect to a specification.
%        CHKGEN(PTYPE,W,Ws,A,B,C,D,G) produces the magnitude plot
%        of the maximal closed-loop response (over arrays) and the
%        specification Ws at the frequency vector w.
%        A, B, C, D and G can be scalars and/or LTI/FRD arrays.
%
%        PTYPE specifies the closed-loop configurations.
%        PTYPE=10 produces (A+BG)/(C+DG)
%        PTYPE=11 produces (|A|+|BG|)/(C+DG)
%
%        When invoked with a left-hand argument
%        [ERR]=CHKGEN(PTYPE,w,Ws,A,B,C,D,G) returns the difference
%        between the maximal closed-loop response (over arrays) and
%        the specification Ws  at the frequency vector w.  No plot
%        is drawn.

% Author: Craig Borghesani
% 9/2/93
% 7/3/03 9:56AM : Updated for LTI support.
% Copyright (c) 2003, Terasoft, Inc.

if prob < 10 | prob > 11,
   error('PTYPE must be 10 or 11');
end

w = w(:)';

if nargin==8,
 pos=[0.333,0.28,0.6620,0.6604];
end

if isa(A,'lti'),
   if length(w) == 1 & prod(size(A)) > 1,
      A = squeeze(freqresp(A, w));
   else
      A = squeeze(freqresp(A, w)).';
   end;%if
end;%if isa(A,'lti')

if isa(B,'lti'),
   if length(w) == 1 & prod(size(B)) > 1,
      B = squeeze(freqresp(B, w));
   else
      B = squeeze(freqresp(B, w)).';
   end;%if
end;%if isa(B,'lti')

if isa(C,'lti'),
   if length(w) == 1 & prod(size(C)) > 1,
      C = squeeze(freqresp(C, w));
   else
      C = squeeze(freqresp(C, w)).';
   end;%if
end;%if isa(C,'lti')

if isa(D,'lti'),
   if length(w) == 1 & prod(size(D)) > 1,
      D = squeeze(freqresp(D, w));
   else
      D = squeeze(freqresp(D, w)).';
   end;%if
end;%if isa(D,'lti')

if isa(W,'lti'),
   if length(w) == 1 & prod(size(W)) > 1,
      W = squeeze(freqresp(W, w));
   else
      W = squeeze(freqresp(W, w)).';
   end;%if
   if any(any(isnan(W))),
      error('Frequency vector for weight inconsistent with w.');
   end
   W = abs(W);
end %if isa(W,'lti')

if isa(G,'lti'),
   if length(w) == 1 & prod(size(G)) > 1,
      G = squeeze(freqresp(G, w));
   else
      G = squeeze(freqresp(G, w)).';
   end;%if
end;%if isa(D,'lti')

[rma,cma]=size(A);[rmc,cmc]=size(C);
[rmb,cmb]=size(B);[rmd,cmd]=size(D);
[rmg,cmg]=size(G);[rmw,cmw]=size(W);
maxr=max([rma,rmb,rmc,rmd]);
maxc=max([cma,cmb,cmc,cmd]);
[rW,cW]=size(W);

% declaring sizes of replicating matricies
if repltest,
 u=ones(maxr,1); v=ones(1,maxc);
else
 u=ones(1,maxr); v=ones(1,maxc);
end

% replicate all matrices to the same size
if rma(1) == 1, A=A(u,:); end
if rmb(1) == 1, B=B(u,:); end
if rmc(1) == 1, C=C(u,:); end
if rmd(1) == 1, D=D(u,:); end
if rmg(1) == 1, G=G(u,:); end
if rmw(1) == 1, W=W(u,:); end

if cma(1) == 1, A=A(:,v); end
if cmb(1) == 1, B=B(:,v); end
if cmc(1) == 1, C=C(:,v); end
if cmd(1) == 1, D=D(:,v); end
if cmg(1) == 1, G=G(:,v); end
if cmw(1) == 1, W=W(:,v); end

upper=zeros(1,maxc);
lower=upper;
pgh=1;
F = 1;
cl = [];

if prob==10,
 for k=1:maxr,
  upper = A(k,:) + B(k,:).*G(k,:);
  lower = C(k,:) + D(k,:).*G(k,:);
  cl(k,:)=F.*(upper./lower);
 end
elseif prob==11,
 for k=1:maxr,
  upper = abs(A(k,:)) + abs(B(k,:)).*abs(G(k,:));
  lower = C(k,:) + D(k,:).*G(k,:);
  cl(k,:)=F.*(upper./lower);
 end
end

[rcl,ccl]=size(cl);
if rcl>1, clmx=max(abs(cl)); clmn=min(abs(cl));
else clmx=abs(cl); clmn=abs(cl); end
a=20*log10(W(1,:));
if nargout==0 & ccl>1,
 scrn=get(0,'screensize');sctlt=['Analysis: Problem Type ',int2str(prob)];
 f1=figure('name',sctlt,'numbertitle','off','units','norm','position',pos,...
        'vis','off','tag','qft');
 axs=[min(w),max(w)];
 if rW==1,
  axs=[axs,min(20*log10([clmx(:);W(:)]))-5,max(20*log10([clmx(:);W(:)]))+5];
 else
  axs=[axs,min(min(W-abs(cl))),max(min(W-abs(cl)))];
 end

 ax=gca;
 set(ax,'box','on',...
       'xgrid','on','ygrid','on',...
       'gridlinestyle',':',...
       'nextplot','add',...
       'xlim',axs(1:2),'ylim',axs(3:4),'xscale','log');

 if rW==1,
  semilogx(w,a,'--',w,20*log10(clmx),'-k');
  title('Weight: --');
  ylabel('Magnitude (dB)');
 else
  ermn = min(W-abs(cl));
  semilogx(w,ermn,'-k');
  title('Maximum Error');
  ylabel('Error')
 end
 xlabel('Frequency (rad/sec)');
 set(f1,'vis','on','userdata',axs);
else
 if rW==1,
  err=W(1,:)-clmx;
 else
  err=min(W-abs(cl));
 end
end
