function [err,cl] = chksiso(prob,w,W,P,R,G,H,F,pos)
% CHKSISO Analyze closed-loop response with respect to specifications.
%         CHKSISO(PTYPE,w,Ws,P,R,G,H,F) produces the magnitude plot
%         of the closed-loop response and the specification W over the frequency vector
%         w.  P, G, H, and F can be either complex or LTI object arrays.
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
%         When invoked with a left-hand argument
%         [ERR]=CHKSISO(PTYPE,W,WS,P,R,G,H,F) returns the difference
%         between the closed-loop and WS in the vector ERR.  No plot
%         is drawn.

% Author: Craig Borghesani
% 9/2/93
% 7/3/03 9:57AM : Updated for LTI support.
% Copyright (c) 2003, Terasoft, Inc.

%%%%%% V5 initialization
err = [];

flag=1;
if nargout==2,
 flag=0;
end
if nargin==4,
 [w,W,P,R,G,H,F,pos]=chkdef(prob,w,W,P,[],[],[],[],[],flag);
elseif nargin==5,
 [w,W,P,R,G,H,F,pos]=chkdef(prob,w,W,P,R,[],[],[],[],flag);
elseif nargin==6,
 [w,W,P,R,G,H,F,pos]=chkdef(prob,w,W,P,R,G,[],[],[],flag);
elseif nargin==7,
 [w,W,P,R,G,H,F,pos]=chkdef(prob,w,W,P,R,G,H,[],[],flag);
elseif nargin==8,
 [w,W,P,R,G,H,F,pos]=chkdef(prob,w,W,P,R,G,H,F,[],flag);
elseif nargin==9,
 [w,W,P,R,G,H,F,pos]=chkdef(prob,w,W,P,R,G,H,F,pos,flag);
else
 error('Improper number of inputs');
end

[rmp,cmp]=size(P); [rmg,cmg]=size(G); [rmh,cmh]=size(H);
[rW,cW]=size(W);
if rW~=0,
 if cW == 1, W=W(:,ones(1,cmp)); end
end

q=all(diff([rmp rmg rmh])==0);
uP=abs(P); uG=abs(G); uH=abs(H);
vP=qatan4(P); vG=qatan4(G); vH=qatan4(H);
if q, upper=zeros(1,cmp); lower=upper; end
pgh=1;
if prob==1,
 if q,
  for k=1:rmp,
   lra=find(R(k,:)==0); lrb=find(R(k,:)~=0 & R(k,:)~=1); lrc=find(R(k,:)==1);
   upper(lra)=P(k,lra).*G(k,lra).*H(k,lra);
   lower(lra)=1+P(k,lra).*G(k,lra).*H(k,lra);
   m=uP(k,lrb).*uG(k,lrb).*uH(k,lrb); mb=1 ./(m.*abs(1-R(k,lrb).^2));
   Rb=R(k,lrb).*mb; phi=pi*(R(k,lrb)>1)-(vP(k,lrb)+vG(k,lrb)+vH(k,lrb));
   upper(lrb)=ones(1,length(lrb)); lower(lrb)=abs(abs(1+mb.*exp(i*phi))-Rb);
   upper(lrc)=2*uP(k,lrc).*uG(k,lrc).*uH(k,lrc); phi=vP(k,lrc)+vG(k,lrc)+vH(k,lrc);
   if length(lrc), lower(lrc)=abs(upper(lrc).*cos(phi)+1); end
   cl(k,:)=F.*(upper./lower);
  end
 else
  for p=1:rmp, for g=1:rmg, for h=1:rmh,
   if any(R(p,:)==0),
    lr=find(R(p,:)==0);
    upper(1,lr)=F.*P(p,lr).*G(g,lr).*H(h,lr); lower(1,lr)=1+upper;
   end
   if any(R(p,:)~=1 & R(p,:)~=0),
    lr=find(R(p,:)~=1 & R(p,:)~=0);
    m=uP(p,lr).*uG(g,lr).*uH(h,lr); mb=1 ./(m.*abs(1-R(p,lr).^2));
    Rb=R(p,lr).*mb; phi=pi*(R(p,lr)>1)-(vP(p,lr)+vG(g,lr)+vH(h,lr));
    upper(1,lr)=ones(1,length(lr)); lower(1,lr)=abs(abs(1+mb.*exp(i*phi))-Rb);
   end
   if any(R(p,:)==1),
    lr=find(R(p,:)==1);
    upper(1,lr)=2*uP(p,lr).*uG(g,lr).*uH(h,lr); phi=vP(p,lr)+vG(g,lr)+vH(h,lr);
    lower(1,lr)=abs(upper.*cos(phi)+1);
   end
   cl(pgh,:)=F.*(upper./lower); pgh=pgh+1;
  end; end; end
 end
elseif prob==2,
 if q,
  for k=1:rmp,
   lra=find(R(k,:)==0); lrb=find(R(k,:)~=0);
   upper(lra)=ones(1,length(lra));
   lower(lra)=1+P(k,lra).*G(k,lra).*H(k,lra);
   upper(lrb)=ones(1,length(lrb));
   lower(lrb)=abs(abs(1+P(k,lrb).*G(k,lrb).*H(k,lrb))-...
              abs(P(k,lrb).*G(k,lrb).*H(k,lrb)).*R(k,lrb));
   cl(k,:)=F.*(upper./lower);
  end
 else
  for p=1:rmp, for g=1:rmg, for h=1:rmh,
   if any(R(p,:)==0),
    lr=find(R(p,:)==0);
    upper(1,lr)=ones(1,length(lr));
    lower(1,lr)=1+P(p,lr).*G(g,lr).*H(h,lr);
   end
   if any(R(p,:)~=0),
    lr=find(R(p,:)~=0);
    upper(1,lr)=ones(1,length(lr));
    lower(1,lr)=abs(abs(1+P(p,lr).*G(g,lr).*H(h,lr))-...
          abs(P(p,lr).*G(g,lr).*H(h,lr)).*R(p,lr));
   end
   cl(pgh,:)=F.*(upper./lower); pgh=pgh+1;
  end; end; end
 end
elseif prob==3,
 if q,
  for k=1:rmp,
   lra=find(R(k,:)==0); lrb=find(R(k,:)~=0 & R(k,:)~=1); lrc=find(R(k,:)==1);
   upper(lra)=P(k,lra);
   lower(lra)=1+P(k,lra).*G(k,lra).*H(k,lra);
   m=uP(k,lrb).*uG(k,lrb).*uH(k,lrb); mb=1 ./(m.*abs(1-R(k,lrb).^2));
   Rb=R(k,lrb).*mb; phi=pi*(R(k,lrb)>1)-(vP(k,lrb)+vG(k,lrb)+vH(k,lrb));
   upper(lrb)=ones(1,length(lrb));
   lower(lrb)=uG(k,lrb).*uH(k,lrb).*abs(abs(1+mb.*exp(i*phi))-Rb);
   upper(lrc)=2*uP(k,lrc); phi=vP(k,lrc)+vG(k,lrc)+vH(k,lrc);
   lower(lrc)=abs(2*uP(k,lrc).*uG(k,lrc).*uH(k,lrc).*cos(phi)+1);
   cl(k,:)=F.*(upper./lower);
  end
 else
  for p=1:rmp, for g=1:rmg, for h=1:rmh,
   if any(R(p,:)==0),
    lr=find(R(p,:)==0);
    upper(1,lr)=P(p,lr);
    lower(1,lr)=1+P(p,lr).*G(g,lr).*H(h,lr);
   end
   if any(R(p,:)~=1 & R(p,:)~=0),
    lr=find(R(p,:)~=1 & R(p,:)~=0);
    m=uP(p,lr).*uG(g,lr).*uH(h,lr); mb=1 ./(m.*abs(1-R(p,lr).^2));
    Rb=R(p,lr).*mb; phi=pi*(R(p,lr)>1)-(vP(p,lr)+vG(g,lr)+vH(h,lr));
    upper(1,lr)=ones(1,length(lr));
    lower(1,lr)=uG(g,lr).*uH(h,lr).*abs(abs(1+mb.*exp(i*phi))-Rb);
   end
   if any(R(p,:)==1),
    lr=find(R(p,:)==1);
    upper(1,lr)=2*uP(p,lr); phi=vP(p,lr)+vG(g,lr)+vH(h,lr);
    lower(1,lr)=abs(2*uP(p,lr).*uG(g,lr).*uH(h,lr).*cos(phi)+1);
   end
   cl(pgh,:)=F.*(upper./lower); pgh=pgh+1;
  end; end; end
 end
elseif prob==4,
 if q,
  for k=1:rmp,
   lra=find(R(k,:)==0); lrb=find(R(k,:)~=0);
   upper(lra)=G(k,lra);
   lower(lra)=1+P(k,lra).*G(k,lra).*H(k,lra);
   upper(lrb)=uG(k,lrb);
   lower(lrb)=abs(abs(1+P(k,lrb).*G(k,lrb).*H(k,lrb))-...
              abs(P(k,lrb).*G(k,lrb).*H(k,lrb)).*R(k,lrb));
   cl(k,:)=F.*(upper./lower);
  end
 else
  for p=1:rmp, for g=1:rmg, for h=1:rmh,
   if any(R(p,:)==0),
    lr=find(R(p,:)==0);
    upper(1,lr)=G(g,lr); lower(1,lr)=1+P(p,lr).*G(g,lr).*H(h,lr);
   end
   if any(R(p,:)~=0),
    lr=find(R(p,:)~=0);
    upper(1,lr)=uG(g,lr);
    lower(1,lr)=abs(abs(1+P(p,lr).*G(g,lr).*H(h,lr))-...
          abs(P(p,lr).*G(g,lr).*H(h,lr)).*R(p,lr));
   end
   cl(pgh,:)=F.*(upper./lower); pgh=pgh+1;
  end; end; end
 end
elseif prob==5,
 if q,
  for k=1:rmp,
   lra=find(R(k,:)==0); lrb=find(R(k,:)~=0);
   upper(lra)=G(k,lra).*H(k,lra);
   lower(lra)=1+P(k,lra).*G(k,lra).*H(k,lra);
   upper(lrb)=uG(k,lrb).*uH(k,lrb);
   lower(lrb)=abs(abs(1+P(k,lrb).*G(k,lrb).*H(k,lrb))-...
              abs(P(k,lrb).*G(k,lrb).*H(k,lrb)).*R(k,lrb));
   cl(k,:)=F.*(upper./lower);
  end
 else
  for p=1:rmp, for g=1:rmg, for h=1:rmh,
   if any(R(p,:)==0),
    lr=find(R(p,:)==0);
    upper(1,lr)=G(g,lr).*H(h,lr);
    lower(1,lr)=1+P(p,lr).*G(g,lr).*H(h,lr);
   end
   if any(R(p,:)~=0),
    lr=find(R(p,:)~=0);
    upper(1,lr)=uG(g,lr).*uH(h,lr);
    lower(1,lr)=abs(abs(1+P(p,lr).*G(g,lr).*H(h,lr))-...
          abs(P(p,lr).*G(g,lr).*H(h,lr)).*R(p,lr));
   end
   cl(pgh,:)=F.*(upper./lower); pgh=pgh+1;
  end; end; end
 end
elseif prob==6,
 if q,
  for k=1:rmp,
   lra=find(R(k,:)==0); lrb=find(R(k,:)~=0 & R(k,:)~=1); lrc=find(R(k,:)==1);
   upper(lra)=P(k,lra).*G(k,lra);
   lower(lra)=1+P(k,lra).*G(k,lra).*H(k,lra);
   m=uP(k,lrb).*uG(k,lrb).*uH(k,lrb); mb=1 ./(m.*abs(1-R(k,lrb).^2));
   Rb=R(k,lrb).*mb; phi=pi*(R(k,lrb)>1)-(vP(k,lrb)+vG(k,lrb)+vH(k,lrb));
   upper(lrb)=ones(1,length(lrb));
   lower(lrb)=uH(k,lrb).*abs(abs(1+mb.*exp(i*phi))-Rb);
   upper(lrc)=2*uP(k,lrc).*uG(k,lrc); phi=vP(k,lrc)+vG(k,lrc)+vH(k,lrc);
   lower(lrc)=abs(2*uP(k,lrc).*uG(k,lrc).*uH(k,lrc).*cos(phi)+1);
   cl(k,:)=F.*(upper./lower);
  end
 else
  for p=1:rmp, for g=1:rmg, for h=1:rmh,
   if any(R(p,:)==0),
    lr=find(R(p,:)==0);
    upper(1,lr)=P(p,lr).*G(g,lr);
    lower(1,lr)=1+P(p,lr).*G(g,lr).*H(h,lr);
   end
   if any(R(p,:)~=1 & R(p,:)~=0),
    lr=find(R(p,:)~=1 & R(p,:)~=0);
    m=uP(p,lr).*uG(g,lr).*uH(h,lr); mb=1 ./(m.*abs(1-R(p,lr).^2));
    Rb=R(p,lr).*mb; phi=pi*(R(p,lr)>1)-(vP(p,lr)+vG(g,lr)+vH(h,lr));
    upper(1,lr)=ones(1,length(lr));
    lower(1,lr)=uH(h,lr).*abs(abs(1+mb.*exp(i*phi))-Rb);
   end
   if any(R(p,:)==1),
    upper(1,lr)=2*uP(p,lr).*uG(g,lr); phi=vP(p,lr)+vG(g,lr)+vH(h,lr);
    lower(1,lr)=abs(2*uP(p,lr).*uG(g,lr).*uH(h,lr).*cos(phi)+1);
   end
   cl(pgh,:)=F.*(upper./lower); pgh=pgh+1;
  end; end; end
 end
elseif prob==7,
 if q,
  for k=1:rmp,
   upper=P(k,:).*G(k,:); lower=1+P(k,:).*G(k,:).*H(k,:);
   cl(k,:)=F.*(upper./lower);
  end
 else
  for p=1:rmp, for g=1:rmg, for h=1:rmh,
   upper=P(p,:).*G(g,:); lower=1+P(p,:).*G(g,:).*H(h,:);
   cl(pgh,:)=F.*(upper./lower); pgh=pgh+1;
  end; end; end;
 end
elseif prob==8,
 if q,
  for k=1:rmp,
   lra=find(R(k,:)==0); lrb=find(R(k,:)~=0);
   upper(lra)=H(k,lra);
   lower(lra)=1+P(k,lra).*G(k,lra).*H(k,lra);
   upper(lrb)=uH(k,lrb);
   lower(lrb)=abs(abs(1+P(k,lrb).*G(k,lrb).*H(k,lrb))-...
              abs(P(k,lrb).*G(k,lrb).*H(k,lrb)).*R(k,lrb));
   cl(k,:)=F.*(upper./lower);
  end
 else
  for p=1:rmp, for g=1:rmg, for h=1:rmh,
   if any(R(p,:)==0),
    lr=find(R(p,:)==0);
    upper(1,lr)=H(h,lr);
    lower(1,lr)=1+P(p,lr).*G(g,lr).*H(h,lr);
   end
   if any(R(p,:)~=0),
    lr=find(R(p,:)~=0);
    upper(1,lr)=uH(h,lr);
    lower(1,lr)=abs(abs(1+P(p,lr).*G(g,lr).*H(h,lr))-...
          abs(P(p,lr).*G(g,lr).*H(h,lr)).*R(p,lr));
   end
   cl(pgh,:)=F.*(upper./lower); pgh=pgh+1;
  end; end; end
 end
else
 if q,
  for k=1:rmp,
   lra=find(R(k,:)==0); lrb=find(R(k,:)~=0 & R(k,:)~=1); lrc=find(R(k,:)==1);
   upper(lra)=P(k,lra).*H(k,lra);
   lower(lra)=1+P(k,lra).*G(k,lra).*H(k,lra);
   m=uP(k,lrb).*uG(k,lrb).*uH(k,lrb); mb=1 ./(m.*abs(1-R(k,lrb).^2));
   Rb=R(k,lrb).*mb; phi=pi*(R(k,lrb)>1)-(vP(k,lrb)+vG(k,lrb)+vH(k,lrb));
   upper(lrb)=ones(1,length(lrb));
   lower(lrb)=uG(k,lrb).*abs(abs(1+mb.*exp(i*phi))-Rb);
   upper(lrc)=2*uP(k,lrc).*uH(k,lrc); phi=vP(k,lrc)+vG(k,lrc)+vH(k,lrc);
   lower(lrc)=abs(2*uP(k,lrc).*uG(k,lrc).*uH(k,lrc).*cos(phi)+1);
   cl(k,:)=F.*(upper./lower);
  end
 else
  for p=1:rmp, for g=1:rmg, for h=1:rmh
   if any(R(p,:)==0),
    lr=find(R(p,:)==0);
    upper(1,lr)=P(p,lr).*H(h,lr);
    lower(1,lr)=1+P(p,lr).*G(g,lr).*H(h,lr);
   end
   if any(R(p,:)~=1 & R(p,:)~=0),
    lr=find(R(p,:)~=1 & R(p,:)~=0);
    m=uP(p,lr).*uG(g,lr).*uH(h,lr); mb=1 ./(m.*abs(1-R(p,lr).^2));
    Rb=R(p,lr).*mb; phi=pi*(R(p,lr)>1)-(vP(p,lr)+vG(g,lr)+vH(h,lr));
    upper(1,lr)=ones(1,length(lr));
    lower(1,lr)=uG(g,lr).*abs(abs(1+mb.*exp(i*phi))-Rb);
   end
   if any(R(p,:)==1),
    lr=find(R(p,:)==1);
    upper(1,lr)=2*uP(p,lr).*uH(h,lr); phi=vP(p,lr)+vG(g,lr)+vH(h,lr);
    lower(1,lr)=abs(2*uP(p,lr).*uG(g,lr).*uH(h,lr).*cos(phi)+1);
   end
   cl(pgh,:)=F.*(upper./lower); pgh=pgh+1;
  end; end; end
 end
end

if length(cl) & nargout<2,
 [rcl,ccl]=size(cl);
 if rcl>1, clmx=max(abs(cl)); clmn=min(abs(cl));
 else clmx=abs(cl); clmn=abs(cl); end
 a=20*log10(W(1,:));
 if nargout==0 & ccl>1,
  scrn=get(0,'screensize');sctlt=['Analysis: Problem Type ',int2str(prob)];
  f1=figure('name',sctlt,'numbertitle','off','units','norm','position',pos,...
         'vis','off','tag','qft');
  axs=[min(w),max(w)];
  if prob==7,
   axs=[axs,min(20*log10([clmn(:);W(:)]))-5,max(20*log10([clmx(:);W(:)]))+5];
  elseif rW==1,
   axs=[axs,min(20*log10([clmx(:);W(:)]))-5,max(20*log10([clmx(:);W(:)]))+5];
  else
   axs=[axs,min(min(W-abs(cl))),max(min(W-abs(cl)))];
  end
  ax=gca;
  set(ax,'box','on','xgrid','on','ygrid','on',...
        'gridlinestyle',':',...
        'nextplot','add',...
        'xlim',axs(1:2),'ylim',axs(3:4),'xscale','log');

  if length(a)~=1,
   if prob~=7,
    if rW==1,
     semilogx(w,a,'b--',w,20*log10(clmx),'-k');
     title('Weight: --');
     ylabel('Magnitude (dB)');
    else
     ermn = min(W-abs(cl));
     semilogx(w,ermn,'-k');
     title('Maximum Error');
     ylabel('Error');
    end
   else
    b=20*log10(W(2,:));
    semilogx(w,a,'r--',w,b,'r--',w,20*log10(clmx),'b-',w,20*log10(clmn),'g-');
    title('Weight: --');
    ylabel('Magnitude (dB)');
   end
   xlabel('Frequency (rad/sec)');
  end
  set(f1,'vis','on','userdata',axs);
 else
  if nargout==1,
   if prob~=7,
    if rW==1,
     err=W(1,:)-clmx;
    else
     err=min(W-abs(cl));
    end
   else
    err=[W(1,:)-clmx;W(2,:)-clmn];
   end
  end
 end
end
