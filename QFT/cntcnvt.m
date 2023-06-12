function cntcnvt(flag)
% CNTCNVT Controller convertor 1. (Utility Function)
%         CNTCNVT converts pole/zero pairs into lead/lag terms.

% Author: Craig Borghesani
% Date: 11/1/93
% Revised: 2/17/96 11:25 PM V1.1 updates
% Copyright (c) 2003, Terasoft, Inc.


if flag,
 f2=gcf;
 f=get(f2,'userdata');
else
 qclswin(0);
 f=gcf;
end

bthan=get(f,'userdata');
infmat=get(bthan(16),'userdata');
proc_str=[]; dec=[];

if infmat(25,2)>1, proc_str=[' (',int2str(infmat(25,2)),')']; end

fig_color=[0.5,0.5,0.5];
if any(flag==[0,2,4]),
 if flag==2,
  cont=get(bthan(19),'userdata');
 else
  cont=get(bthan(3),'userdata');
  set(bthan(28),'userdata',cont);
  set(bthan(27),'userdata',[]);
  set(infmat(8,1),'enable','off');
 end
 set(infmat(8,1),'enable','off');

 T=get(bthan(13),'userdata');
 lp=length(find(cont(:,4)==1));
 lz=length(find(cont(:,4)==2));
 if sum([lp,lz])~=0 & lp>=1 & lz>=1,
  set(bthan(19),'userdata',cont);
  cntdisp(f,cont,1);
 else
  errordlg('No first order pairs present to convert','Message','on');
 end
elseif flag==1,
 t=get(bthan(23),'userdata');
 cont=get(bthan(19),'userdata');
 T=get(bthan(13),'userdata');

 for tv=1:length(t),
  vv_sty=get(t(tv),'style');
  vv_val=get(t(tv),'value');
  vv_vec(tv)=vv_val;
  if strcmp(vv_sty,'checkbox'),
   if vv_val,
    dec=[dec vv_val];
   end
  end
 end
 if length(dec),
  contzp=cont(dec,:);
  cont(dec,:)=[];
  [contldlg,msg]=zp2ldlg(contzp,T);
  if length(contldlg),
   contnew=[cont;contldlg];
   set(bthan(19),'userdata',contnew);
   cntdisp(f,contnew,3);
  else
   if msg==1,
    errordlg('More ZEROS than POLES selected','Message','on');
   elseif msg==2,
    errordlg('More POLES than ZEROS selected','Message','on');
   elseif msg==3,
    errordlg('Unstable ZERO selected','Message','on');
   elseif msg==4,
    errordlg('Unstable POLE selected','Message','on');
   end
  end
 end
elseif flag==3,
 contnew=get(bthan(19),'userdata');
 set(bthan(3),'userdata',contnew);
 set(infmat(8,1),'enable','on');
 qclswin(1);
end
