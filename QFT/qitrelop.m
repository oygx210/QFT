function qitrelop(flag,flag2)
% QITRELOP Iterate operations. (Utility Function)
%          QITRELOP executes the iteration slider presses.

% Author: Craig Borghesani
% 9/6/93
% Copyright (c) 2003, Terasoft, Inc.


f2=gcf;
f=get(f2,'userdata');
bthan=get(f,'userdata');
infmat=get(bthan(16),'userdata');

lomat=get(bthan(20),'userdata');
T=get(bthan(13),'userdata');
cont=get(bthan(19),'userdata');
set(infmat(21,3),'enable','on');

flag4=infmat(9,1);
if any(flag4==[1 3]),
 q=1; s=0;
elseif flag4==2,
 q=[1;1]; s=1;
end

butn_sel = get(infmat(21,1),'userdata');
cur_str = get(butn_sel,'string');
brac = find(cur_str=='[');
new_str = cur_str(1:brac);

if flag==0, % Gain
 oldk=cont(flag2,1);
 gain_val = get(infmat(17,1),'value');
 cont(flag2,1)=sign(cont(flag2,1))*gain_val;
 rtnv=ones(1,length(lomat(1,:)))*cont(flag2,1);
 lomat(2:2+s,:)=lomat(2:2+s,:).*(cont(flag2,1)/oldk);
 slider_min=get(infmat(17,1),'min');
 slider_max=get(infmat(17,1),'max');
 if gain_val <= slider_min | gain_val >= slider_max,
  set(infmat(17,1),'min',10^((20*log10(gain_val)-5)/20),...
                   'max',10^((20*log10(gain_val)+5)/20));
 end
 v1=num2str(cont(flag2,1),4);
 set(butn_sel,'string',[new_str,v1,']']);
elseif any(flag==[1 2]), % First order
 oldp=cont(flag2,1);
 newp=get(infmat(17,1),'value');
 if imag(oldp)~=0,
  sign_imag = sign(imag(oldp));
  newp = newp + sign_imag*pi/T*i;
 end
 rt=rlroot(oldp,lomat(1,:),[(flag==2)-(flag==1) T]);
 rtnv=rlroot(newp,lomat(1,:),[(flag==2)-(flag==1) T]);
 lomat(2:2+s,:)=lomat(2:2+s,:).*(rtnv(q,:)./rt(q,:));

 slider_min=get(infmat(17,1),'min');
 slider_max=get(infmat(17,1),'max');
 if real(newp) <= slider_min | real(newp) >= slider_max,
  set(infmat(17,1),'min',real(newp)*0.5,...
                   'max',real(newp)*1.5);
 end
 cont(flag2,1) = newp;
 v1=num2str(cont(flag2,1),4);
 set(butn_sel,'string',[new_str,v1,']']);

elseif any(flag==[3 4]), % Second order
 oldz=cont(flag2,1); oldw=cont(flag2,2);
 cont(flag2,1)=get(infmat(17,2),'value');
 cont(flag2,2)=get(infmat(17,1),'value');
 rt=cproot(oldz,oldw,lomat(1,:),[(flag==4)-(flag==3) T]);
 rtnv=cproot(cont(flag2,1),cont(flag2,2),lomat(1,:),[(flag==4)-(flag==3) T]);
 lomat(2:2+s,:)=lomat(2:2+s,:).*(rtnv(q,:)./rt(q,:));
 rtnv=lomat(2,:)./rtnv;

 slider_min=get(infmat(17,2),'min');
 slider_max=get(infmat(17,2),'max');
 if cont(flag2,1) <= slider_min | cont(flag2,1) >= slider_max,
  set(infmat(17,2),'min',cont(flag2,1)*0.5,...
                   'max',cont(flag2,1)*1.5);
 end

 slider_min=get(infmat(17,1),'min');
 slider_max=get(infmat(17,1),'max');
 if cont(flag2,2) <= slider_min | cont(flag2,2) >= slider_max,
  if cont(flag2,2)*1.5 <= lomat(1,length(lomat(1,:))) & ...
     cont(flag2,2)*0.5 >= lomat(1,1),
   set(infmat(17,1),'min',cont(flag2,2)*0.5,...
                    'max',cont(flag2,2)*1.5);
  elseif cont(flag2,2)*1.5 > lomat(1,length(lomat(1,:))),
   set(infmat(17,1),'min',cont(flag2,2)*0.5,...
                    'max',lomat(1,length(lomat(1,:))));
  elseif cont(flag2,2)*0.5 < lomat(1,1),
   set(infmat(17,1),'min',lomat(1,1),...
                    'max',cont(flag2,2)*1.5);

  end
 end
 v1=num2str(cont(flag2,1),4);
 v2=num2str(cont(flag2,2),4);
 set(butn_sel,'string',[new_str,v1,', ',v2,']']);

elseif flag==5, % Lead/Lag
 oldp=cont(flag2,1); oldw=cont(flag2,2);
 rt=ldlgcplx(oldp,oldw,lomat(1,:),T);
 cont(flag2,1)=get(infmat(17,2),'value');
 cont(flag2,2)=get(infmat(17,1),'value');
 rtnv=ldlgcplx(cont(flag2,1),cont(flag2,2),lomat(1,:),T);
 lomat(2,:)=lomat(2,:).*(rtnv./rt);

 slider_min=get(infmat(17,1),'min');
 slider_max=get(infmat(17,1),'max');
 if cont(flag2,2) <= slider_min | cont(flag2,2) >= slider_max,
  if cont(flag2,2)*1.5 <= lomat(1,length(lomat(1,:))) & ...
     cont(flag2,2)*0.5 >= lomat(1,1),
   set(infmat(17,1),'min',cont(flag2,2)*0.5,...
                    'max',cont(flag2,2)*1.5);
  elseif cont(flag2,2)*1.5 > lomat(1,length(lomat(1,:))),
   set(infmat(17,1),'min',cont(flag2,2)*0.5,...
                    'max',lomat(1,length(lomat(1,:))));
  elseif cont(flag2,2)*0.5 < lomat(1,1),
   set(infmat(17,1),'min',lomat(1,1),...
                    'max',cont(flag2,2)*1.5);

  end
 end
 v1=num2str(cont(flag2,1),4);
 v2=num2str(cont(flag2,2),4);
 set(butn_sel,'string',[new_str,v1,', ',v2,']']);

elseif flag==6, % Notch
 zta1=cont(flag2,1); zta2=cont(flag2,2);
 rt=ntchcplx(zta1,zta2,cont(flag2,3),lomat(1,:),T);
 cont(flag2,1)=get(infmat(17,3),'value');
 cont(flag2,2)=get(infmat(17,2),'value');
 cont(flag2,3)=get(infmat(17,1),'value');
 rtnv=ntchcplx(cont(flag2,1),cont(flag2,2),cont(flag2,3),lomat(1,:),T);
 lomat(2:2+s,:)=lomat(2:2+s,:).*(rtnv(q,:)./rt(q,:));

 slider_min=get(infmat(17,3),'min');
 slider_max=get(infmat(17,3),'max');
 if cont(flag2,1) <= slider_min | cont(flag2,1) >= slider_max,
  set(infmat(17,3),'min',cont(flag2,1)*0.5,...
                   'max',cont(flag2,1)*1.5);
 end

 slider_min=get(infmat(17,2),'min');
 slider_max=get(infmat(17,2),'max');
 if cont(flag2,2) <= slider_min | cont(flag2,2) >= slider_max,
  set(infmat(17,2),'min',cont(flag2,2)*0.5,...
                   'max',cont(flag2,2)*1.5);
 end

 slider_min=get(infmat(17,1),'min');
 slider_max=get(infmat(17,1),'max');
 if cont(flag2,3) <= slider_min | cont(flag2,3) >= slider_max,
  if cont(flag2,3)*1.5 <= lomat(1,length(lomat(1,:))) & ...
     cont(flag2,3)*0.5 >= lomat(1,1),
   set(infmat(17,1),'min',cont(flag2,3)*0.5,...
                    'max',cont(flag2,3)*1.5);
  elseif cont(flag2,3)*1.5 > lomat(1,length(lomat(1,:))),
   set(infmat(17,1),'min',cont(flag2,3)*0.5,...
                    'max',lomat(1,length(lomat(1,:))));
  elseif cont(flag2,3)*0.5 < lomat(1,1),
   set(infmat(17,1),'min',lomat(1,1),...
                    'max',cont(flag2,3)*1.5);

  end
 end
 v1=num2str(cont(flag2,1),4);
 v2=num2str(cont(flag2,2),4);
 v3=num2str(cont(flag2,3),4);
 set(butn_sel,'string',[new_str,v1,', ',v2,', ',v3,']']);

end
set(bthan(19),'userdata',cont);
set(bthan(20),'userdata',lomat);

if flag4==1, qnicplt(f);
elseif flag4==2, qmagplt(f);
elseif flag4==3, mgphplot(f);
end
