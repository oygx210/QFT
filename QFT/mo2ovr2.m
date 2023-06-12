function mo2ovr2(qflag)
% MO2OVR2 Mouse-implemented 2/2 element. (Utility Function)
%         MO2OVR2 determines the 2/2 value to be added from various
%         mouse motions within the IDEs.

% Author: Craig Borghesani
% 3/25/94
% Copyright (c) 2003, Terasoft, Inc.


f=gcf;
bthan=get(f,'userdata');
infmat=get(bthan(16),'userdata');
nom=get(bthan(2),'userdata');
conts = get(bthan(31),'userdata');
if isempty(conts{1}),
   cont=get(bthan(19),'userdata');
   lomat = get(bthan(20),'userdata');
   if isempty(cont),
      cont = get(bthan(3),'userdata');
      lomat = get(bthan(1),'userdata');
      set(bthan(19),'userdata', cont);
      set(bthan(20),'userdata', lomat);
   end
   cont2 = [];

else
   cont = conts{1};
   cont2 = conts{2};
   lomat = get(bthan(20),'userdata');

end
hint_bar = get(bthan(36),'userdata');
QFTToolData = getappdata(f,'QFTToolData');
ElementsListbox = QFTToolData.Elements(17);
ElementsPopup = infmat(16,4);

w=lomat(1,:);
T=get(bthan(13),'userdata');
delay=infmat(10,1);
if infmat(9,1)==3, axs=infmat(2,3:4);
else axs=infmat(1,:); end
am = infmat(24,1);
ap = infmat(24,2);
txt1 = infmat(30,1);
txt2 = infmat(30,2);

if qflag==0,
 if ~length(cont2),
  obj=get(f,'currentobject');
  p_obj=get(obj,'parent');
  if strcmp(get(p_obj, 'type'), 'root'),
   a=get(f,'currentaxes');
  elseif p_obj==f,
   a=obj;
  else
   a=p_obj;
  end
  pt=get(a,'currentpoint');
  if infmat(9,1)==1,
   kw=qfindfrq(pt(1,1),pt(1,2),lomat,a);
   if length(kw),
    set(infmat(6,1),'xdata',qfixfase(lomat(2,:),axs,kw),...
                    'ydata',20*log10(abs(lomat(2,kw))),...
                    'vis','on');
   end
  elseif a==am | a==ap,
   if a==am,
    kw=qfindfrq(pt(1,1),pt(1,2),lomat,am);
   else
    kw=qfindfrq(pt(1,1),pt(1,2),lomat(2,:),axs,w,2);
   end
   if length(kw),
    set(infmat(6,1),'xdata',w(kw),'ydata',20*log10(abs(lomat(2,kw))),'vis','on');
    set(infmat(6,2),'xdata',w(kw),'ydata',qfixfase(lomat(2,:),axs,kw),'vis','on');
   end
  end
  if length(kw),
   set(f,'windowbuttonmotionfcn','mo2ovr2(1)',...
         'windowbuttonupfcn','mo2ovr2(2)');
   lomat(3,:)=lomat(2,:);
   infmat(3,3)=kw;
   infmat(4,3)=lomat(2,kw);
   set(hint_bar,'string','Now, drag the marker to ADD the 2/2 element based on mag and phase differences');
   set(bthan(20),'userdata',lomat);
   set(bthan(16),'userdata',infmat);
  end
 else
  set(f,'windowbuttonmotionfcn','mo2ovr2(1)',...
        'windowbuttonupfcn','mo2ovr2(2)');
  set(hint_bar,'string','Now, drag the marker to EDIT the 2/2 element based on mag and phase differences');
 end
elseif qflag==1,
 kw=infmat(3,3);
 wval=w(kw);
 p=infmat(4,3);
 go_for_it=0;
 obj=get(f,'currentobject');
 p_obj=get(obj,'parent');
 if strcmp(get(p_obj, 'type'), 'root'),
  a=get(f,'currentaxes');
 elseif p_obj==f,
  a=obj;
 else
  a=p_obj;
 end
 pt=get(a,'currentpoint');
 if infmat(9,1)==1,
  delmag=10^(pt(1,2)/20)/abs(p);
  delph=pt(1,1)-qfixfase(lomat(3,:),axs,kw);
  axislims = infmat(1,:);
 else
  if a==ap,
   delph=pt(1,2)-qfixfase(lomat(2,:),axs,kw);
   delmag = infmat(24,4);
   axislims = infmat(2,:);
  elseif a==am,
   delmag=10^(pt(1,2)/20)/abs(p);
   delph = infmat(24,3);
   axislims = infmat(1,:);
  end
 end
 if abs(delph)<176 & delph~=0,
  [contnew,state]=seconsec(delph,delmag,wval,T);
  if any(state==[1,2]), go_for_it=1; end
 else go_for_it = 2; end
 if any(pt(1,1:2)<axislims([1,3]) | pt(1,1:2)>axislims([2,4])),
  go_for_it = 3;
 end
 if ~go_for_it,
  if ~length(cont2),  % adding with the mouse
   if length(nom),
% if there are more than one second order computed, find smallest zeta
    ztaloc=find(contnew(:,4)==3 | contnew(:,4)==4);
    zta=contnew(ztaloc,1);
    wn=contnew(ztaloc,2);
    [zta,loc]=min(zta);
    wn=wn(loc);
    wlt=qfrqenh(wn,zta,w,T);
    if infmat(9,1)==1,
     nlot=squeeze(freqresp(nom,wlt)).';
     clot=qcntbode(cont,wlt,T);
%     clot=qcntbode(cont,wlt,T).*exp(-i*wlt*delay);
     cpnew=qcntbode(contnew,wlt,T);
     lot=nlot.*clot;
    else
     lot=qcntbode(cont,wlt,T).*exp(-i*wlt*delay);
     cpnew=qcntbode(contnew,wlt,T);
    end
    lomat=[wlt;lot.*cpnew;lot];
   else
    cpnv=qcntbode(contnew,w,T);
    lomat(2,:)=lomat(2,:).*cpnv;
   end
  else  % iterating with the mouse or continuation of adding
   cp=qcntbode(cont2,w,T);
   cpnv=qcntbode(contnew,w,T);
   lomat(2,:)=lomat(2,:).*(cpnv./cp);
  end
  kw=find(lomat(1,:)>=w(kw)); kw=kw(1);
  if infmat(9,1)==1,
   set(infmat(6,1),'xdata',qfixfase(lomat(2,:),axs,kw),...
                   'ydata',20*log10(abs(lomat(2,kw))));
  else
   set(infmat(6,1),'xdata',lomat(1,kw),'ydata',20*log10(abs(lomat(2,kw))));
   set(infmat(6,2),'xdata',lomat(1,kw),'ydata',qfixfase(lomat(2,:),axs,kw));
  end

% add 2/2 element to controller
  if any(infmat(9,1)==[1,3]),
   if length(contnew),
    cont(1,1) = cont(1,1)*contnew(1,1);
    if T > 0,
     cont(3,1) = cont(3,1)+contnew(3,1);
     contnew(1:3,:) = [];
    else
     contnew(1:2,:) = [];
    end
    contv = [cont; contnew];
   else
      contv = cont;
   end
   ElementLocation = size(contv,1);
  end

  infmat(3,3)=kw;
  infmat(4,3)=lomat(3,kw);
  infmat(24,3:4)=[delph,delmag];
  set(bthan(19),'userdata',contv);
  set(bthan(20),'userdata',lomat);
  set(bthan(16),'userdata',infmat);
  set(bthan(31),'userdata',{cont, contnew});
  if infmat(9,1)==1, qnicplt(f); setptr(f,'fleur');
  else mgphplot(f); end
  [ControllerString,ListboxInfo,ListboxValue] = cntstr(f,contv,ElementLocation);
  set(ElementsListbox, 'string', ControllerString,...
                       'userdata', ListboxInfo,...
                       'value',ListboxValue);
  set(hint_bar,'string','2/2 element being implemented ...','fore','k');

%  drawnow;
 elseif go_for_it == 1,
  set(hint_bar,'string','Implementation not possible.','fore','r');
  setptr(f,'forbidden');
 elseif go_for_it == 2,
  set(hint_bar,'string','Phase difference must be > 0 and < 176 degrees.','fore','r');
  setptr(f,'forbidden');
 elseif go_for_it == 3,
  set(hint_bar,'string','Must remain within axis limits.','fore','r');
  setptr(f,'forbidden');
 end

elseif qflag==2,

 set(f,'windowbuttonmotionfcn','qmouse(''Floating'')','windowbuttondownfcn','',...
       'windowbuttonupfcn','');
 set(infmat(6,1),'userdata','mo2ovr2(0)');
 set(hint_bar,'string','To continue with EDIT of the 2/2 element re-select the marker');
 set(infmat(13,2),'callback','qelmts(-1,0,''Done'')','enable','on');
 set(infmat(13,3),'callback','qelmts(-1,0,''Cancel'')','enable','on');
 qelmtlistbox('Iterate');
 drawnow;

end
