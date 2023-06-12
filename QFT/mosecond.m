function mosecond(qflag,qflag2)
% MOSECNOD Mouse-implemented second order element. (Utility Function)
%          MOSECNOD determines the second order value to be added from
%          various mouse motions within the IDEs.

% Author: Craig Borghesani
% 9/5/93
% Copyright (c) 2003, Terasoft, Inc.


f=gcf;
bthan=get(f,'userdata');
infmat=get(bthan(16),'userdata');
cont = get(bthan(19),'userdata');
lomat = get(bthan(20),'userdata');
if isempty(cont),
   cont = get(bthan(3),'userdata');
   lomat=get(bthan(1),'userdata');
   set(bthan(19),'userdata',cont);
   set(bthan(20),'userdata',lomat);
end
hint_bar = get(bthan(36),'userdata');
QFTToolData = getappdata(f,'QFTToolData');
ElementsListbox = QFTToolData.Elements(17);
ElementsPopup = infmat(16,4);
nom=get(bthan(2),'userdata');

w=lomat(1,:);
q=1; r=2;
if infmat(9,1)==2, q=[1;1]; r=2:3; end
T=get(bthan(13),'userdata');
delay=infmat(10,1);
if infmat(9,1)==3, axs=infmat(2,3:4);
else axs=infmat(1,:); end
am = infmat(24,1);
ap = infmat(24,2);
txt1 = infmat(30,1);
txt2 = infmat(30,2);
lcont = length(cont(:,1));

if qflag2==0,
 if qflag == 0,
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
  elseif infmat(9,1)==2,
   kw=qfindfrq(pt(1,1),pt(1,2),lomat,a);
   if length(kw),
    set(infmat(6,1),'xdata',w(kw),'ydata',20*log10(abs(lomat(2,kw))),'vis','on');
   end
  elseif infmat(9,1)==3,
   if a==am,
    kw=qfindfrq(pt(1,1),pt(1,2),lomat,a);
   elseif a==ap,
    kw=qfindfrq(pt(1,1),pt(1,2),lomat,a,1);
   end
   if length(kw),
    set(infmat(6,1),'xdata',w(kw),'ydata',20*log10(abs(lomat(2,kw))),'vis','on');
    set(infmat(6,2),'xdata',w(kw),'ydata',qfixfase(lomat(2,:),axs,kw),'vis','on');
   end
  end
  if length(kw),
   lomat(3+(infmat(9,1)==2),:)=lomat(2,:);
   infmat(3,3)=kw;
   infmat(4,3)=lomat(2,kw);
   set(f,'windowbuttonmotionfcn',['mosecond(',int2str(lcont+1),',1)'],...
         'windowbuttonupfcn','mosecond(0,2)');
   if a==am & infmat(9,1)==1,
    set(hint_bar,'string','Now, drag the marker to ADD the Second Order element based on mag and phase differences');
   elseif a==am,
    set(hint_bar,'string','Now, drag the marker to ADD the Second Order element based on mag difference');
   else
    set(hint_bar,'string','Now, drag the marker to ADD the Second Order element based on phase difference');
   end
   set(bthan(20),'userdata',lomat);
   set(bthan(16),'userdata',infmat);
  end
 else
  set(f,'windowbuttonmotionfcn',['mosecond(',int2str(lcont),',1)'],...
        'windowbuttonupfcn','mosecond(0,2)');
  set(hint_bar,'string','Now, drag the marker to EDIT the Second Order element');
 end
elseif qflag2==1,
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
   delph=pt(1,2)-qfixfase(lomat(3,:),axs,kw);
   axislims = infmat(2,:);
  elseif a==am,
   delmag=10^(pt(1,2)/20)/abs(p);
   axislims = infmat(1,:);
  end
 end
 if infmat(9,1)==1,
  if abs(delph)<176 & delph~=0,
   delval=delph;
   val=0;
   if T > 0, [zta,wn,state]=dsecond(delph,delmag,wval,T);
   else [zta,wn,state]=csecond(delph,delmag,wval); end
   if any(state==[1,2]), go_for_it=1; end
  else go_for_it = 2; end
 elseif infmat(9,1)==2,
  delval=delmag;
  if delmag > 1, delmag = 1/delmag; end
  [zta,wn,state]=csecond(0,delmag,wval);
  val=1;
 elseif infmat(9,1)==3,
  if a==ap,
   if abs(delph)<176 & delph~=0,
    [zta,wn,state]=csecond(delph,0,wval);
    delval=delph;
    val=0;
   else go_for_it = 2; end
  else
   delval=delmag;
   if delmag > 1, delmag = 1/delmag; end
   [zta,wn,state]=csecond(0,delmag,wval);
   val=1;
  end
 end
 if any(pt(1,1:2)<axislims([1,3]) | pt(1,1:2)>axislims([2,4])),
  go_for_it = 3;
 end
 if ~go_for_it,
  ElementLocation = qflag;
  if lcont~=qflag,  % adding with the mouse
   cont=[cont;0,0,0,0];
   cont(qflag,1:4)=[zta,wn,NaN,3*(delval<val)+4*(delval>val)];
   if length(nom) & zta<1 & infmat(9,1)~=2,
    wlt=qfrqenh(wn,zta,w,T);
    if infmat(9,1)~=3,
     nlot=squeeze(freqresp(nom,wlt)).';
%     clot=qcntbode(cont,wlt,T).*exp(-i*wlt*delay);
     clot=qcntbode(cont,wlt,T);
     lot=nlot.*clot;
     lot2=lot./qcntbode(cont(qflag,:),wlt,T);
    else
     lot=qcntbode(cont,wlt,T).*exp(-i*wlt*delay);
     lot2=lot./qcntbode(cont(qflag,:),wlt,T);
    end
    lomat=[wlt;lot;lot2];
   else
    cpnv=cproot(zta,wn,w,[(cont(qflag,4)==4)-(cont(qflag,4)==3),T]);
    lomat(r,:)=lomat(r,:).*cpnv(q,:);
   end
  else  % iterating with the mouse or continuation of adding
   cp=cproot(cont(qflag,1),cont(qflag,2),w,[(cont(qflag,4)==4)-(cont(qflag,4)==3),T]);
   cont(qflag,1:4)=[zta,wn,NaN,3*(delval<val)+4*(delval>val)];
   cpnv=cproot(zta,wn,w,[(cont(qflag,4)==4)-(cont(qflag,4)==3),T]);
   lomat(r,:)=lomat(r,:).*(cpnv(q,:)./cp(q,:));
  end
  set(ElementsPopup,'value',cont(ElementLocation,4)+1);
  set(infmat(16,1),'string',num2str(cont(ElementLocation,1)));
  set(infmat(16,2),'string',num2str(cont(ElementLocation,2)));

  kw=find(lomat(1,:)>=w(kw)); kw=kw(1);
  if infmat(9,1)==1,
   set(infmat(6,1),'xdata',qfixfase(lomat(2,:),axs,kw),...
                   'ydata',20*log10(abs(lomat(2,kw))));
  elseif infmat(9,1)==2,
   set(infmat(6,1),'xdata',w(kw),'ydata',20*log10(abs(lomat(2,kw))));
  elseif infmat(9,1)==3,
   set(infmat(6,1),'xdata',lomat(1,kw),'ydata',20*log10(abs(lomat(2,kw))));
   set(infmat(6,2),'xdata',lomat(1,kw),'ydata',qfixfase(lomat(2,:),axs,kw));
  end
  infmat(3,3)=kw;
  infmat(4,3)=lomat(3+(infmat(9,1)==2),kw);
  set(bthan(20),'userdata',lomat);
  set(bthan(19),'userdata',cont);
  set(bthan(16),'userdata',infmat);
  if infmat(9,1)==1, qnicplt(f); setptr(f,'fleur');
  elseif infmat(9,1)==2, qmagplt(f); setptr(f,'uddrag');
  elseif infmat(9,1)==3, mgphplot(f);
  end
  [ControllerString,ListboxInfo,ListboxValue] = cntstr(f,cont,ElementLocation);
  set(ElementsListbox, 'string', ControllerString,...
                       'userdata', ListboxInfo,...
                       'value',ListboxValue);
  set(hint_bar,'string','Second Order being implemented ...','fore','k');
%  drawnow;
 elseif go_for_it == 1,
  set(hint_bar,'string','Implementation not possible.','fore','r');
  setptr(f,'forbidden');
 elseif go_for_it == 2,
  set(hint_bar,'string','Phase difference must be > 0 and < 176 degrees.','fore','r');
  setptr(f,'forbidden')
 elseif go_for_it == 3,
  set(hint_bar,'string','Must remain within axis limits.','fore','r');
  setptr(f,'forbidden')
 end

elseif qflag2==2,

 set(f,'windowbuttonmotionfcn','qmouse(''Floating'')','windowbuttondownfcn','',...
       'windowbuttonupfcn','');
 set(infmat(6,1),'userdata',['mosecond(',int2str(lcont),',0)']);
 set(hint_bar,'string','To continue with EDIT of the Second Order element re-select the marker');
 set(infmat(13,2),'callback','qelmts(-1,0,''Done'')','enable','on');
 set(infmat(13,3),'callback','qelmts(-1,0,''Cancel'')','enable','on');
 qelmtlistbox('Iterate');
 drawnow;

end
