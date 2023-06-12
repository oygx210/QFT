function mofirst(qflag,qflag2)
% MOFIRST Mouse-implemented first order term. (Utility Function)
%         MOFIRST determines the first order value to be added from
%         various mouse motions within the IDEs.

% Author: Craig Borghesani
% 9/5/93
% Copyright (c) 2003, Terasoft, Inc.


f=gcf;
str=['pole';'zero'];
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

w=lomat(1,:);
q=1; r=2;
if infmat(9,1)==2, q=[1;1]; r=2:3; end
lo=lomat(r,:);
T=get(bthan(13),'userdata');
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
  elseif p_obj == f,
   a=obj;
  else
   a=p_obj;
  end

  pt=get(a,'currentpoint');
  if infmat(9,1)==1,
   kw=qfindfrq(pt(1,1),pt(1,2),lomat,a);
   if length(kw),
    set(infmat(6,1),'xdata',qfixfase(lo,axs,kw),...
                    'ydata',20*log10(abs(lo(kw))),...
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
   set(f,'windowbuttonmotionfcn',['mofirst(',int2str(lcont+1),',1)'],...
         'windowbuttonupfcn','mofirst(0,2)');
   infmat(3,3)=kw;
   if a==am,
    set(hint_bar,'string','Now, drag the marker to ADD the First Order element based on magnitude difference');
   else
    set(hint_bar,'string','Now, drag the marker to ADD the First Order element based on phase difference');
   end
  end
  set(bthan(16),'userdata',infmat);
 else
  set(f,'windowbuttonmotionfcn',['mofirst(',int2str(lcont),',1)'],...
        'windowbuttonupfcn','mofirst(0,2)');
  set(hint_bar,'string','Now, drag the marker to EDIT the First Order element');
 end

elseif qflag2==1,
 kw=infmat(3,3);
 wval=w(kw);
 p=infmat(4,3);
 am = infmat(24,1);
 ap = infmat(24,2);
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
 if infmat(9,1)==1 | a==ap,
  pt=get(a,'currentpoint');
  if infmat(9,1)==1,
   axislims = infmat(1,:);
  elseif a==ap,
   axislims = infmat(2,:);
  end
  val=0;
  delval=180/pi*atan2(imag(p),real(p))-(qfixfase(lo,axs,kw)-pt(1,1+(a==ap)));

  if abs(delval)<88 & delval~=0,
   if T > 0,
    phi=pi-(atan2(sin(wval*T),cos(wval*T))+pi/180*abs(delval));
    if lcont~=qflag,
     rt=cos(wval*T)+sin(wval*T)/tan(phi);
    else
     rt=sign(cont(qflag,1))*cos(wval*T)+sin(wval*T)/tan(phi);
    end
   else
    if lcont~=qflag,
     rt=wval/tan(abs(delval*pi/180));
    else
     rt=sign(cont(qflag,1))*wval/tan(abs(delval*pi/180));
    end
   end
   if T > 0 & rt>1, go_for_it = 1; end
  else go_for_it = 2; end
 elseif infmat(9,1)==2 | a==am,
  pt=get(a,'currentpoint');
  axislims = infmat(1,:);
  val=1;
  if p~=0,
   delval=(abs(p)*10^(pt(1,2)/20))/abs(lo(1,kw));
  else
   delval=10^(pt(1,2)/20)/abs(lo(1,kw));
  end
  if T > 0,
   if delval<1,
    delmax=1/sqrt(2*(1-cos(wval*T))); delmin=1/sqrt(2*(1+cos(wval*T)));
   else
    delmin=sqrt(2*(1-cos(wval*T))); delmax=sqrt(2*(1+cos(wval*T)));
   end
   if (delval<1 & delval>=delmin & delval<=delmax) | ...
      (delval>1 & delval>=delmin & delval<=delmax),
    if delval<1, rts=roots([1-1/delval^2 -2*(cos(wval*T)-1/delval^2) 1-1/delval^2]);
    else rts=roots([1-delval^2 -2*(cos(wval*T)-delval^2) 1-delval^2]); end
    rt=rts(find(abs(rts)<0.999999));
   else go_for_it = 1; end
  else
   if delval<1, rt=sqrt(wval^2/(1/delval^2-1));
   else rt=sqrt(wval^2/(delval^2-1)); end
  end
 end
 if any(pt(1,1:2)<axislims([1,3]) | pt(1,1:2)>axislims([2,4])),
  go_for_it = 3;
 end
 if ~go_for_it,
  if lcont~=qflag,
   cp=ones(size(lomat(r,:)));
   cont=[cont;0,0,0,0];
  else
   cp=rlroot(cont(qflag,1),w,[(cont(qflag,4)==2)-(cont(qflag,4)==1),T]);
  end
  if T > 0, rt = -log(rt)/T; end
  cont(qflag,1:4)=[rt,NaN,NaN,(delval<val)+2*(delval>val)];
  cpnv=rlroot(rt,w,[(cont(qflag,4)==2)-(cont(qflag,4)==1),T]);
  lo=lo.*(cpnv(q,:)./cp(q,:));
  lomat(r,:)=lo;
  ElementLocation = qflag;
  set(ElementsPopup,'value',cont(qflag,4)+1);
  if cont(qflag,4) == 1,
   set(infmat(19,1),'string','pole');
  else
   set(infmat(19,1),'string','zero');
  end
  set(infmat(16,1),'string',num2str(cont(ElementLocation,1)));

  if infmat(9,1)==1,
   set(infmat(6,1),'xdata',qfixfase(lo(1,:),axs,kw),...
                   'ydata',20*log10(abs(lo(1,kw))));
  elseif infmat(9,1)==2,
   set(infmat(6,1),'xdata',w(kw),'ydata',20*log10(abs(lo(1,kw))));
  elseif infmat(9,1)==3,
   set(infmat(6,1),'xdata',w(kw),'ydata',20*log10(abs(lo(1,kw))));
   set(infmat(6,2),'xdata',w(kw),'ydata',qfixfase(lo(1,:),axs,kw));
  end
  infmat(4,3)=cpnv(kw);
  set(bthan(20),'userdata',lomat);
  set(bthan(19),'userdata',cont);
  set(bthan(16),'userdata',infmat);
  if infmat(9,1)==1, qnicplt(f); setptr(f,'lrdrag');
  elseif infmat(9,1)==2, qmagplt(f); setptr(f,'uddrag');
  elseif infmat(9,1)==3, mgphplot(f);
  end
  [ControllerString,ListboxInfo,ListboxValue] = cntstr(f,cont,ElementLocation);
  set(ElementsListbox, 'string', ControllerString,...
                       'userdata', ListboxInfo,...
                       'value',ListboxValue);

  set(hint_bar,'string','First Order element being implemented ...','fore','k');
%  drawnow;
 elseif go_for_it == 1,
  set(hint_bar,'string','Implementation not possible','fore','r');
  setptr(f,'forbidden');
 elseif go_for_it == 2,
  set(hint_bar,'string','Phase difference must be > 0 and < 88 degrees','fore','r');
  setptr(f,'forbidden');
 elseif go_for_it == 3,
  set(hint_bar,'string','Must remain within axis limits','fore','r');
  setptr(f,'forbidden');
 end

elseif qflag2==2,

 set(f,'windowbuttonmotionfcn','qmouse(''Floating'')','windowbuttondownfcn','',...
       'windowbuttonupfcn','');
 set(infmat(6,1),'userdata',['mofirst(',int2str(lcont),',0)']);
 set(hint_bar,'string','To continue with EDIT of the First Order element re-select the marker');
 set(infmat(13,2),'callback','qelmts(-1,0,''Done'')','enable','on');
 set(infmat(13,3),'callback','qelmts(-1,0,''Cancel'')','enable','on');
 qelmtlistbox('Iterate');
 drawnow;

end
