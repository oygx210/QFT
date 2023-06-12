function mocpld(qflag,qflag2)
% MOCPLD Mouse-implemented complex lead element. (Utility Function)
%        MOCPLD determines the complex lead value to be added from various
%        mouse motions within the IDEs.

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

w=lomat(1,:);
q=1; r=2;
lo=lomat(r,:);
T=get(bthan(13),'userdata');
ap = infmat(24,2);
if infmat(9,1)==3, axs=infmat(2,3:4);
else axs=infmat(1,:); end
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
    set(infmat(6,1),'xdata',qfixfase(lo,axs,kw),...
                    'ydata',20*log10(abs(lo(kw))),...
                    'vis','on');
    set(f,'windowbuttonmotionfcn',['mocpld(',int2str(lcont+1),',1)'],...
          'windowbuttonupfcn','mocpld(0,2)');
    set(hint_bar,'string','Now, drag the marker to ADD the Complex Lead/Lag element based on phase/mag difference');
    infmat(3,3)=kw;
   end
  elseif a==ap
   kw=qfindfrq(pt(1,1),pt(1,2),lomat,a,1);
   if length(kw),
    set(infmat(6,2),'xdata',w(kw),'ydata',qfixfase(lo(1,:),axs,kw),'vis','on');
    set(f,'windowbuttonmotionfcn',['mocpld(',int2str(lcont+1),',1)'],...
          'windowbuttonupfcn','mocpld(0,2)');
    infmat(3,3)=kw;
    set(hint_bar,'string','Now, drag the marker to ADD the Complex Lead/Lag element based on phase/mag difference');
   end
  end
  set(bthan(16),'userdata',infmat);
 else
  set(f,'windowbuttonmotionfcn',['mocpld(',int2str(lcont),',1)'],...
        'windowbuttonupfcn','mocpld(0,2)');
  set(hint_bar,'string','Now, drag the marker to EDIT the Complex Lead/Lag element');
 end
elseif qflag2==1,
 kw=infmat(3,3);
 wval=w(kw);
 p=infmat(4,3);
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
 if (infmat(9,1)==1 | a==ap),
  pt=get(a,'currentpoint');
  if a==ap,
   delph=pt(1,2)-qfixfase(lo,axs,kw);
   axislims = infmat(2,:);
  else
   delph=pt(1,1)-qfixfase(lo,axs,kw);
   axislims = infmat(1,:);
   if p~=0,
    delmag=abs(p)*10^(pt(1,2)/20)/abs(lomat(2,kw));
   else
    delmag=10^(pt(1,2)/20)/abs(lomat(2,kw));
   end
  end
  if lcont~=qflag,
   init_ph = 0;
  else
   init_ph = cont(qflag,1);
  end
  if any(pt(1,1:2)<axislims([1,3]) | pt(1,1:2)>axislims([2,4])),
   go_for_it = 2;
  end
  if abs(init_ph+delph)>178, go_for_it = 1; end
  if ~go_for_it,
   if lcont~=qflag,
    cp=ones(size(lomat(r,:)));
    cont=[cont;0,0,0,0];
   else
    cp=complexlead(cont(qflag,1), cont(qflag,2), cont(qflag,3), w, T);
   end
   cont(qflag,1:4)=[cont(qflag,1)+delph, 0.45, wval, 7];
%   cont(qflag,1:4)=[delph, delmag, wval, 7];
   cpnv=complexlead(cont(qflag,1), cont(qflag,2), cont(qflag,3), w, T);
   lo=lo.*(cpnv./cp);
   lomat(r,:)=lo;
   ElementLocation = qflag;
   set(ElementsPopup,'value',cont(ElementLocation,4)+2);
   set(infmat(16,1),'string',num2str(cont(ElementLocation,1)));
   set(infmat(16,2),'string',num2str(cont(ElementLocation,2)));
   set(infmat(16,3),'string',num2str(cont(ElementLocation,3)));
   if infmat(9,1)==1,
    set(infmat(6,1),'xdata',qfixfase(lo(1,:),axs,kw),...
                    'ydata',20*log10(abs(lo(1,kw))));
   else
    set(infmat(6,2),'xdata',w(kw),'ydata',qfixfase(lo(1,:),axs,kw));
   end
   infmat(4,3)=cpnv(kw);
   set(bthan(20),'userdata',lomat);
   set(bthan(19),'userdata',cont);
   set(bthan(16),'userdata',infmat);
   if infmat(9,1)==1, qnicplt(f); setptr(f,'lrdrag');
   else mgphplot(f); end
   set(hint_bar,'string','Complex Lead\Lag being implemented ...','fore','k');
   [ControllerString,ListboxInfo,ListboxValue] = cntstr(f,cont,ElementLocation);
   set(ElementsListbox, 'string', ControllerString,...
                        'userdata', ListboxInfo,...
                        'value',ListboxValue);
%   drawnow;
  elseif go_for_it == 1,
   set(hint_bar,'string','Phase difference must be > 0 and < 178 degrees.','fore','r');
   setptr(f,'forbidden');
  elseif go_for_it == 2,
   set(hint_bar,'string','Must remain within axis limits.','fore','r');
   setptr(f,'forbidden');
  end
 end

elseif qflag2==2,

 set(f,'windowbuttonmotionfcn','qmouse(''Floating'')','windowbuttondownfcn','',...
       'windowbuttonupfcn','');
 set(infmat(6,1),'userdata',['mocpld(',int2str(lcont),',0)']);
 set(hint_bar,'string','To continue with EDIT of the Complex Lead/Lag element re-select the marker');
 set(infmat(13,2),'callback','qelmts(-1,0,''Done'')','enable','on');
 set(infmat(13,3),'callback','qelmts(-1,0,''Cancel'')','enable','on');
 qelmtlistbox('Iterate');
 drawnow;

end
