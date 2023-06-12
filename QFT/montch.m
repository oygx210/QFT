function montch(qflag,qflag2)
% MONTCH Mouse-implemented notch element. (Utility Function)
%        MONTCH determines the notch value to be added from various
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
nom=get(bthan(2),'userdata');

w=lomat(1,:);
q=1; r=2;
if infmat(9,1)==2, q=[1;1]; r=2:3; end
lo=lomat(r,:);
T=get(bthan(13),'userdata');
delay=infmat(10,1);
am = infmat(24,1);
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
    set(f,'windowbuttonmotionfcn',['montch(',int2str(lcont+1),',1)'],...
          'windowbuttonupfcn','montch(0,2)');
    set(hint_bar,'string','Now, drag the marker to ADD the Notch element based on mag difference');
    infmat(3,3)=kw;
   end
  else
   kw=qfindfrq(pt(1,1),pt(1,2),lomat,a);
   if length(kw),
    set(infmat(6,1),'xdata',w(kw),'ydata',20*log10(abs(lo(1,kw))),'vis','on');
    set(f,'windowbuttonmotionfcn',['montch(',int2str(lcont+1),',1)'],...
          'windowbuttonupfcn','montch(0,2)');
    set(hint_bar,'string','Now, drag the marker to ADD the Notch element based on mag difference');
    infmat(3,3)=kw;
   end
  end
  set(bthan(16),'userdata',infmat);
 else
  set(f,'windowbuttonmotionfcn',['montch(',int2str(lcont),',1)'],...
        'windowbuttonupfcn','montch(0,2)');
  set(hint_bar,'string','Now, drag the marker to EDIT the Notch element');
 end
elseif qflag2==1,
 kw=infmat(3,3);
 wval=w(kw);
 p=infmat(4,3);
 am = infmat(24,1);
 obj=get(f,'currentobject');
 p_obj=get(obj,'parent');
 if strcmp(get(p_obj, 'type'), 'root'),
  a=get(f,'currentaxes');
 elseif p_obj==f,
  a=obj;
 else
  a=p_obj;
 end
 if a==am,
  pt=get(a,'currentpoint');
  axislims = infmat(1,:);
  if p~=0,
   delmag=abs(p)*10^(pt(1,2)/20)/abs(lomat(2,kw));
  else
   delmag=10^(pt(1,2)/20)/abs(lomat(2,kw));
  end
  if all(pt(1,1:2)>axislims([1,3]) & pt(1,1:2)<axislims([2,4])),
   ElementLocation = qflag;
   if lcont~=qflag,  % adding with the mouse
    if delmag>1, zeta(1)=0.5; zeta(2)=zeta(1)/delmag;
    else zeta(2)=0.5; zeta(1)=zeta(2)*delmag; end
    cont=[cont;0,0,0,0];
    cont(qflag,1:4)=[zeta(1),zeta(2),wval,6];
    if length(nom) & infmat(9,1)~=2,
     wlt=qfrqenh(wval,zeta(1+(zeta(1)>zeta(2))),w,T);
     [ncon,dcon]=cntextr(cont,T);
     if infmat(9,1)~=3,
      nlot=squeeze(freqresp(nom,wlt)).';
      clot=qcntbode(cont,wlt,T);
%      clot=qcntbode(cont,wlt,T).*exp(-i*wlt*delay);
      lot=nlot.*clot;
     else
      lot=qcntbode(cont,wlt,T).*exp(-i*wlt*delay);
     end
     lomat=[wlt;lot;ones(1,length(lot))];
     cpnv=ntchcplx(zeta(1),zeta(2),wval,lomat(1,:),T);
    else
     cpnv=ntchcplx(zeta(1),zeta(2),wval,w,T);
     lomat(r,:)=lomat(r,:).*cpnv(q,:);
    end
   else  % iterating with the mouse or continuation of adding
    zeta=cont(qflag,1:2);
    if delmag>1, zeta(1)=0.5; zeta(2)=zeta(1)/delmag;
    else zeta(2)=0.5; zeta(1)=zeta(2)*delmag; end
    cp=ntchcplx(cont(qflag,1),cont(qflag,2),wval,w,T);
    cont(qflag,1:4)=[zeta(1),zeta(2),wval,6];
    cpnv=ntchcplx(zeta(1),zeta(2),wval,w,T);
    lomat(r,:)=lomat(r,:).*(cpnv(q,:)./cp(q,:));
   end
   set(ElementsPopup,'value',cont(ElementLocation,4)+1);
   set(infmat(16,1),'string',num2str(cont(ElementLocation,1)));
   set(infmat(16,2),'string',num2str(cont(ElementLocation,2)));
   set(infmat(16,3),'string',num2str(cont(ElementLocation,3)));

   kw=find(lomat(1,:)>=w(kw)); kw=kw(1);
   if infmat(9,1)==1,
    set(infmat(6,1),'xdata',qfixfase(lomat(2,:),axs,kw),...
                    'ydata',20*log10(abs(lomat(2,kw))));
   elseif infmat(9,1)==2,
    set(infmat(6,1),'xdata',w(kw),'ydata',20*log10(abs(lomat(2,kw))));
   elseif infmat(9,1)==3,
    set(infmat(6,1),'xdata',lomat(1,kw),'ydata',20*log10(abs(lomat(2,kw))));
   end
   infmat(3,3)=kw;
   infmat(4,3)=cpnv(kw);
   set(bthan(20),'userdata',lomat);
   set(bthan(19),'userdata',cont);
   set(bthan(16),'userdata',infmat);
   if infmat(9,1)==1, qnicplt(f); setptr(f,'uddrag');
   elseif infmat(9,1)==2, qmagplt(f); setptr(f,'uddrag');
   elseif infmat(9,1)==3, mgphplot(f);
   end
   [ControllerString,ListboxInfo,ListboxValue] = cntstr(f,cont,ElementLocation);
   set(ElementsListbox, 'string', ControllerString,...
                        'userdata', ListboxInfo,...
                        'value',ListboxValue);
   set(hint_bar,'string','Notch element being implemented ...','fore','k');
%   drawnow;
  else
   set(hint_bar,'string','Must remain within axis limits','fore','r');
   setptr(f,'forbidden');
  end
 end

elseif qflag2==2,

 set(f,'windowbuttonmotionfcn','qmouse(''Floating'')','windowbuttondownfcn','',...
       'windowbuttonupfcn','');
 set(infmat(6,1),'userdata',['montch(',int2str(lcont),',0)']);
 set(hint_bar,'string','To continue with EDIT of the Notch element re-select the marker');
 set(infmat(13,2),'callback','qelmts(-1,0,''Done'')','enable','on');
 set(infmat(13,3),'callback','qelmts(-1,0,''Cancel'')','enable','on');
 qelmtlistbox('Iterate');
 drawnow;

end
