function qedtelmt(flag1,flag2)
% QEDTELMT Set up the edit interface. (Utility Function)
%          QEDTELMT sets up the user interface for editing elements in the
%          controller matrix.  CNTDISP sets up the uicontrols that allow
%          the user to select which element needs to be edited.

% Author: Craig Borghesani
% Date: 9/5/93
% Revised: 2/18/96 12:17 AM V1.1 updates
% Copyright (c) 2003, Terasoft, Inc.


if all([flag1,flag2]==0),
 qclswin(0);
 f=gcf;

 bthan=get(f,'userdata');
 cont=get(bthan(3),'userdata');
 hint_bar = get(bthan(36),'userdata');

 set(bthan([1,8,29]),'enable','off');

 h=cntdisp(f,cont,4);

 infmat=get(bthan(16),'userdata');
 f2=infmat(15,1);
 pos = get(f2,'pos');
 T=get(bthan(13),'userdata');
 lomat=get(bthan(1),'userdata');
 set(bthan(19),'userdata',cont);
 set(bthan(20),'userdata',lomat);
 set(bthan(27),'userdata',lomat);
 set(bthan(28),'userdata',cont);

 fig_color=[0.5,0.5,0.5];
 proc_str=[];
 set(f2,'pos',[pos(1:3),h+86]);

 if infmat(25,2)>1, proc_str=[' (',int2str(infmat(25,2)),')']; end

 if infmat(17,1)==0,
  infmat(17,1)=uicontrol(f2,'style','edit','pos',[152,h+5,100,20]);
  infmat(17,2)=uicontrol(f2,'style','edit','pos',[152,h+32,100,20]);
  infmat(17,3)=uicontrol(f2,'style','edit','pos',[152,h+59,100,20]);
  set(infmat(17,1:3),'background',[1,1,1],'enable','off','horiz','right');
  infmat(20,1)=uicontrol(f2,'style','text','pos',[45,h+5,105,17]);
  infmat(20,2)=uicontrol(f2,'style','text','pos',[45,h+32,105,17]);
  infmat(20,3)=uicontrol(f2,'style','text','pos',[45,h+59,105,17]);
  set(infmat(20,1:3),'horizontalalignment','right','back',fig_color);
  set(f2,'vis','on');
 else
  set(infmat(17,1),'style','edit','pos',[152,h+5,100,20]);
  set(infmat(17,2),'style','edit','pos',[152,h+32,100,20]);
  set(infmat(17,3),'style','edit','pos',[152,h+59,100,20]);
  set(infmat(20,1),'pos',[45,h+5,105,17]);
  set(infmat(20,2),'pos',[45,h+32,105,17]);
  set(infmat(20,3),'pos',[45,h+59,105,17]);
  set(infmat(17,1:3),'string','','enable','off','backgroundcolor',[1,1,1],...
                     'callback','','vis','on');
  set(infmat(20,1:3),'string','','vis','on');
 end
 set(hint_bar,'string','Edit elements of the current design');
 drawnow;
 set(bthan(16),'userdata',infmat);
 figure(f2);
end

if flag2~=0,
 f2=gcf;
 f=get(f2,'userdata');
 bthan=get(f,'userdata');
 infmat=get(bthan(16),'userdata');
 cont=get(bthan(19),'userdata');
 hint_bar = get(bthan(36),'userdata');

 if flag1~=8,
  rad = gco;
  other_rads = get(rad,'userdata');
  set(other_rads,'value',0);
  set(infmat(17,1:3),'enable','off','string','');
  set(infmat(20,1:3),'string','');
  set(infmat(21,[1,4]),'enable','on');
  set(infmat(21,[1,4]),'callback',['qelmts(',num2str(flag1),',',num2str(flag2),')']);
  set(infmat(21,1),'userdata',rad);
  set(hint_bar,'string','Enter values into the appropriate edit boxes');
 end

end

if flag1==0 & flag2~=0,
 set(infmat(17,1),'string',num2str(cont(flag2,1)),'enable','on');
 set(infmat(20,1),'string','gain=');

elseif any(flag1==[1 2]) & flag2~=0,
 str=['pole=';'zero='];
 set(infmat(17,1),'string',num2str(cont(flag2,1)),'enable','on');
 set(infmat(20,1),'string',str(flag1,:));

elseif any(flag1==[3 4]) & flag2~=0,
 set(infmat(17,2),'string',num2str(cont(flag2,1)),'enable','on');
 set(infmat(17,1),'string',num2str(cont(flag2,2)),'enable','on');
 if flag1==3,
  set(infmat(20,2),'string','zeta(pole)=');
  set(infmat(20,1),'string','wn(pole)=');
 else
  set(infmat(20,2),'string','zeta(zero)=');
  set(infmat(20,1),'string','wn(zero)=');
 end

elseif flag1==5 & flag2~=0,
 set(infmat(17,2),'string',num2str(cont(flag2,1)),'enable','on');
 set(infmat(17,1),'string',num2str(cont(flag2,2)),'enable','on');
 set(infmat(20,2),'string','phase=');
 set(infmat(20,1),'string','w=');

elseif flag1==6 & flag2~=0,
 set(infmat(17,3),'string',num2str(cont(flag2,1)),'enable','on');
 set(infmat(17,2),'string',num2str(cont(flag2,2)),'enable','on');
 set(infmat(17,1),'string',num2str(cont(flag2,3)),'enable','on');
 set(infmat(20,3),'string','zeta1(zero)=');
 set(infmat(20,2),'string','zeta2(pole)=');
 set(infmat(20,1),'string','wn=');

elseif flag1==8 & flag2~=0, % Cancel/Close
 set(bthan([1,8,29]),'enable','on');
 if length(get(bthan(30),'userdata')),
  v2=get(bthan(21),'userdata');
  vo2=get(bthan(22),'userdata');
  set(bthan(30),'userdata',[]);
  if infmat(9,1)==1,
   set([v2(:);vo2(:)],'vis','off');
  else
   set(v2,'vis','off');
  end
 end
 set(bthan([19,20]),'userdata',[]);
 set(f2,'vis','off');

end
