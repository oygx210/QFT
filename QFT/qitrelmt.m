function qitrelmt(flag1,flag2)
% QITRELMT Set up user interface for iterating. (Utility Function)
%          QITRELMT sets up the user interface for iterating elements in
%          the controller matrix.  CNTDISP sets up the uicontrols that
%          allow the user to select which element is to be iterated.

% Author: Craig Borghesani
% Date: 9/5/93
% Revision: 2/18/96 12:17 AM V1.1 updates
% Copyright (c) 2003, Terasoft, Inc.


if all([flag1,flag2]==0),
 qclswin(0);
 f=gcf;

 bthan=get(f,'userdata');
 cont=get(bthan(3),'userdata');
 hint_bar = get(bthan(36),'userdata');

 set(bthan([1,8,29]),'enable','off');

 h=cntdisp(f,cont,6);

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
 set(f2,'pos',[pos(1:3),h+65]);

 if infmat(25,2)>1, proc_str=[' (',int2str(infmat(25,2)),')']; end

 if infmat(17,1)==0,
  infmat(17,1)=uicontrol(f2,'style','slider','pos',[100,h+5,180,15]);
  infmat(17,2)=uicontrol(f2,'style','slider','pos',[100,h+25,180,15]);
  infmat(17,3)=uicontrol(f2,'style','slider','pos',[100,h+45,180,15]);
  set(infmat(17,1:3),'enable','off');
  infmat(20,1)=uicontrol(f2,'style','text','pos',[0,h+5,95,17]);
  infmat(20,2)=uicontrol(f2,'style','text','pos',[0,h+25,95,17]);
  infmat(20,3)=uicontrol(f2,'style','text','pos',[0,h+45,95,17]);
  set(infmat(20,1:3),'horizontalalignment','right','back',fig_color);
  set(f2,'vis','on');
 else
  set(infmat(17,1),'style','slider','pos',[100,h+5,180,15]);
  set(infmat(17,2),'style','slider','pos',[100,h+25,180,15]);
  set(infmat(17,3),'style','slider','pos',[100,h+45,180,15]);
  set(infmat(20,1),'pos',[0,h+5,95,17]);
  set(infmat(20,2),'pos',[0,h+25,95,17]);
  set(infmat(20,3),'pos',[0,h+45,95,17]);
  set(infmat(17,1:3),'enable','off','vis','on','callback','');
  set(infmat(20,1:3),'string','','vis','on');
 end
 set(hint_bar,'string','Iterate elements of the current design');
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

 if all(flag1~=[8,9]),
  rad = gco;
  other_rads = get(rad,'userdata');
  set(other_rads,'value',0);
  set(infmat(17,1:3),'enable','off');
  set(infmat(20,1:3),'string','');
  set(infmat(21,1),'enable','on','userdata',rad);
  set(infmat(21,3),'enable','on','userdata',[flag2,rad]);
  set(hint_bar,'string','Use sliders to update selected element');
 end

end

if flag1==0 & flag2~=0,
 gain_val = abs(cont(flag2,1));
 set(infmat(17,1),'callback',['qitrelop(',int2str(flag1),',',int2str(flag2),')'],...
                  'value',gain_val,'enable','on');
 set(infmat(17,1),'min',10^((20*log10(gain_val)-5)/20),...
                  'max',10^((20*log10(gain_val)+5)/20));
 set(infmat(20,1),'string','gain=');

elseif any(flag1==[1 2]) & flag2~=0,
 str=['pole=';'zero='];
 set(infmat(17,1),'callback',['qitrelop(',int2str(flag1),',',int2str(flag2),')'],...
       'min',real(cont(flag2,1))*0.5,'max',real(cont(flag2,1))*1.5,...
       'value',real(cont(flag2,1)),'enable','on');
 set(infmat(20,1),'string',str(flag1,:));

elseif any(flag1==[3 4]) & flag2~=0,
 val=cont(flag2,1); val2=cont(flag2,2);
 set(infmat(17,2),'callback',['qitrelop(',int2str(flag1),',',int2str(flag2),')'],...
        'min',val/2,'max',val*1.5,'value',val,'enable','on');
 set(infmat(17,1),'callback',['qitrelop(',int2str(flag1),',',int2str(flag2),')'],...
        'min',val2/2,'max',val2*1.5,'value',val2,'enable','on');
 if flag1==3,
  set(infmat(20,2),'string','zeta(pole)=');
  set(infmat(20,1),'string','wn(pole)=');
 else
  set(infmat(20,2),'string','zeta(zero)=');
  set(infmat(20,1),'string','wn(zero)=');
 end

elseif flag1==5 & flag2~=0,
 set(infmat(17,2),'callback',['qitrelop(',int2str(flag1),',',int2str(flag2),')'],...
        'min',-87,'max',87,'value',cont(flag2,1),'enable','on');
 set(infmat(17,1),'callback',['qitrelop(',int2str(flag1),',',int2str(flag2),')'],...
        'min',cont(flag2,2)/2,'max',cont(flag2,2)*1.5,...
        'value',cont(flag2,2),'enable','on');
 set(infmat(20,2),'string','phase=');
 set(infmat(20,1),'string','w=');

elseif flag1==6 & flag2~=0,
 val=cont(flag2,1); val2=cont(flag2,2); val3=cont(flag2,3);
 set(infmat(17,3),'callback',['qitrelop(',int2str(flag1),',',int2str(flag2),')'],...
        'min',val/2,'max',val*1.5,'value',val,'enable','on');
 set(infmat(17,2),'callback',['qitrelop(',int2str(flag1),',',int2str(flag2),')'],...
        'min',val2/2,'max',val2*1.5,'value',val2,'enable','on');
 set(infmat(17,1),'callback',['qitrelop(',int2str(flag1),',',int2str(flag2),')'],...
        'min',val3/2,'max',val3*1.5,'value',val3,'enable','on');
 set(infmat(20,3),'string','zeta1(zero)=');
 set(infmat(20,2),'string','zeta2(pole)=');
 set(infmat(20,1),'string','wn=');

elseif flag1==8, % Cancel
 set(bthan([1,8,29]),'enable','on');
 v2=get(bthan(21),'userdata');
 vo2=get(bthan(22),'userdata');
 if infmat(9,1)==1,
  set([v2(:);vo2(:)],'vis','off');
 else
  set(v2,'vis','off');
 end
 set(bthan([19,20]),'userdata',[]);
 set(f2,'vis','off');

elseif flag1==9, % Done
 set(bthan([1,8,29]),'enable','on');
 cont_new = get(bthan(19),'userdata');
 lomat_new = get(bthan(20),'userdata');
 set(bthan(1),'userdata',lomat_new);
 set(bthan(3),'userdata',cont_new);
 set(bthan([19,20]),'userdata',[]);
 v=get(bthan(10),'userdata');
 v2=get(bthan(21),'userdata');
 if infmat(9,1)==1,
  vo2=get(bthan(22),'userdata');
  vo=get(bthan(17),'userdata');
  set([v,vo],'vis','off');
  set(bthan(22),'userdata',vo);
  set(bthan(17),'userdata',vo2);
 else
  set(v,'vis','off');
 end
 set(v2,'linestyle','-');
 set(v,'linestyle','--');
 set(bthan(10),'userdata',v2);
 set(bthan(21),'userdata',v);
 set(infmat(8,1),'enable','on');
 set(f2,'vis','off');

end
