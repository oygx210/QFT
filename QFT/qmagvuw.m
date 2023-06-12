function qmagvuw(flag)
% QMAGVUW QView pulldown menu for PFSHAPE. (Utility Function)
%         QMAGVUW implements all the plot viewing commands: Full, In, Out,
%         Move, and Zoom.

% Author: Craig Borghesani
% 9/6/93
% Copyright (c) 2003, Terasoft, Inc.


str=str2mat('Continuous-time Filter Shaping','Discrete-time Filter Shaping');

f=gcf; a=gca;
bthan=get(f,'userdata');
infmat=get(bthan(16),'userdata');
lomat=get(bthan(1),'userdata');
lomat2=get(bthan(20),'userdata');
hint_bar = get(bthan(36),'userdata');
lo2=[];
if length(lomat2), lo2=lomat2(2:3,:); end
w=lomat(1,:); lo=lomat(2:3,:);
cont=get(bthan(3),'userdata');
ab=get(bthan(18),'userdata');
axs=infmat(1,:);
T=get(bthan(13),'userdata');
fig_color=[0.5,0.5,0.5];
proc_num = int2str(infmat(25,2));

if flag==1, % Full
 infmat(26,:)=axs;
 axs(1:2)=infmat(27,1:2);
 axs(3)=min([ab(:);20*log10(abs([lo(:);lo2(:)]))])-5;
 axs(4)=max([ab(:);20*log10(abs([lo(:);lo2(:)]))])+5;
elseif flag==2, % In
 infmat(26,:)=axs;
 amt=(axs(4)-axs(3))*0.2; axs=axs+[0 0 amt -amt];
elseif flag==3, % Out
 infmat(26,:)=axs;
 amt=(axs(4)-axs(3))*0.2; axs=axs+[0 0 -amt amt];
elseif flag==4, % Last
 tempaxs=axs;
 axs=infmat(26,:);
 infmat(26,:)=tempaxs;
elseif flag==11, % Move-initialization
 set(f,'pointer','crosshair');
 set(f,'windowbuttondownfcn','qmagvuw(10)');
elseif flag==10, % Move-getting first point
 pt=get(a,'currentpoint');
 set(f,'windowbuttondownfcn','qmagvuw(5)');
 set(a,'userdata',pt);
elseif flag==5, % Move-getting second point and setting new axis limits
 pt2=get(a,'currentpoint'); pt = get(a,'userdata');
 dely=pt2(1,2)-pt(1,2);
 infmat(26,:)=axs;
 axs=[axs(1:2),axs(3:4)+dely];
elseif flag==13, % Zoom-initialization
 set(f,'pointer','crosshair','windowbuttonmotionfcn','');
 set(f,'windowbuttondownfcn','qmagvuw(6)','windowbuttonupfcn','1;');
elseif flag==6, % Zoom
 pt0=get(a,'currentpoint');
 rbbox([get(f,'currentpoint'),0,0],get(f,'currentpoint'));
 drawnow;
 pt1=get(a,'currentpoint');
 minx=min([pt0(1,1),pt1(1,1)]); miny=min([pt0(1,2),pt1(1,2)]);
 maxx=max([pt0(1,1),pt1(1,1)]); maxy=max([pt0(1,2),pt1(1,2)]);
 minx=max([minx,axs(1)]); miny=max([miny,axs(3)]);
 maxx=min([maxx,axs(2)]); maxy=min([maxy,axs(4)]);
 infmat(26,:)=axs;
 axs=[minx,maxx,miny,maxy];
 set(f,'windowbuttonupfcn','');
elseif flag==9, % View/See
 cntdisp(f,cont,5);
elseif flag==15, % Refresh
 v1 = get(bthan(10),'userdata'); v2 = get(bthan(21),'userdata');
 init_vis = get(v2(1),'vis');
 set(v1,'vis','off');
 set(a,'vis','off');
 set(f,'color','k');
 drawnow;
 set(a,'vis','on');
 set(v1,'vis','on');
 set(v2,'vis',init_vis);
elseif flag == 17, % Zoom mode on/off
   zoomMenuItem = findobj(f,'tag','zoommode');
   if strcmp(get(zoomMenuItem,'checked'),'on'),
      set(zoomMenuItem,'checked','off');
   else
      set(zoomMenuItem,'checked','on');
   end

end
if any(flag==[1:7]),
 set(f,'windowbuttonmotionfcn','qmouse(''Floating'')','windowbuttonupfcn','',...
       'windowbuttondownfcn','');
 if diff(axs(1:2)) <= 0 | diff(axs(3:4)) <= 0,
  set(hint_bar,'string','Invalid region selected.  Please try again.');
 else
  set(a,'xlim',axs(1:2),'ylim',axs(3:4));
  infmat(1,:)=axs; set(bthan(16),'userdata',infmat);
  qmagplt(f);
  set(hint_bar,'string','Ready');
 end
 if length(findobj('tag',['qft3',proc_num])),
  if strcmp(get(infmat(15,1),'vis'),'on'),
   figure(infmat(15,1));
  end
 end
 if length(findobj('tag',['qft4',proc_num])),
  if strcmp(get(infmat(15,2),'vis'),'on'),
   figure(infmat(15,2));
  end
 end
 if length(findobj('tag',['qft1',proc_num])),
  if strcmp(get(infmat(14,1),'vis'),'on'),
   figure(infmat(14,1));
  end
 end
end
