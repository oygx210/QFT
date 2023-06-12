function qrelmts(UIOperation)
% QRELMTS Compute and reduce terms. (Utility Function)
%        QDELMTS computes the reduced frequency response of the elements
%        that are selected from within the Add or Edit window.

% Author: Craig Borghesani
% Date: 1/26/02 2:02PM
% Copyright (c) 2003, Terasoft, Inc.
%       $Revision: 1.5 $

if ~strcmp(UIOperation,'Setup'),
   f2=gcf;
   f=get(f2,'userdata');
else
   f=gcf;
end
bthan=get(f,'userdata');
infmat=get(bthan(16),'userdata');
hint_bar = get(bthan(36),'userdata');
QFTToolData = getappdata(f,'QFTToolData');

ElementEnvironment=infmat(9,1);
cont=get(bthan(19),'userdata');
lomat=get(bthan(20),'userdata');
if isempty(cont),
	cont = get(bthan(3),'userdata');
	lomat=get(bthan(1),'userdata');
	set(bthan(19),'userdata',cont);
	set(bthan(20),'userdata',lomat);
end
T=get(bthan(13),'userdata');
ElementsListbox = QFTToolData.Elements(17);

ListboxString = get(ElementsListbox,'string');
ListboxValue = get(ElementsListbox,'value');
ListboxInfo = get(ElementsListbox,'userdata');

% setting up for whether FSHAPE/DFSHAPE is being used
if any(ElementEnvironment==[1 3]), % SHAPE/DSHAPE/BODPLOT/DBODPLOT
 q=1; s=0;
elseif ElementEnvironment==2, % FSHAPE/DFSHAPE
 q=[1;1]; s=1;
end

proc_str = [];
if infmat(25,2)>1, proc_str=[' (',int2str(infmat(25,2)),')']; end

if strcmp(UIOperation, 'Setup') | strcmp(UIOperation, 'Use Weights'),
   red = [1,0,0];
   blue = [0,0,1];
   cyan = [0,1,1];

% store originally selected elements
   if strcmp(UIOperation,'Setup'),
      setappdata(f,'OriginalCont',cont);
      setappdata(f,'OriginalLomat',lomat);
      SelectedElements = ListboxInfo(ListboxValue);
      SelectedCont = cont(SelectedElements, :);

      if T == 0,
         SelectedCont=[1,NaN,NaN,0;
                       SelectedCont];
      else
         SelectedCont = [1,NaN,NaN,0;
                         0,0,NaN,0.5;
                         0,NaN,NaN,0.6;
                         SelectedCont];
      end

      SelectedCont=cntcvrt(SelectedCont,T);
      dcgain=cntdcgn(SelectedCont,T);
      SelectedCont(1,1)=dcgain;
      setappdata(f,'SelectedCont',SelectedCont);
      cont(SelectedElements, :) = [];
      LeftOverCont = cont;
      setappdata(f,'LeftOverCont',LeftOverCont);

   else
      SelectedCont = getappdata(f,'SelectedCont');
      LeftOverCont = getappdata(f,'LeftOverCont');

   end

   [z,p,k]=cnt2zpk(SelectedCont,T);

% convert z,p,k controller to continuous time (use modal dialog)
   if T > 0,
      sysCT = d2c(zpk(z, p, k, T), 'zoh');
      [z,p,k] = zpkdata(sysCT, 'v');
   end

   [do_it,repeat] = chkzp(z,p,T);
   wfunc=[];

   if strcmp(UIOperation,'Use Weights') & do_it,
      wfunc_han = get(bthan(26),'userdata');
      objects = get(wfunc_han(3),'userdata');
      obj_len = length(objects);
      if obj_len,
         ct = 1;
         y_scale = get(wfunc_han(11),'userdata');
         for k = 2:2:obj_len,
            if strcmp(get(objects(k),'vis'),'on'),
               ydata = get(objects(k),'ydata');
               ydata = ydata / max(ydata);
               w = get(objects(k),'xdata');
               a0 = ydata(1);
               a1 = (ydata(2) - ydata(1))/(w(2)-w(1));
               clr = get(objects(k),'userdata');
               wtype = all(clr==blue)+all(clr==cyan)*2+all(clr==red)*3;
               wfunc(ct,:) = [a0,a1,w(1),w(2),wtype];
               ct = ct + 1;
            end
         end
         wfunc = [wfunc;0.001,0,0,inf,3];
      end
   end

   if do_it,
      if repeat,
         [sysb,hsv]=qfwbal(z,p,k,wfunc,'z');
      else
         [r,p,k]=qzp2rp(z,p,k);
         [sysb,hsv]=qfwbal(r,p,k,wfunc,'r');
      end
      if ~sysb,
         errordlg('Weight function needs more points','Message');
         return;
      end
     proc_num = int2str(infmat(25,2));
     f2 = findobj('tag',['qft4',proc_num]);
     if ~length(f2),
         ScreenSize     = get(0,'screensize');
         FigureWidth    = 360;
         FigureHeight   = 365;
         FigureLeft     = (ScreenSize(3) - FigureWidth)/2;
         FigureBottom   = (ScreenSize(4) - (FigureHeight+20))/2;
         FigurePosition = [FigureLeft, FigureBottom, FigureWidth, FigureHeight];
         infmat(15,2) = figure('name',['Hankel Singular Values',proc_str],...
                               'numbertitle','off',...
                               'pos',FigurePosition,...
                               'vis','off',...
                               'menubar','none',...
                               'userdata',f,'tag',['qft4',proc_num],...
                               'closerequestfcn','qrelmts(''Cancel'')',...
                               'handlevisibility','callback');
         f2 = infmat(15,2);
         SubFrameHeight = 170;
         SubFrameBorder = 5;
         SubFrameWidth = FigureWidth - 2*SubFrameBorder;
         SubFrameLeft = SubFrameBorder;
         SubFrameBottom = 5;
         SubFrameUIBottom = SubFrameBottom + SubFrameHeight - 25;
         SubFrameUILeft = SubFrameLeft + 5;
         SubFrameUIRight = SubFrameLeft + SubFrameWidth - 5;
         AxesLeft = 35;
         AxesHeight = FigureHeight - SubFrameHeight - 45;
         AxesBottom = SubFrameBottom + SubFrameHeight + 30;
         AxesWidth = FigureWidth - 2*AxesLeft;

         xticks = 0:length(hsv)+1;
         h(1)=gca;
         set(h(1),'units','pixel',...
                  'pos',[AxesLeft, AxesBottom, AxesWidth, AxesHeight],...
                  'xtick',xticks,...
                  'yscale','log',...
                  'ygrid','on',...
                  'xgrid','on',...
                  'box','on',...
                  'nextplot','add',...
                  'xlim',[0,length(hsv)+1],...
                  'ylim',[10 .^[floor(log10(min(hsv))),ceil(log10(max(hsv)))]]);
         h(2)=line('ydata',hsv,...
                   'xdata',[1:length(hsv)],...
                   'color','m',...
                   'userdata',hsv);
         h(3)=line('ydata',hsv,...
                   'xdata',[1:length(hsv)],...
                   'color','r',...
                   'marker','x',...
                   'userdata',[AxesLeft, AxesBottom, AxesWidth, AxesHeight],...
                   'linestyle','none');
         h(10) = uicontrol(f2,'style','push',...
                      'enable','inactive',...
                      'pos',[SubFrameLeft, SubFrameBottom, SubFrameWidth, SubFrameHeight]);

         h(11) = uicontrol('style','text',...
                     'pos',[SubFrameUILeft, SubFrameUIBottom, 150, 17],...
                     'string','D2C Conversion',...
                     'horiz','left');

         SubFrameUIBottom = SubFrameUIBottom - 20;
         h(12) = uicontrol('style','text',...
                   'pos',[SubFrameUILeft, SubFrameUIBottom, 60, 17],...
                   'string','Method',...
                   'horiz','left');

         h(13) = uicontrol('style','popup',...
                     'pos',[SubFrameUILeft+60, SubFrameUIBottom, 70, 20],...
                     'string',{'ZOH','Tustin','Matched','Prewarp'}, ...
                     'callback','qrelmts(''Update D2C'')',...
                     'background','w');

         h(14)=uicontrol('style','text',...
                  'pos',[SubFrameUILeft+135, SubFrameUIBottom, 135, 17],...
                  'string','Critical Frequency (rad/sec)',...
                  'horiz','left',...
                  'enable','off');
         h(15)=uicontrol('style','edit',...
                  'pos',[SubFrameUILeft+270, SubFrameUIBottom, 60, 20],...
                  'background','w',...
                  'horiz','right',...
                  'string','1',...
                  'callback','qrelmts(''Update D2C'')',...
                  'enable','off');

         SubFrameUIBottom = SubFrameUIBottom - 25;
         h(16) = uicontrol('style','text',...
                     'pos',[SubFrameUILeft, SubFrameUIBottom, 150, 17],...
                     'string','C2D Conversion',...
                     'horiz','left');

         SubFrameUIBottom = SubFrameUIBottom - 20;
         h(17) = uicontrol('style','text',...
                   'pos',[SubFrameUILeft, SubFrameUIBottom, 60, 17],...
                   'string','Method',...
                   'horiz','left');

         h(18) = uicontrol('style','popup',...
                     'pos',[SubFrameUILeft+60, SubFrameUIBottom, 70, 20],...
                     'string',{'ZOH','FOH','IMP','Tustin','Matched','Prewarp'}, ...
                     'callback','qrelmts(''Update C2D'')',...
                     'background','w');

         h(19)=uicontrol('style','text',...
                  'pos',[SubFrameUILeft+135, SubFrameUIBottom, 135, 17],...
                  'string','Critical Frequency (rad/sec)',...
                  'horiz','left',...
                  'enable','off');
         h(20)=uicontrol('style','edit',...
                  'pos',[SubFrameUILeft+270, SubFrameUIBottom, 60, 20],...
                  'background','w',...
                  'horiz','right',...
                  'string','1',...
                  'enable','off');

         SubFrameUIBottom = SubFrameUIBottom - 25;
         h(20) = uicontrol('style','text',...
                     'pos',[SubFrameUILeft, SubFrameUIBottom, 150, 17],...
                     'string','Order Selection',...
                     'horiz','left');

         SubFrameUIBottom = SubFrameUIBottom - 20;
         h(4)=uicontrol(f2,'style','text',...
                           'pos',[SubFrameUILeft, SubFrameUIBottom, 60, 17],...
                           'string','Order',...
                           'horiz','left');
         pstring = sprintf('%d|', 1:length(hsv));
         pstring(end) = [];
         h(5)=uicontrol(f2,'style','popup',...
                           'pos',[SubFrameUILeft+60, SubFrameUIBottom, 70, 20],...
                           'string',pstring,...
                           'background','w','horiz','right');

         SubFrameUIBottom = SubFrameUIBottom - 30;
         h(6)=uicontrol(f2,'style','push',...
                           'pos',[SubFrameUILeft, SubFrameUIBottom, 60, 20],...
                           'string','Reduce',...
                           'callback','qrelmts(''Reduce'')');
         h(7)=uicontrol(f2,'style','push','pos',[SubFrameUILeft+65, SubFrameUIBottom, 100, 20],...
                           'string','Define Weights',...
                           'callback','qwatecad(0)');
         h(8)=uicontrol(f2,'style','push',...
                           'pos',[SubFrameUIRight-125, SubFrameUIBottom, 60, 20],...
                           'string','Done',...
                           'callback','qrelmts(''Done'')');
         h(9)=uicontrol(f2,'style','push',...
                           'pos',[SubFrameUIRight-60, SubFrameUIBottom, 60, 20],...
                           'string','Cancel',...
                           'callback','qrelmts(''Cancel'')');
         set(bthan(25),'userdata',h);

         h2 = uimenu('label','Line','enable','off');
         h2_sub(1)=uimenu(h2,'label','Add Line',...
                             'callback',...
                             'set(gcf,''pointer'',''crosshair'',''windowbuttonupfcn'',''qaddobj(0)'',''windowbuttondownfcn'','''',''windowbuttonmotionfcn'','''')');
         h2_sub(2)=uimenu(h2,'label','Add Point',...
                             'callback',...
                             'set(gcf,''pointer'',''crosshair'',''windowbuttonupfcn'',''qedtobj(0)'',''windowbuttondownfcn'','''',''windowbuttonmotionfcn'','''')');
         h2_sub(3)=uimenu(h2,'label','Move',...
                             'callback',...
                             'qbtnkill;set(gcf,''pointer'',''fleur'',''windowbuttondownfcn'',''set(gca,''''userdata'''',get(gca,''''currentpoint''''));set(gcf,''''windowbuttonmotionfcn'''',''''qmoveobj'''')'')');
         h2_sub(4)=uimenu(h2,'label','Delete',...
                              'callback','qdelobj(0)');
         h2_sub(5)=uimenu(h2,'label','Break',...
                              'callback','qdelobj(1)');
         h2_sub(6)=uimenu(h2,'label','Connect',...
                             'callback','qdelobj(2)');
         h2_sub(7)=uimenu(h2,'label','Type');
         uimenu(h2_sub(7),'label','Input',...
                          'callback','qdelobj(6)');
         uimenu(h2_sub(7),'label','Output',...
                          'callback','qdelobj(7)');
         uimenu(h2_sub(7),'label','Both',...
                          'callback','qdelobj(8)');
         h2_sub(8)=uimenu(h2,'label','Select');
         uimenu(h2_sub(8),'label','All',...
                          'callback','qdelobj(5)');
         uimenu(h2_sub(8),'label','None',...
                          'callback','qdelobj(4)');
         set(h2_sub(2:8),'enable','off');

         o=uimenu('label','Options','enable','off');
         o_sub(1)=uimenu(o,'label','Full',...
                           'callback','qzoomaxs(3)');
         o_sub(2)=uimenu(o,'label','Zoom',...
                           'callback','qzoomaxs(0)');
         o_sub(3)=uimenu(o,'label','Clear',...
                           'callback','qdelobj(3)');
         o_sub(4)=uimenu(o,'label','Open...',...
                           'callback','qputobj(0)',...
                           'separator','on');
         o_sub(5)=uimenu(o,'label','Save...',...
                           'callback','qputobj(1)');
         han = [infmat(15,2),h(1),h2,h2_sub,o,o_sub];
         set(bthan(26),'userdata',han);
         set(infmat(15,2),'vis','on','resizefcn','qrelmts(''Resize'')');
         set(bthan(16),'userdata',infmat);
         if T == 0,
            set(h([11:13,16:18]),'enable','off');
         end

      else
         set(f2,'pointer','arrow',...
                'windowbuttondownfcn','',...
                'windowbuttonupfcn','',...
                'windowbuttonmotionfcn','',...
                'name',['Hankel Singular Values',proc_str]);
         h = get(bthan(25),'userdata');
         h2 = get(bthan(26),'userdata');
         set(h2([3,12]),'enable','off');
         set(get(h2(3),'userdata'),'color',[0,0,0]);
         pos=get(h(3),'userdata');
         xticks = 0:length(hsv)+1;
         set(h(1),'pos',pos,...
                  'xtick',xticks,...
                  'xlim',[0,length(hsv)+1],...
                  'ylim',[10 .^[floor(log10(min(hsv))),ceil(log10(max(hsv)))]],...
                  'yscale','log',...
                  'xscale','linear');
         set(h(2),'ydata',hsv,...
                  'xdata',1:length(hsv));
         set(h(3),'ydata',hsv,...
                  'xdata',1:length(hsv),...
                  'vis','on');

         pstring = sprintf('%d|', 1:length(hsv));
         pstring(end) = [];
         set(h(5),'value',1,'string',pstring);
         set(h([4,5,6,7]),'enable','on');
         set(h(6),'callback','qrelmts(''Reduce'')');
         set(h(8),'callback','qrelmts(''Done'')');
         set(h(9),'callback','qrelmts(''Cancel'')');
         figure(infmat(15,2));
      end

      set(bthan(29),'userdata',sysb);
      set(hint_bar,'string','Plotting Hankel Singular Values');
   end

elseif strcmp(UIOperation,'Resize'),

   h = get(bthan(25),'userdata');

   f2 = infmat(15,2);
   FigurePosition = get(f2,'pos');
   FigureLeft = FigurePosition(1);
   FigureBottom = FigurePosition(2);
   FigureWidth = FigurePosition(3);
   FigureHeight = FigurePosition(4);

   SubFrameHeight = 170;
   SubFrameBorder = 5;
   SubFrameWidth = FigureWidth - 2*SubFrameBorder;
   SubFrameLeft = SubFrameBorder;
   SubFrameBottom = 5;
   SubFrameUIBottom = SubFrameBottom + SubFrameHeight - 25;
   SubFrameUILeft = SubFrameLeft + 5;
   SubFrameUIRight = SubFrameLeft + SubFrameWidth - 5;
   AxesLeft = 35;
   AxesHeight = FigureHeight - SubFrameHeight - 45;
   AxesBottom = SubFrameBottom + SubFrameHeight + 30;
   AxesWidth = FigureWidth - 2*AxesLeft;

   set(h(1),'pos',[AxesLeft, AxesBottom, AxesWidth, AxesHeight]);
   set(h(3),'userdata',[AxesLeft, AxesBottom, AxesWidth, AxesHeight]);
   set(h(10),'pos',[SubFrameLeft, SubFrameBottom, SubFrameWidth, SubFrameHeight]);

   set(h(11),'pos',[SubFrameUILeft, SubFrameUIBottom, 150, 17]);
   SubFrameUIBottom = SubFrameUIBottom - 20;
   set(h(12),'pos',[SubFrameUILeft, SubFrameUIBottom, 60, 17]);
   set(h(13),'pos',[SubFrameUILeft+60, SubFrameUIBottom, 70, 20]);
   set(h(14),'pos',[SubFrameUILeft+135, SubFrameUIBottom, 135, 17]);
   set(h(15),'pos',[SubFrameUILeft+270, SubFrameUIBottom, 60, 20]);
   SubFrameUIBottom = SubFrameUIBottom - 25;
   set(h(16),'pos',[SubFrameUILeft, SubFrameUIBottom, 150, 17]);
   SubFrameUIBottom = SubFrameUIBottom - 20;
   set(h(17),'pos',[SubFrameUILeft, SubFrameUIBottom, 60, 17]);
   set(h(18),'pos',[SubFrameUILeft+60, SubFrameUIBottom, 70, 20]);
   set(h(19),'pos',[SubFrameUILeft+135, SubFrameUIBottom, 135, 17]);
   set(h(20),'pos',[SubFrameUILeft+270, SubFrameUIBottom, 60, 20]);
   SubFrameUIBottom = SubFrameUIBottom - 25;
   set(h(20),'pos',[SubFrameUILeft, SubFrameUIBottom, 150, 17]);
   SubFrameUIBottom = SubFrameUIBottom - 20;
   set(h(4),'pos',[SubFrameUILeft, SubFrameUIBottom, 60, 17]);
   set(h(5),'pos',[SubFrameUILeft+60, SubFrameUIBottom, 70, 20]);
   SubFrameUIBottom = SubFrameUIBottom - 30;
   set(h(6),'pos',[SubFrameUILeft, SubFrameUIBottom, 60, 20]);
   set(h(7),'pos',[SubFrameUILeft+65, SubFrameUIBottom, 100, 20]);
   set(h(8),'pos',[SubFrameUIRight-125, SubFrameUIBottom, 60, 20]);
   set(h(9),'pos',[SubFrameUIRight-60, SubFrameUIBottom, 60, 20]);

elseif strcmp(UIOperation,'Reduce'),
   T=get(bthan(13),'userdata');
   h = get(bthan(25),'userdata');
   lomat=getappdata(f,'OriginalLomat');
   OriginalCont=getappdata(f,'OriginalCont');
   delay=infmat(10,1);
   q=1; s=0;
   if infmat(9,1)==2, q=[1;1]; s=1; end
   loc=qcntbode(OriginalCont,lomat(1,:),T);
   sysb=get(bthan(29),'userdata');
   [sr,sc] = size(sysb);
   LeftOverCont = getappdata(f,'LeftOverCont');
   val=get(h(5),'value');
   if val>0 & val<=(sc-2) & (~rem(val,1)),
      set(hint_bar,'string','Performing model order reduction');
      ar = sysb(1:val,1:val);
      br = sysb(1:val,sc-1);
      cr = sysb(sr-1,1:val);
      dr = sysb(sr-1,sc-1);
      [z,p]=ss2zp(ar,br,cr,dr);
      sysDT = zpk(z, p, 1, T);

% obtain discrete-time representation of the continuous z, p
      if T > 0,
         h18string = get(h(18),'string');
         h18value = get(h(18),'value');
         if strcmp(lower(h18string{h18value}), 'prewarp'),
            sysDT = c2d(zpk(z, p, 1), T, lower(h18string{h18value}), str2num(get(h(20),'string')));
         else
            sysDT = c2d(zpk(z, p, 1), T, lower(h18string{h18value}));
         end
      end

      ReducedCont=cntpars(sysDT, LeftOverCont);
      locmd=qcntbode(ReducedCont,lomat(1,:),T);
      lomat(2:2+s,:)=lomat(2:2+s,:).*(locmd(q,:)./loc(q,:));
      ElementLocation = size(ReducedCont,1);
      set(bthan(20),'userdata',lomat);
      set(bthan(19),'userdata',ReducedCont);
      setappdata(f,'ReducedCont',ReducedCont);
      if infmat(9,1)==1, qnicplt(f);
      elseif infmat(9,1)==2, qmagplt(f);
      elseif infmat(9,1)==3, mgphplot(f);
      end
      [ControllerString,ListboxInfo,ListboxValue] = cntstr(f,ReducedCont,ElementLocation);
      set(ElementsListbox, 'string', ControllerString,...
                           'userdata', ListboxInfo,...
                           'value',ListboxValue);

      qelmtlistbox('Edit');

   else
      errordlg('Improper input.  Check values.','Message','on');

   end

elseif strcmp(UIOperation,'Update D2C'),

   h = get(bthan(25),'userdata');
   h2 = get(bthan(26),'userdata');

   SelectedCont = getappdata(f,'SelectedCont');

   [z,p,k]=cnt2zpk(SelectedCont,T);

% convert z,p,k controller to continuous time (use modal dialog)
   h13string = get(h(13),'string');
   h13value = get(h(13),'value');
   if strcmp(lower(h13string{h13value}), 'prewarp'),
      set(h(14:15),'enable','on');
      sysCT = d2c(zpk(z, p, k, T), lower(h13string{h13value}), str2num(get(h(15),'string')));
   else
      set(h(14:15),'enable','off');
      sysCT = d2c(zpk(z, p, k, T), lower(h13string{h13value}));
   end
   [z,p,k] = zpkdata(sysCT, 'v');

   [do_it,repeat] = chkzp(z,p,T);
   wfunc=[];

   if repeat,
      [sysb,hsv]=qfwbal(z,p,k,wfunc,'z');
   else
      [r,p,k]=qzp2rp(z,p,k);
      [sysb,hsv]=qfwbal(r,p,k,wfunc,'r');
   end

   set(h2([3,12]),'enable','off');
   set(get(h2(3),'userdata'),'color',[0,0,0]);
   xticks = 0:length(hsv)+1;
   set(h(1),'xtick',xticks,...
            'xlim',[0,length(hsv)+1],...
            'ylim',[10 .^[floor(log10(min(hsv))),ceil(log10(max(hsv)))]]);
   set(h(2),'ydata',hsv,...
            'xdata',1:length(hsv));
   set(h(3),'ydata',hsv,...
            'xdata',1:length(hsv));

   set(bthan(29),'userdata',sysb);
   set(hint_bar,'string','Plotting Hankel Singular Values');

elseif strcmp(UIOperation,'Update C2D'),

   h = get(bthan(25),'userdata');

   h18string = get(h(18),'string');
   h18value = get(h(18),'value');
   if strcmp(lower(h18string{h18value}), 'prewarp'),
      set(h(19:20),'enable','on');
   else
      set(h(19:20),'enable','off');
   end

elseif strcmp(UIOperation,'Done'), % Done
   ReducedCont = getappdata(f,'ReducedCont');
   han2 = get(bthan(26),'userdata');
   if ~isempty(ReducedCont),
      set(bthan(3),'userdata',get(bthan(19),'userdata'));
      set(bthan(1),'userdata',get(bthan(20),'userdata'));
      w_func = get(han2(3),'userdata');
      delete(w_func);
      set(han2(3),'userdata',[]);
      set(bthan([1,8,29]),'enable','on');
      set(infmat(15,2),'vis','off');
   else
      set(han2(3),'userdata',[]);
      set(bthan([1,8,29]),'enable','on');
      set(infmat(15,2),'vis','off');
      return;
   end

   v=get(bthan(10),'userdata');
   v2=get(bthan(21),'userdata');
   if infmat(9,1)==1,
      vo2=get(bthan(22),'userdata');
      vo=get(bthan(17),'userdata');
      set([v,vo],'xdata',0,'ydata',0);
      set(bthan(22),'userdata',vo);
      set(bthan(17),'userdata',vo2);
      set(v2,'linestyle','-');
      set(v,'linestyle',':');
      set(bthan(10),'userdata',v2);
      set(bthan(21),'userdata',v);
   else
      if infmat(9,1)==1, qnicplt(f);
      elseif infmat(9,1)==2, qmagplt(f);
      elseif infmat(9,1)==3, mgphplot(f);
      end
   end

   ElementLocation = size(ReducedCont,1);
   [ControllerString,ListboxInfo,ListboxValue] = cntstr(f,ReducedCont,ElementLocation);
   set(ElementsListbox, 'string', ControllerString,...
                        'userdata', ListboxInfo,...
                        'value',ListboxValue);

   qelmtlistbox('Edit');

elseif strcmp(UIOperation, 'Cancel'), % Cancel
   OriginalCont=getappdata(f,'OriginalCont');
   OriginalLomat =getappdata(f,'OriginalLomat');

   set(bthan(19),'userdata',[]);
   set(bthan(20),'userdata',[]);
   set(bthan(1),'userdata',OriginalLomat);
   set(bthan(3),'userdata',OriginalCont);

   v2=get(bthan(21),'userdata');
   vo2=get(bthan(22),'userdata');
   if infmat(9,1)==1,
      set([v2(:);vo2(:)],'vis','off');
   else
      set(v2,'xdata',0,'ydata',0);
   end

   han2 = get(bthan(26),'userdata');
   w_func = get(han2(3),'userdata');
   delete(w_func);
   set(han2(3),'userdata',[]);
   set(bthan([1,8,29]),'enable','on');
   set(infmat(15,2),'vis','off');

   ElementLocation = size(OriginalCont,1);
   [ControllerString,ListboxInfo,ListboxValue] = cntstr(f,OriginalCont,ElementLocation);
   set(ElementsListbox, 'string', ControllerString,...
                        'userdata', ListboxInfo,...
                        'value',ListboxValue);

   qelmtlistbox('Edit');

end
