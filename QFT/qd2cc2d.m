function sel = qd2cc2d(UIOperation, convDirection)
%
% Utility function: QD2CC2D
%
% The purpose of this function is to provide a dialog
% for converting a discrete controller to continuous
% for use within reduction.

% Author: Craig Borghesani
% Date: 4/14/03 2:32PM
% Copyright (c) 2003, Terasoft, Inc.

if strcmp(UIOperation,'Setup'),
   f=findobj(allchild(0),'tag','CAD Window');
   bthan=get(f,'userdata');
   infmat=get(bthan(16),'userdata');
   T=get(bthan(13),'userdata');

   ScreenSize     = get(0,'screensize');
   FigureWidth    = 355;
   FigureHeight   = 85;
   FigureLeft     = (ScreenSize(3) - FigureWidth)/2;
   FigureBottom   = (ScreenSize(4) - FigureHeight)/2;
   FigurePosition = [FigureLeft, FigureBottom, FigureWidth, FigureHeight];

   if strcmp(convDirection, 'd2c'),
      FigureTitle = 'Conversion Method: Discrete to Continuous';
   else
      FigureTitle = 'Conversion Method: Continuous to Discrete';
   end

   f2 =figure('name',FigureTitle,...
                        'numbertitle','off',...
                        'position',FigurePosition,...
                        'menubar','none',...
                        'resize','off',...
                        'vis','off',...
                        'userdata',f,...
                        'handlevisibility','callback',...
                        'windowstyle','modal');

% determine frame sizes depending upon figure size
   BorderWidth = 5;
   SectionBorderWidth = 3;
   MainFrameWidth = FigureWidth - BorderWidth*2;
   MainFramePosition = [FigureWidth - (MainFrameWidth + 5), ...
                        5, ...
                        MainFrameWidth, ...
                        FigureHeight-10];

   MainFrameLeft = MainFramePosition(1) + 5;
   MainFrameRight = FigureWidth - 2*BorderWidth;

   uicontrol(f2,'style','push',...
                  'pos',MainFramePosition,...
                  'enable','off');

   UIBottom = FigureHeight - 5;
   UILeft = 10;

   UIBottom = UIBottom - 20;
   uicontrol('style','text',...
                     'pos',[UILeft, UIBottom, 130, 17],...
                     'string','Select Conversion Method',...
                     'horiz','left');

   UIBottom = UIBottom - 25;
   uicontrol('style','text',...
                     'pos',[UILeft, UIBottom, 60, 17],...
                     'string','Method',...
                     'horiz','left');


   if strcmp(convDirection, 'd2c'),
      popupString = {'ZOH','Tustin','Matched','Prewarp'};
   else
      popupString = {'ZOH','FOH','IMP','Tustin','Matched','Prewarp'};
   end
   pop = uicontrol('style','popup',...
                     'pos',[UILeft+60, UIBottom, 60, 20],...
                     'string',popupString, ...
                     'callback','qd2cc2d(''Method Selection'')',...
                     'background','w');

   txt=uicontrol('style','text',...
                  'pos',[UILeft+125, UIBottom, 130, 17],...
                  'string','Critical Frequency (rad/sec)',...
                  'horiz','left',...
                  'enable','off');
   edt=uicontrol('style','edit',...
                  'pos',[UILeft+255, UIBottom, 60, 20],...
                  'background','w',...
                  'horiz','right',...
                  'string','1',...
                  'enable','off');
   set(pop,'userdata',[txt, edt]);

   uicontrol('style','push',...
               'pos',[MainFrameLeft,10,60,20],...
               'string','OK',...
               'callback','uiresume(gcf);');
   drawnow;
   set(f2,'vis','on');

   uiwait(f2);

   if ishandle(f2),

      pop = findobj(gcf,'style','popup');
      edt = findobj(gcf,'style','edit');
      popString = get(pop,'string');
      popValue = get(pop,'value');

      sel.Method = lower(popString{popValue});
      sel.Wc = str2num(get(edt,'string'));

   else
      sel.Method = 'zoh';
      sel.Wc = 1;

   end

   close(f2,'force');
   drawnow;

elseif strcmp(UIOperation, 'Method Selection'),

   curObject = gco;
   curString = get(curObject,'string');
   curValue = get(curObject,'value');
   if strcmp(curString{curValue}, 'Prewarp'),
      set(get(curObject,'userdata'),'enable','on');
   else
      set(get(curObject,'userdata'),'enable','off');
   end

elseif strcmp(UIOperation,'Closing'),

   close(gcf,'force');

end
