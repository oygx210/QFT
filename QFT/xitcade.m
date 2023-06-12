function xitcade(f,flag)
% XITCADE Close all windows associated with shaping environments. (Utility)
%         XITCADE closes all windows associated with LPSHAPE and PFSHAPE.
%         It also offers the user a chance to save a design or cancel the close.

% Author: Craig Borghesani
% 9/6/93
% Copyright (c) 2003, Terasoft, Inc.

bthan=get(f,'userdata');
infmat=get(bthan(16),'userdata');
proc_num = int2str(infmat(25,2));
proc_str = [];
if infmat(25,2)>1, proc_str = [' (',proc_num,')']; end

if flag > 0,
   cont_r = get(bthan(3),'userdata');
   T = get(bthan(13),'userdata');
   cade = infmat(9,1)*2-(T > 0==0);
end

if flag==0,

   qclswin(0);
   f2 = findobj('tag',['qft7',proc_num]);
   if ~length(f2),

      ScreenSize     = get(0,'screensize');
      FigureWidth    = 211;
      FigureHeight   = 70;
      FigureLeft     = (ScreenSize(3) - FigureWidth)/2;
      FigureBottom   = (ScreenSize(4) - FigureHeight)/2;
      FigurePosition = [FigureLeft, FigureBottom, FigureWidth, FigureHeight];
      f2 = figure('name',['Exit IDE',proc_str],...
                  'numbertitle','off',...
                  'position',FigurePosition,...
                  'menubar','none',...
                  'vis','off',...
                  'resize','off',...
                  'userdata',f,...
                  'tag',['qft7',proc_num],...
                  'windowstyle','modal',...
                  'closerequestfcn',['xitcade(',int2str(f),',2)'],...
                  'handlevisibility','callback');

      infmat(32,1) = f2;

      BorderWidth = 5;
      SectionBorderWidth = 3;
      MainFrameWidth = FigureWidth - BorderWidth*2;
      MainFramePosition = [FigureWidth - (MainFrameWidth + 5), ...
                           5, ...
                           MainFrameWidth, ...
                           FigureHeight-10];
      SubFramePosition = [MainFramePosition(1) + BorderWidth, ...
                           MainFramePosition(2) + BorderWidth, ...
                           MainFramePosition(3) - 2*BorderWidth, ...
                           FigureHeight - 2*BorderWidth];

      MainFrameLeft = MainFramePosition(1) + 5;
      MainFrameRight = FigureWidth - 2*BorderWidth;

      uicontrol('style','push',...
                'pos',MainFramePosition,...
                'enable','off');

      uicontrol('style','text',...
                'string','Save Design?',...
                'pos',[20,35,100,20],...
                'horizontalalignment','left');

      uicontrol('style','push','string','Yes','pos',[MainFrameLeft,10,60,20],...
               'callback','cntsave(3)');
      uicontrol('style','push','string','No','pos',[MainFrameLeft+65,10,60,20],...
               'callback',['xitcade(',int2str(f),',1)']);
      uicontrol('style','push','string','Cancel','pos',[MainFrameLeft+131,10,60,20],...
               'callback',['xitcade(',int2str(f),',2)']);

      set(bthan(16),'userdata',infmat);
      drawnow;
   end
   set(infmat(32,1),'vis','on');

elseif flag==2, % Cancel
 set(infmat(32,1),'vis','off');

%elseif flag==3, % Save controller into tempfile and exit

% if cade==1, fname='shape.shp';
% elseif cade==2, fname='dshape.dsh';
% elseif cade==3, fname='fshape.fsh';
% elseif cade==4, fname='dfshape.dfs';
% elseif cade==5, fname='bodplot.bod';
% elseif cade==6, fname='dbodplot.dbo';
% end

% eval(['save ',fname,' cont_r T;']);

end

if any(flag==[1,3]),

   close(findobj(allchild(0),'tag',['qft1',proc_num]),'force');
   close(findobj(allchild(0),'tag',['qft2',proc_num]),'force');
   close(findobj(allchild(0),'tag',['qft3',proc_num]),'force');
   close(findobj(allchild(0),'tag',['qft4',proc_num]),'force');
   close(findobj(allchild(0),'tag',['qft5',proc_num]),'force');
   close(findobj(allchild(0),'tag',['qft6',proc_num]),'force');
   close(findobj(allchild(0),'tag',['qft7',proc_num]),'force');
   close(findobj(allchild(0),'tag',['qft8',proc_num]),'force');
   close(findobj(allchild(0),'tag',['qft9',proc_num]),'force');
   set(f,'vis','off');
   drawnow;
   close(f,'force');
   leftover = allchild(0);
   for k=1:length(leftover),
      if strcmp(get(leftover(k),'name'),'Message'),
         close(leftover(k),'force');
      end
   end

end
