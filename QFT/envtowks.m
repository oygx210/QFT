function envtowks(UIOperation)
% ENVTOWKS Environment to workspace. (Utility Function)
%          ENVTOWKS provides a gateway from the environment to the workspace
%          for the user during an IDE session.

% Author: Craig Borghesani
% 7/9/94
% Copyright (c) 2003, Terasoft, Inc.
%       $Revision: 1.6 $

if ~strcmp(UIOperation, 'Setup'),
   f2 = gcf;
   f = get(f2,'userdata');
   bthan=get(f,'userdata');
   infmat=get(bthan(16),'userdata');
   hint_bar = get(bthan(36),'userdata');
end


switch UIOperation,
   case('Setup'),
      f = gcf;
      bthan = get(f,'userdata');
      infmat = get(bthan(16),'userdata');
      del=25;

      hint_bar = get(bthan(36),'userdata');
      proc_str=[];
      if infmat(25,2)>1, proc_str=['(',int2str(infmat(25,2)),')']; end

      proc_num = int2str(infmat(25,2));
      f2 = findobj('tag',['qft6',proc_num]);

      if ~length(f2),
         ScreenSize     = get(0,'screensize');
         FigureWidth    = 225;
         FigureHeight   = 90;
         FigureLeft     = (ScreenSize(3) - FigureWidth)/2;
         FigureBottom   = (ScreenSize(4) - FigureHeight)/2;
         FigurePosition = [FigureLeft, FigureBottom, FigureWidth, FigureHeight];
         f2 = figure('name',['Export',proc_str],...
                     'numbertitle','off',...
                     'position',FigurePosition,...
                     'menubar','none',...
                     'vis','off',...
                     'userdata',f,...
                     'tag',['qft6',proc_num],...
                     'windowstyle','modal',...
                     'closerequestfcn','envtowks(''Cancel'')',...
                     'handlevisibility','callback');
         infmat(31,1) = f2;

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

         SubFrame1Height = 35;

         SubFrame1Position = [SubFramePosition(1), ...
                              MainFramePosition(4) - (SubFrame1Height + 5), ...
                              SubFramePosition(3), ...
                              SubFrame1Height];
         SubFrame1UILeft   = SubFrame1Position(1) + BorderWidth;
         SubFrame1UIBottom = sum(SubFrame1Position([2,4])) - 10;

         uicontrol(f2,'style','push',...
                      'pos',MainFramePosition,...
                      'enable','off');

         uicontrol('style','frame',...
                   'pos',SubFrame1Position,...
                   'fore','k');
         uicontrol(f2,'style','text',...
                      'string',' LTI Name',...
                      'pos',[SubFrame1UILeft, SubFrame1UIBottom, 150, 17],...
                      'horiz','left');

   % Numerator/Denominator
         SubFrame1UIBottom = SubFrame1UIBottom - 20;
         uicontrol(f2,'style','text',...
                      'string','Name',...
                      'pos',[SubFrame1UILeft, SubFrame1UIBottom, 30, 17],'horizontalalignment','left');
         edt(1) = uicontrol(f2,'style','edit',...
                               'pos',[SubFrame1UILeft+35, SubFrame1UIBottom, 160, 20],...
                               'backgroundcolor','w');

         set(edt,'horiz','left');
         uicontrol(f2,'style','push',...
                     'pos',[MainFrameLeft, 10, 60, 20],...
                     'string','OK',...
                     'callback','envtowks(''OK'')');
         uicontrol(f2,'style','push',...
                     'pos',[MainFrameLeft+65, 10, 60, 20],...
                     'string','Cancel',...
                     'callback','envtowks(''Cancel'')');

         set(bthan(35),'userdata',edt);
         set(bthan(16),'userdata',infmat);
         set(f2,'vis','on');

      else
         set(infmat(31,1),'vis','on');

      end

      set(hint_bar,'string','Enter variable name you wish to be passed to the workspace.');
      drawnow;

   case('OK'), % OK

      edt = get(bthan(35),'userdata');
      T = get(bthan(13),'userdata');

      cont=get(bthan(19),'userdata');
      lomat=get(bthan(20),'userdata');
      if isempty(cont),
         cont = get(bthan(3),'userdata');
         lomat=get(bthan(1),'userdata');
      end
%      if length(cont2),
%         cont_r(1,1) = cont_r(1,1)*cont2(1,1);
%         if T > 0,
%            cont_r(3,1) = cont_r(3,1)+cont2(3,1);
%            cont2(1:3,:) = [];
%         else
%            cont2(1:2,:) = [];
%         end
%         cont_r = [cont_r;cont2];
%      end

      edtObjString = get(edt(1),'string');
      if ~isempty(edtObjString),
         [var1,var2,var3]=cnt2zpk(cont,T);
         var4 = zpk(var1, var2, var3, T);
         assignin('base',edtObjString,var4);
      end

      set(hint_bar,'string','LTI object passed to workspace');

      set(infmat(31,1),'vis','off');

   case('Cancel'),

      set(infmat(31,1),'vis','off');

end
