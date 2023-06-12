function qaxschng(UIOperation)
% QAXSCHNG Set up axis change interface. (Utility Function)
%          QAXSCHNG sets up the uicontrols for changing the axis vector
%          while still within the IDE environment.

% Author: Craig Borghesani
% 5/12/94
% Copyright (c) 2003, Terasoft, Inc.
%       $Revision: 1.5 $

if ~strcmp(UIOperation,'Setup'),
   f2 = gcf;
   f = get(f2,'userdata');

else
   f = gcf;

end
bthan = get(f,'userdata');
infmat = get(bthan(16),'userdata');
hint_bar = get(bthan(36),'userdata');

switch UIOperation,
   case('Setup'),

      proc_str=[];
      if infmat(25,2)>1, proc_str=['(',int2str(infmat(25,2)),')']; end

      if infmat(9,1)==1,
         phmin=infmat(1,1); phmins=num2str(phmin,3);
         phmax=infmat(1,2); phmaxs=num2str(phmax,3);
      elseif infmat(9,1)==2,
         phmins='N/A'; phmaxs='N/A'; phmin=0; phmax=5;
      elseif infmat(9,1)==3,
         phmin=infmat(2,3); phmins=num2str(phmin,3);
         phmax=infmat(2,4); phmaxs=num2str(phmax,3);
      end
      mgmin=infmat(1,3); mgmins=num2str(mgmin,3);
      mgmax=infmat(1,4); mgmaxs=num2str(mgmax,3);

      proc_num = int2str(infmat(25,2));
      f2 = findobj('tag',['qft5',proc_num]);
      if ~length(f2),
         ScreenSize     = get(0,'screensize');
         FigureWidth    = 300;
         FigureHeight   = 140;
         FigureLeft     = (ScreenSize(3) - FigureWidth)/2;
         FigureBottom   = (ScreenSize(4) - FigureHeight)/2;
         FigurePosition = [FigureLeft, FigureBottom, FigureWidth, FigureHeight];
         f2 = figure('name',['Axis Limits ',proc_str],...
                     'numbertitle','off',...
                     'position',FigurePosition,...
                     'menubar','none',...
                     'vis','off',...
                     'userdata',f,...
                     'tag',['qft5',proc_num],...
                     'windowstyle','modal',...
                     'closerequestfcn','qaxschng(''Cancel'')',...
                     'handlevisibility','callback');
         infmat(28,1) = f2;

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
         SubFrame2Height = 35;

         SubFrame1Position = [SubFramePosition(1), ...
                              MainFramePosition(4) - (SubFrame1Height + 5), ...
                              SubFramePosition(3), ...
                              SubFrame1Height];
         SubFrame1UILeft   = SubFrame1Position(1) + BorderWidth;
         SubFrame1UIBottom = sum(SubFrame1Position([2,4])) - 10;

         SubFrame2Position = [SubFramePosition(1), ...
                              SubFrame1Position(2) - (SubFrame2Height + 10), ...
                              SubFramePosition(3), ...
                              SubFrame2Height];
         SubFrame2UILeft   = SubFrame2Position(1) + BorderWidth;
         SubFrame2UIBottom = sum(SubFrame2Position([2,4])) - 10;

         uicontrol(f2,'style','push',...
                      'pos',MainFramePosition,...
                      'enable','off');

         uicontrol('style','frame',...
                           'pos',SubFrame1Position,...
                           'fore','k');
         t(1)=uicontrol(f2,'style','text',...
                           'string',' Phase Limits',...
                           'pos',[SubFrame1UILeft, SubFrame1UIBottom, 100, 17]);

         SubFrame1UIBottom = SubFrame1UIBottom - 20;
         t(2)=uicontrol(f2,'style','text',...
                           'string','Minimum',...
                           'pos',[SubFrame1UILeft, SubFrame1UIBottom, 60, 17]);
         infmat(29,1)=uicontrol(f2,'style','edit',...
                                   'pos',[SubFrame1UILeft+65, SubFrame1UIBottom, 60, 20],...
                                   'background','w',...
                                   'string',phmins,...
                                   'horiz','right');
         t(3)=uicontrol(f2,'style','text',...
                           'string','Maximum',...
                           'pos',[SubFrame1UILeft+135, SubFrame1UIBottom, 60, 17]);
         infmat(29,2)=uicontrol(f2,'style','edit',...
                                   'pos',[SubFrame1UILeft+200, SubFrame1UIBottom, 60, 20],...
                                   'background','w',...
                                   'string',phmaxs,...
                                   'horiz','right');

         uicontrol('style','frame',...
                           'pos',SubFrame2Position,...
                           'fore','k');
         t(4)=uicontrol(f2,'style','text',...
                           'string',' Magnitude Limits',...
                           'pos',[SubFrame2UILeft, SubFrame2UIBottom, 100, 17]);

         SubFrame2UIBottom = SubFrame2UIBottom - 20;
         t(5)=uicontrol(f2,'style','text',...
                           'string','Minimum',...
                           'pos',[SubFrame2UILeft, SubFrame2UIBottom, 60, 17]);
         infmat(29,3)=uicontrol(f2,'style','edit',...
                                   'pos',[SubFrame2UILeft+65, SubFrame2UIBottom, 60, 20],...
                                   'background','w',...
                                   'string',mgmins,...
                                   'horiz','right');
         t(6)=uicontrol(f2,'style','text',...
                           'string','Maximum',...
                           'pos',[SubFrame2UILeft+135, SubFrame2UIBottom, 60, 17]);
         infmat(29,4)=uicontrol(f2,'style','edit',...
                                   'pos',[SubFrame2UILeft+200, SubFrame2UIBottom, 60, 20],...
                                   'background','w',...
                                   'string',mgmaxs,...
                                   'horiz','right');

         h(1)=uicontrol(f2,'style','push',...
                           'pos',[MainFrameLeft,10,60,20],...
                           'string','Apply',...
                           'callback','qaxschng(''Apply'')');
         h(4)=uicontrol(f2,'style','push',...
                           'pos',[MainFrameLeft+65,10,60,20],...
                           'string','Done',...
                           'callback','qaxschng(''Done'')');
         h(3)=uicontrol(f2,'style','push',...
                           'pos',[MainFrameRight-125,10,60,20],...
                           'enable','off',...
                           'string','UnDo',...
                           'callback','qaxschng(''UnDo'')');
         h(2)=uicontrol(f2,'style','push',...
                           'pos',[MainFrameRight-60,10,60,20],...
                           'string','Cancel',...
                           'callback','qaxschng(''Cancel'')');
         set(bthan(16),'userdata',infmat);
         set(t([2,3,5,6]),'horizontalalignment','left');
         set(t([1,4]),'horizontalalignment','left');
         if infmat(9,1)==2,
            set([infmat(29,1:2),t(1:3)],'enable','off');
         end
         set(infmat(28,1),'vis','on');
         set(infmat(29,1),'userdata',h);

      else
         set(infmat(29,1),'string',phmins);
         set(infmat(29,2),'string',phmaxs);
         set(infmat(29,3),'string',mgmins);
         set(infmat(29,4),'string',mgmaxs);
         set(infmat(28,1),'vis','on');
         h=get(infmat(29,1),'userdata');
         set(h(3),'enable','off');
      end
      set(h(1),'userdata',[phmin,phmax]);
      set(h(2),'userdata',[mgmin,mgmax]);
      set(hint_bar,'string','Enter new values for desired axis limits');
      drawnow;

   case({'Apply', 'Done'}),
      if any(infmat(9,1)==[1,3]),
         phmin=str2num(get(infmat(29,1),'string'));
         phmax=str2num(get(infmat(29,2),'string'));
         mgmin=str2num(get(infmat(29,3),'string'));
         mgmax=str2num(get(infmat(29,4),'string'));
      else
         phmin=0; phmax=5;
         mgmin=str2num(get(infmat(29,3),'string'));
         mgmax=str2num(get(infmat(29,4),'string'));
      end
      if length([phmin,phmax,mgmin,mgmax])==4,
         if phmin<phmax & mgmin<mgmax,
            if infmat(9,1)==1,
               if any(abs(infmat(1,:)-[phmin,phmax,mgmin,mgmax])>0.5),
                  infmat(26,:)=infmat(1,:);
                  infmat(1,:)=[phmin,phmax,mgmin,mgmax];
                  infmat(27,:)=infmat(1,:);
                  set(bthan(16),'userdata',infmat);
                  qnicplt(f);

% redraw nichols chart
                  GridHandles = findall(infmat(24,1),'Tag','CSTgridLines');
                  if ~isempty(GridHandles),
                     delete(GridHandles);
                     axes(infmat(24,1));
                  end

                  if diff([phmin,phmax]) >= 360,
                     XTickValues = qaxesadjust([phmin,phmax]);
                     set(infmat(24,1),'xlim',[phmin,phmax],'ylim',[mgmin,mgmax],...
                              'xtick',XTickValues);

                  else
                     set(infmat(24,1),'xlim',[phmin,phmax],'ylim',[mgmin,mgmax],...
                              'xtickmode','auto');

                  end
                  copybnds(f);

                  if ~isempty(GridHandles),
                     ngrid;
                     set(infmat(24,1),'xlim',[phmin,phmax]);
                  end
                  drawnow;
               end
            elseif infmat(9,1)==2,
               if any(infmat(1,3:4)~=[mgmin,mgmax]),
                  infmat(26,:)=infmat(1,:);
                  infmat(1,3:4)=[mgmin,mgmax];
                  infmat(27,:)=infmat(1,:);
                  set(bthan(16),'userdata',infmat);
                  qmagplt(f);
                  set(infmat(24,1),'ylim',[mgmin,mgmax]);
               end
            elseif infmat(9,1)==3,
               if any([infmat(1,3:4),infmat(2,3:4)]~=[phmin,phmax,mgmin,mgmax]),
                  infmat(26,:)=[infmat(1,3:4),infmat(2,3:4)];
                  infmat(1,3:4)=[mgmin,mgmax];
                  infmat(2,3:4)=[phmin,phmax];
                  infmat(27,3:4)=infmat(2,3:4);
                  set(bthan(16),'userdata',infmat);
                  mgphplot(f);
                  set(infmat(24,1),'ylim',[mgmin,mgmax]);
                  set(infmat(24,2),'ylim',[phmin,phmax]);
               end
            end
            if strcmp(UIOperation,'Done'),
               set(infmat(28,1),'vis','off');

            else
               h = get(infmat(29,1),'userdata');
               set(h(2),'callback','qaxschng(''Cancel'')');
               set(h(3),'enable','on');
               figure(infmat(28,1));

            end
         else
            errordlg('MAX must be greater than MIN','Message','on');
         end
      else
         errordlg('Invalid input data','Message','on');
      end

   case({'UnDo', 'Cancel'}),

      h = get(infmat(29,1),'userdata');
      if infmat(9,1)==1,
         if strcmp(UIOperation,'Cancel'),
            infmat(26,:)=infmat(1,:);
            infmat(1,:)=[get(h(1),'userdata'),get(h(2),'userdata')];
         else
            tempaxs=infmat(1,:);
            infmat(1,:)=infmat(26,:);
            infmat(26,:)=tempaxs;
         end

% redraw nichols chart
         GridHandles = findall(infmat(24,1),'Tag','CSTgridLines');
         if ~isempty(GridHandles),
            delete(GridHandles);
            axes(infmat(24,1));
            ngrid;
            set(infmat(24,1),'xlim',infmat(1,1:2));
         end

         if diff(infmat(1,1:2)) >= 360,
            XTickValues = qaxesadjust(infmat(1,1:2));
            set(infmat(24,1),'xlim',infmat(1,1:2),'ylim',infmat(1,3:4),...
                     'xtick',XTickValues);

         else
            set(infmat(24,1),'xlim',infmat(1,1:2),'ylim',infmat(1,3:4),...
                     'xtickmode','auto');

         end
         set(bthan(16),'userdata',infmat);
         qnicplt(f);

      elseif infmat(9,1)==2,
         if strcmp(UIOperation,'Cancel'),
            infmat(26,:)=infmat(1,:);
            infmat(1,3:4)=get(h(2),'userdata');
         else
            tempaxs=infmat(1,:);
            infmat(1,:)=infmat(26,:);
            infmat(26,:)=tempaxs;
         end
         set(infmat(24,1),'ylim',infmat(1,3:4));
         set(bthan(16),'userdata',infmat);
         qmagplt(f);

      elseif infmat(9,1)==3,
         if strcmp(UIOperation,'Cancel'),
            infmat(26,:)=[infmat(1,3:4),infmat(2,3:4)];
            infmat(28,2:3) = infmat(1,1:2);
            axsm=get(h(2),'userdata');
            axsp=get(h(1),'userdata');

         else
            axsm=infmat(1,:);
            axsp=infmat(2,:);
            infmat(26,:)=[axsm(3:4),axsp(3:4)];
            infmat(28,2:3)=axsm(1,1:2);
            axsm=[infmat(28,2:3),infmat(26,1:2)];
            axsp=[infmat(28,2:3),infmat(26,3:4)];

         end
         set(infmat(24,1),'xlim',axsm(1:2),'ylim',axsm(3:4));
         set(infmat(24,1),'xlim',axsp(1:2),'ylim',axsp(3:4));
         infmat(1,:)=axsm; infmat(2,3:4)=axsp(3:4);
         set(bthan(16),'userdata',infmat);
         mgphplot(f);

      end

      if strcmp(UIOperation,'Cancel'),
         set(infmat(28,1),'vis','off');

      else
         if infmat(9,1)==1,
            set(infmat(29,1),'string',num2str(infmat(1,1),3));
            set(infmat(29,2),'string',num2str(infmat(1,2),3));
            set(infmat(29,3),'string',num2str(infmat(1,3),3));
            set(infmat(29,4),'string',num2str(infmat(1,4),3));

         elseif infmat(9,1)==2,
            set(infmat(29,3),'string',num2str(infmat(1,3),3));
            set(infmat(29,4),'string',num2str(infmat(1,4),3));

         elseif infmat(9,1)==3,
            set(infmat(29,1),'string',num2str(infmat(2,3),3));
            set(infmat(29,2),'string',num2str(infmat(2,4),3));
            set(infmat(29,3),'string',num2str(infmat(1,3),3));
            set(infmat(29,4),'string',num2str(infmat(1,4),3));

         end
         set(h(3),'enable','off');
         figure(infmat(28,1));
      end
end
