function qautoshp(UIOperation)
%
% Utility function: QAUTOSHP
%
% The purpose of this function is to view the effects of auto loop
% shaping the continuous controller within LPSHAPE.

% Author: Craig Borghesani
% Date: 1/26/97 8:14AM
% Copyright (c) 1994 by C. Borghesani, Y. Chait, O. Yaniv.

if strcmp(UIOperation,'Setup'),
   f=gcf;

else
   f2 = gcf;
   f = get(f2,'userdata');

end

bthan=get(f,'userdata');
infmat=get(bthan(16),'userdata');
cont=get(bthan(19),'userdata');
lomat=get(bthan(20),'userdata');
if isempty(cont),
   cont = get(bthan(3),'userdata');
   lomat=get(bthan(1),'userdata');
   set(bthan(19),'userdata',cont);
   set(bthan(20),'userdata',lomat);
end
hint_bar = get(bthan(36),'userdata');
QFTToolData = getappdata(f,'QFTToolData');
ElementsListbox = QFTToolData.Elements(17);

switch UIOperation
   case('Setup'),

      T=get(bthan(13),'userdata');
      proc_str=[];
      if infmat(25,2)>1, proc_str=['(',int2str(infmat(25,2)),')']; end
      proc_num = int2str(infmat(25,2));
      f2 = findobj('tag',['qft9',proc_num]);

      if ~length(f2),

         ScreenSize     = get(0,'screensize');
         FigureWidth    = 215;
         FigureHeight   = 160;
         FigureLeft     = (ScreenSize(3) - FigureWidth)/2;
         FigureBottom   = (ScreenSize(4) - FigureHeight)/2;
         FigurePosition = [FigureLeft, FigureBottom, FigureWidth, FigureHeight];
         f2 = figure('name',['Auto-Shape ',proc_str],...
                     'numbertitle','off',...
                     'position',FigurePosition,...
                     'menubar','none',...
                     'vis','off',...
                     'userdata',f,...
                     'tag',['qft9',proc_num],...
                     'windowstyle','modal',...
                     'closerequestfcn','qautoshp(''Cancel'')',...
                     'handlevisibility','callback');
         infmat(5,3) = f2;

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

         SubFrame1Height = 110;

         SubFrame1Position = [SubFramePosition(1), ...
                              MainFramePosition(4) - (SubFrame1Height + 5), ...
                              SubFramePosition(3), ...
                              SubFrame1Height];
         SubFrame1UILeft   = SubFrame1Position(1) + BorderWidth;
         SubFrame1UIBottom = sum(SubFrame1Position([2,4])) - 10;

         uicontrol(f2,'style','push',...
                        'pos',MainFramePosition,...
                        'enable','off');

         uicontrol(f2,'style','frame',...
                           'pos',SubFrame1Position,...
                           'fore','k');
         uicontrol(f2,'style','text',...
                           'string',' Auto-Shape Settings',...
                           'pos',[SubFrame1UILeft, SubFrame1UIBottom, 150, 17],...
                           'horiz','left');

         SubFrame1UIBottom = SubFrame1UIBottom - 20;
         txt(1)=uicontrol('style','text',...
                           'pos',[SubFrame1UILeft, SubFrame1UIBottom, 60, 17],...
                           'string','Zeta',...
                           'horiz','left');
         edt(1)=uicontrol('style','edit',...
                           'pos',[SubFrame1UILeft+65, SubFrame1UIBottom, 60, 20],...
                           'background','w',...
                           'string','0.6',...
                           'horiz','right');

         SubFrame1UIBottom = SubFrame1UIBottom - 25;
         txt(2)=uicontrol('style','text',...
                           'pos',[SubFrame1UILeft, SubFrame1UIBottom, 60, 17],...
                           'string','IC',...
                           'horiz','left');
         edt(2)=uicontrol('style','edit',...
                           'pos',[SubFrame1UILeft+65, SubFrame1UIBottom, 60, 20],...
                           'background','w',...
                           'string','5',...
                           'horiz','right');

         SubFrame1UIBottom = SubFrame1UIBottom - 25;
         txt(3)=uicontrol('style','text',...
                           'pos',[SubFrame1UILeft, SubFrame1UIBottom, 60, 17],...
                           'string','Iteration',...
                           'horiz','left');
         edt(3)=uicontrol('style','edit',...
                           'pos',[SubFrame1UILeft+65, SubFrame1UIBottom, 60, 20],...
                           'background','w',...
                           'string','4',...
                           'horiz','right');

         SubFrame1UIBottom = SubFrame1UIBottom - 25;
         rad(1)=uicontrol('style','radio',...
                           'pos',[SubFrame1UILeft, SubFrame1UIBottom, 80, 20],...
                           'string','Linear',...
                           'value',1,...
                           'callback','qautoshp(''RadioButton'')');
         rad(2)=uicontrol('style','radio',...
                           'pos',[SubFrame1UILeft+85, SubFrame1UIBottom, 80, 20],...
                           'string','Quadratic',...
                           'callback','qautoshp(''RadioButton'')');

         for c=1:2,
            set(rad(c),'userdata',rad(:,[1:(c-1),(c+1):2]));
         end

         btn(1) = uicontrol('style','push',...
                              'pos',[MainFrameLeft, 10, 60, 20],...
                              'string','Apply',...
                              'callback','qautoshp(''Apply'')');

         btn(2) = uicontrol('style','push',...
                              'pos',[MainFrameLeft+65, 10, 60, 20],...
                              'string','Done',...
                              'callback','qautoshp(''Done'')',...
                              'enable','off');

         btn(3) = uicontrol('style','push',...
                              'pos',[MainFrameLeft+131, 10, 60, 20],...
                              'string','Cancel',...
                              'callback','qautoshp(''Cancel'')');

         set(bthan(16),'userdata',infmat);
         set(bthan(37),'userdata',[edt, rad, btn]);
         drawnow;
         set(infmat(5,3),'vis','on');

      else
         set(infmat(5,3),'vis','on');

      end

   case('Apply'),

      han=get(bthan(37),'userdata');
      edt=han(1:3); rad=han(4:5);
      btn = han(6:8);
      curobj = get(gcf,'currentobject');

      nom = get(bthan(2),'userdata');
      wbs=get(bthan(11),'userdata');

      % obtain continuous controller and nominal loop
      contcp=qcntbode(cont,lomat(1,:));
      nom_loop = lomat(2,:)./contcp;

      phs  =get(bthan(33),'userdata');
      ubdb = get(bthan(32),'userdata');

      [index,li,ri,m,b0,ix] = convbd(ubdb,phs);

%      np0  = nom(1,:);
%      dp0  = nom(2,:);
      [np0, dp0] = tfdata(nom, 'v');

      % get settings
      zeta = str2num(get(edt(1),'string'));
      ic   = str2num(get(edt(2),'string'));
      iter = str2num(get(edt(3),'string'));

      linear_rad = get(rad(1),'value');
      quad_rad   = get(rad(2),'value');
      optype = (quad_rad==1) + (linear_rad==1)*2;

      [numc,denc]=autot3d(np0,dp0,wbs,ic,zeta,iter,optype,index,li,ri,m,b0,ix,ubdb);

      if length(numc),
         set(hint_bar,'string','Solution found!!');
         set(btn(2),'enable','on','userdata',{numc,denc});
         w = lomat(1,:);
         contcp2=freqcp(numc,denc,w);
         lomat2=[w;nom_loop.*contcp2;ones(1,length(w))];

         set(bthan(20),'userdata',lomat2);
         qnicplt(f);
         cont_new = cntpars(tf(numc,denc));
         ElementLocation = 1;
         [ControllerString, ListboxInfo, ListboxValue] = cntstr(f,cont_new,ElementLocation);
         set(ElementsListbox, 'string', ControllerString,...
                              'userdata', ListboxInfo,...
                              'value', ListboxValue);

      else
         set(hint_bar,'string','No solution found.');
         set(btn(2),'enable','off');
         set(bthan(19:20),'userdata',[]);

      end

   case('RadioButton'), % excluscivity of radio buttons

      h = gco;
      set(h,'value',1);
      set(get(h,'userdata'),'value',0);

   case('Done'), % Accepting new design

      han=get(bthan(37),'userdata');
      edt=han(1:3); rad=han(4:5);
      btn = han(6:8);
      data = get(btn(2),'userdata');
      cont_new = cntpars(data{1},data{2});
      ElementLocation = 1;
      [ControllerString, ListboxInfo, ListboxValue] = cntstr(f,cont_new,ElementLocation);
      set(ElementsListbox, 'string', ControllerString,...
                           'userdata', ListboxInfo,...
                           'value', ListboxValue);
      lomat = get(bthan(20),'userdata');
      set(bthan(1),'userdata',lomat);
      set(bthan(3),'userdata',cont_new);
      set(bthan([19,20]),'userdata',[]);
      set(btn(2),'enable','off');

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

      set(infmat(5,3),'vis','off');

   case('Cancel'),

      qclswin(1);

end
