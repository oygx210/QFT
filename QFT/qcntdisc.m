function qcntdisc(UIOperation)
%
% Utility function: QCNTDISC
%
% The purpose of this function is to view the effects of discretizing the
% continuous controller within SHAPE

% Author: Craig Borghesani
% Date: 4/15/95
% Copyright (c) 2003, Terasoft, Inc.

if strcmp(UIOperation,'Setup'),
   qclswin(0);
   f=gcf;
   bthan=get(f,'userdata');
   infmat=get(bthan(16),'userdata');
   T=get(bthan(13),'userdata');
   proc_str=[];
   if infmat(25,2)>1, proc_str=['(',int2str(infmat(25,2)),')']; end

   lomat=get(bthan(1),'userdata');
   w=lomat(1,:);
   proc_num = int2str(infmat(25,2));
   f2 = findobj('tag',['qft8',proc_num]);
   if ~length(f2),
      ScreenSize     = get(0,'screensize');
      FigureWidth    = 240;
      FigureHeight   = 195;
      FigureLeft     = (ScreenSize(3) - FigureWidth)/2;
      FigureBottom   = (ScreenSize(4) - FigureHeight)/2;
      FigurePosition = [FigureLeft, FigureBottom, FigureWidth, FigureHeight];
      f2 =figure('name',['Discretize ',proc_str],...
                         'numbertitle','off',...
                         'position',FigurePosition,...
                         'menubar','none',...
                         'resize','off',...
                         'vis','off',...
                         'userdata',f,...
                         'tag',['qft8',proc_num],...
                         'handlevisibility','callback');
      infmat(25,3) = f2;

% determine frame sizes depending upon figure size
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

      SubFrame1Height = 65;
      SubFrame2Height = 65;

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

      uicontrol('style','text',...
                'pos',[SubFrame1UILeft, SubFrame1UIBottom, 150, 17],...
                'string',' Conversion Information',...
                'horiz','left');

      SubFrame1UIBottom = SubFrame1UIBottom - 20;
      txt(1)=uicontrol('style','text',...
                       'pos',[SubFrame1UILeft, SubFrame1UIBottom, 130, 17],...
                       'string','Sample Frequency (Hz)',...
                       'horiz','left');
      edt(1)=uicontrol('style','edit',...
                       'pos',[SubFrame1UILeft+135, SubFrame1UIBottom, 75, 20],...
                       'horiz','right',...
                       'back','w');

      SubFrame1UIBottom = SubFrame1UIBottom - 25;
      txt(2)=uicontrol('style','text',...
                       'pos',[SubFrame1UILeft, SubFrame1UIBottom, 130, 17],...
                       'string','Sample Period (sec)',...
                       'horiz','left');
      edt(2)=uicontrol('style','edit',...
                       'pos',[SubFrame1UILeft+135, SubFrame1UIBottom, 75, 20],...
                       'horiz','right',...
                       'background','w');

      uicontrol('style','frame',...
                'pos',SubFrame2Position,...
                'fore','k');
      uicontrol('style','text',...
                'pos',[SubFrame2UILeft, SubFrame2UIBottom, 150, 17],...
                'string',' Conversion Method',...
                'horiz','left');

      SubFrame2UIBottom = SubFrame2UIBottom - 20;
      txt(4) = uicontrol('style','text',...
                   'pos',[SubFrame2UILeft, SubFrame2UIBottom, 135, 17],...
                   'string','Method',...
                   'horiz','left');

      pop = uicontrol('style','popup',...
                  'pos',[SubFrame2UILeft+135, SubFrame2UIBottom, 75, 20],...
                  'string',{'ZOH','FOH','IMP','Tustin','Matched','Prewarp'}, ...
                  'callback','qcntdisc(''Display'')',...
                  'background','w');

      SubFrame2UIBottom = SubFrame2UIBottom - 25;
      txt(3)=uicontrol('style','text',...
                       'pos',[SubFrame2UILeft, SubFrame2UIBottom, 130, 17],...
                       'string','Critical Frequency (Hz)',...
                       'horiz','left');
      edt(3)=uicontrol('style','edit',...
                       'pos',[SubFrame2UILeft+135, SubFrame2UIBottom, 75, 20],...
                       'background','w',...
                       'horiz','right',...
                       'string','1',...
                       'enable','off');

      set([edt,pop],'callback','qcntdisc(''Display'')');
      uicontrol('style','push',...
                'pos',[MainFrameLeft,10,60,20],...
                'string','Display',...
                'callback','qcntdisc(''Display'')');
      uicontrol('style','push',...
                'pos',[MainFrameLeft+65,10,60,20],...
                'string','Close',...
                'callback','qclswin(1)');
      set(bthan(16),'userdata',infmat);
      set(bthan(38),'userdata',[edt,pop]);
      drawnow;
      set(infmat(25,3),'vis','on');

   else
      set(infmat(25,3),'vis','on');

   end

elseif strcmp(UIOperation,'Display'),
   f2=gcf;
   f=get(f2,'userdata');
   bthan=get(f,'userdata');
   infmat=get(bthan(16),'userdata');
   han=get(bthan(38),'userdata');
   edt=han(1:3); pop=han(4);
   curobj = get(gcf,'currentobject');

   % make radio buttons mutually exclusive
   if ~strcmp(get(curobj,'style'),'popup'),
      Tf=str2num(get(edt(1),'string'));
      Tp=str2num(get(edt(2),'string'));
      if curobj==edt(1),
         if length(Tf),
            set(edt(2),'string',num2str(1/Tf));
         else
            set(edt(1),'string',num2str(1/Tp));
         end
      elseif curobj==edt(2),
         if length(Tp),
            set(edt(1),'string',num2str(1/Tp));
         else
            set(edt(2),'string',num2str(1/Tf));
         end
      end
   end

   % determine which method has been selected
   value=get(pop,'value');
   popstring = get(pop,'string');
   method = lower(popstring{value});

   if strcmp(method,'prewarp'),
      set(edt(3),'enable','on');
   else
      set(edt(3),'enable','off');
   end

   % obtain sampling period
   Ts=str2num(get(edt(2),'string'));

   if ~length(Ts), return; end

   drawnow;

   % obtain continuous controller and nominal loop
   cont=get(bthan(3),'userdata');
   lomat=get(bthan(1),'userdata');
   contcp=qcntbode(cont,lomat(1,:));
   nom_loop = lomat(2,:)./contcp;

   % get discrete transfer function of controller
   [z,p,k]=cnt2zpk(cont);
   if strcmp(method,'prewarp'),
      sysDT=c2d(zpk(z,p,k), Ts, method,str2num(get(edt(3),'string')));
   else
      sysDT=c2d(zpk(z,p,k), Ts, method);
   end

   locz=find(lomat(1,:)<=(1/Ts)*pi);
   wz = lomat(1,locz);
   contcpz=squeeze(freqresp(sysDT, wz)).';
   lomatz=[wz;nom_loop(locz).*contcpz];

   set(bthan(20),'userdata',lomatz);
   qnicplt(f);
   figure(infmat(25,3));

end
