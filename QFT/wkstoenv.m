function wkstoenv(UIOperation)
%
% Utility function: WKSTOENV
%
% The purpose of this function is to import
% LTI objects into the shaping environment.

% Author: Craig Borghesani
% Date: 4/18/03 2:46PM
% Copyright (c) 2003, Terasoft, Inc.

if strcmp(UIOperation,'Setup'),
   f = gcf;
   ScreenSize     = get(0,'screensize');
   FigureWidth    = 240;
   FigureHeight   = 200;
   FigureLeft     = (ScreenSize(3) - FigureWidth)/2;
   FigureBottom   = (ScreenSize(4) - FigureHeight)/2;
   FigurePosition = [FigureLeft, FigureBottom, FigureWidth, FigureHeight];
   f2 =figure('name','Import',...
                        'numbertitle','off',...
                        'position',FigurePosition,...
                        'menubar','none',...
                        'resize','off',...
                        'vis','off',...
                        'userdata',f,...
                        'windowstyle','modal');

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

   SubFrame1Height = FigureHeight - 50;

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

   uicontrol('style','text',...
               'pos',[SubFrame1UILeft, SubFrame1UIBottom, 150, 17],...
               'string',' Select LTI Object',...
               'horiz','left');

   lb = uicontrol('style','listbox',...
               'pos',[SubFrame1UILeft, SubFrame1Position(2)+5, SubFrame1Position(3)-10, SubFrame1Position(4)-15],...
               'back','w',...
               'callback','wkstoenv(''Select'')',...
               'fontname','courier');

   wksData = evalin('base','whos');
   wksString = {};
   for k = 1:length(wksData),
      if (strcmp(wksData(k).class, 'zpk') | ...
         strcmp(wksData(k).class, 'tf') | ...
         strcmp(wksData(k).class, 'ss') | ...
         strcmp(wksData(k).class, 'frd')) & prod(wksData(k).size) == 1,
         wksString = [wksString, {sprintf('%-12s-%-3s', wksData(k).name, wksData(k).class)}];
      end
   end %for
   set(lb,'string',wksString);

   uicontrol('style','push',...
               'pos',[10, 10, 60, 20],...
               'string','OK',...
               'callback','wkstoenv(''OK'')');
   uicontrol('style','push',...
               'pos',[74, 10, 60, 20],...
               'string','Cancel',...
               'callback','close(gcf)');

   set(f2,'vis','on');
   setappdata(f2,'listbox',lb);

elseif strcmp(UIOperation, 'Select'),

   selType = get(gcf,'selectiontype');
   if strcmp(selType,'open'),
      wkstoenv('OK');
   end

elseif strcmp(UIOperation, 'OK'),

   f2 = gcf;
   f = get(f2,'userdata');

   lb = getappdata(f2,'listbox');
   lbValue = get(lb,'value');
   lbString = get(lb,'string');
   sel = lbString{lbValue};
   locdash = findstr('-', sel);
   selLTI = evalin('base', sel(1:locdash-1));

   bthan=get(f,'userdata');
   infmat=get(bthan(16),'userdata');
   QFTToolData = getappdata(f,'QFTToolData');

   old_cont=get(bthan(3),'userdata');
   old_lomat=get(bthan(1),'userdata');
   hint_bar = get(bthan(36),'userdata');

   set(bthan(27),'userdata',old_lomat);
   set(bthan(28),'userdata',old_cont);
   set(infmat(8,1),'enable','on');
   lomat=get(bthan(1),'userdata');

   T2=selLTI.Ts;

   q=1;
   if infmat(9,1)==2, q=[1;1]; end

   wl=lomat(1,:);
   nom=get(bthan(2),'userdata');
   delay=nom.ioDelay;
   T=nom.Ts;
   go_for_it=1;
   if T ~= T2,
      go_for_it=2;
   end

   cont_r = cntpars(selLTI);

   if go_for_it==1,
      cp=qcntbode(cont_r,wl,T);
      if length(nom) & infmat(9,1)==1,
         ncp=squeeze(freqresp(nom,wl)).';

      elseif infmat(9,1)==1,
         conL0=get(bthan(4),'userdata');
         L0=get(bthan(8),'userdata');
         cp0=qcntbode(conL0,wl,T); ncp=L0./cp0;

      elseif infmat(9,1)==2,
         ncp=get(bthan(17),'userdata');

      elseif infmat(9,1)==3,
         if length(nom),
            ncp=ones(1,length(wl));
         else
            ncp=get(bthan(8),'userdata');
         end
      end
      lomat(2:2+(infmat(9,1)==2),:)=cp(q,:).*ncp;
      set(bthan(3),'userdata',cont_r);
      set(bthan(1),'userdata',lomat);
      if infmat(9,1)==1, qnicplt(f);
      elseif infmat(9,1)==2, qmagplt(f);
      elseif infmat(9,1)==3, mgphplot(f);
      end

      [ControllerString,ListboxInfo] = cntstr(f,cont_r);
      set(QFTToolData.Elements(17), 'string', ControllerString,...
                                    'userdata', ListboxInfo);
      close(gcf);

   elseif go_for_it==2,
      if infmat(9,1)==1,
         errordlg(['Sampling time mismatch. Controller (',num2str(T2),') Environment (',num2str(T),')'],'Message','on');
      elseif infmat(9,1)==2,
         errordlg(['Sampling time mismatch. Pre-filter (',num2str(T2),') Environment (',num2str(T),')'],'Message','on');
      end
   end

end
