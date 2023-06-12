function qfrqchng(UIOperation)
% QFRQCHNG Set up frequency change interface. (Utility Function)
%          QFRQCHNG sets up the uicontrols to change the frequency vector
%          while still within the CAD environment.

% Author: Craig Borghesani
% 9/5/93
% Copyright (c) 2003, Terasoft, Inc.
%       $Revision: 1.5 $

if ~strcmp(UIOperation,'Setup'),
 f2=gcf;
 f = get(f2,'userdata');
else
 f = gcf;
end
bthan = get(f,'userdata');
infmat = get(bthan(16),'userdata');

T=get(bthan(13),'userdata');
hint_bar = get(bthan(36),'userdata');

proc_str=[];
if infmat(25,2)>1, proc_str=['(',int2str(infmat(25,2)),')']; end

switch UIOperation,
   case('Setup'),

      nom=get(bthan(2),'userdata');
      if length(nom),
         lomat=get(bthan(1),'userdata');
         w=lomat(1,:);
         first=num2str(w(1),3);
         last=num2str(w(length(w)),3);
         len=int2str(length(w));
         proc_num = int2str(infmat(25,2));
         f2 = findobj('tag',['qft2',proc_num]);

         if ~length(f2),
            ScreenSize     = get(0,'screensize');
            FigureWidth    = 410;
            FigureHeight   = 115;
            FigureLeft     = (ScreenSize(3) - FigureWidth)/2;
            FigureBottom   = (ScreenSize(4) - FigureHeight)/2;
            FigurePosition = [FigureLeft, FigureBottom, FigureWidth, FigureHeight];
            f2 = figure('name',['Frequency ',proc_str],...
                        'numbertitle','off',...
                        'position',FigurePosition,...
                        'menubar','none',...
                        'vis','off',...
                        'userdata',f,...
                        'tag',['qft2',proc_num],...
                        'windowstyle','modal',...
                        'closerequestfcn','qfrqchng(''Cancel'')',...
                        'handlevisibility','callback');
            infmat(12,1) = f2;

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

            SubFrame1Height = 60;

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
                              'string',' Frequency Limits',...
                              'pos',[SubFrame1UILeft, SubFrame1UIBottom, 100, 17],...
                              'horiz','left');

            SubFrame1UIBottom = SubFrame1UIBottom - 20;
            t(1)=uicontrol(f2,'style','text',...
                              'string','Minimum',...
                              'pos',[SubFrame1UILeft, SubFrame1UIBottom, 60, 17],...
                              'horiz','left');
            infmat(12,2)=uicontrol(f2,'style','edit',...
                                      'pos',[SubFrame1UILeft+60, SubFrame1UIBottom, 60, 20],...
                                      'background','w',...
                                      'string',first,...
                                      'horiz','right');
            t(2)=uicontrol(f2,'style','text',...
                              'string','Maximum',...
                              'pos',[SubFrame1UILeft+130, SubFrame1UIBottom, 60, 17],...
                              'horiz','left');
            infmat(12,3)=uicontrol(f2,'style','edit',...
                                      'pos',[SubFrame1UILeft+190, SubFrame1UIBottom, 60, 20],...
                                      'background','w',...
                                      'string',last,...
                                      'horiz','right');
            t(3)=uicontrol(f2,'style','text',...
                              'string','Length',...
                              'pos',[SubFrame1UILeft+260, SubFrame1UIBottom, 60, 17],...
                              'horiz','left');
            infmat(12,4)=uicontrol(f2,'style','edit',...
                                      'pos',[SubFrame1UILeft+320, SubFrame1UIBottom, 60, 20],...
                                      'background','w',...
                                      'string',len,...
                                      'horiz','right');

            SubFrame1UIBottom = SubFrame1UIBottom - 25;
            h(5)=uicontrol(f2,'style','check',...
                              'string','Pad Frequency Vector',...
                              'value',0,...
                              'pos',[SubFrame1UILeft, SubFrame1UIBottom, 200, 20]);

            h(1)=uicontrol(f2,'style','push','pos',[MainFrameLeft, 10, 60, 20],...
                              'string','Apply','callback','qfrqchng(''Apply'')');
            h(2)=uicontrol(f2,'style','push','pos',[MainFrameLeft+65, 10, 60, 20],...
                              'string','Cancel','callback','qfrqchng(''Cancel'')');
            h(3)=uicontrol(f2,'style','push','pos',[MainFrameLeft+131, 10, 60, 20],...
                              'string','UnDo','callback','qfrqchng(''UnDo'')');
            h(4)=uicontrol(f2,'style','push','pos',[MainFrameLeft+196, 10, 60, 20],...
                              'string','Done','callback','qfrqchng(''Done'')');
            set(h(3),'enable','off');
            set(infmat(12,2),'userdata',h);
            set(bthan(16),'userdata',infmat);
            set(infmat(12,1),'vis','on');

         else
            h=get(infmat(12,2),'userdata');
            set(h(3),'enable','off');
            set(infmat(12,2),'string',first);
            set(infmat(12,3),'string',last);
            set(infmat(12,4),'string',len);
            set(infmat(12,1),'vis','on');

         end
         set(h(1),'userdata',lomat);
         set(hint_bar,'string','Enter new frequency range');
         drawnow;

      else
         if infmat(9,1)==1,
            errordlg('Plant in complex format','Message','on');
         else
            errordlg('Initial loop in complex format','Message','on');
         end
      end

   case({'Apply', 'Done'}),

      h = get(infmat(12,2),'userdata');
      lomat = get(bthan(20),'userdata');
      cont = get(bthan(19),'userdata');

      if infmat(9,1)==3, axs=infmat(2,3:4);
      else axs=infmat(1,:); end

      if isempty(lomat),
         lomat = get(bthan(1),'userdata');
         cont = get(bthan(3),'userdata');
      end

%      if infmat(3,3) ~= 0,
%        wkw = lomat(1, infmat(3,3));
%     end

      T = get(bthan(13),'userdata');
      cur_fir = log10(lomat(1,1));
      cur_las = log10(lomat(1,length(lomat(1,:))));
      cur_len = length(lomat(1,:));
      if str2num(get(infmat(12,2),'string')) <= 0,
         set(infmat(12,2),'string',num2str(lomat(1,1),3));
         errordlg('Frequency minimum cannot be zero.', 'Input Error');
         return;
      end;%if
      fir=log10(str2num(get(infmat(12,2),'string')));
      las=log10(str2num(get(infmat(12,3),'string')));
      len=str2num(get(infmat(12,4),'string'));
      go_for_it = 1;
      if T > 0,
         if las > log10(pi/T),
            las = log10(pi/T);
         end
      end

      if go_for_it,
         set(h(2),'userdata',lomat);
         delay=infmat(10,1);
         nom=get(bthan(2),'userdata');
         wbs=get(bthan(11),'userdata');
         w=logspace(real(fir),real(las),len);
         w=sort([w,wbs]);
         w(find(diff(w)==0))=[];
         nomt=cntpars(nom);
         contnom=[cont;nomt];

      % pad frequency vector
         if get(h(5),'value'),
            for k=1:length(contnom(:,1)),
               if any(contnom(k,4)==[3,4]),
                  if contnom(k,1)<0.4,
                     w=qfrqenh(contnom(k,2),contnom(k,1),w,T);
                  end
               elseif contnom(k,4)==6,
                  ztas=contnom(k,1:2);
                  zta=ztas((ztas(1)>ztas(2))+1);
                  w=qfrqenh(contnom(k,3),zta,w,T);
               end
            end
            if T > 0,
               len_w = length(w);
               w_extra = logspace(log10(w(len_w-1)),log10(w(len_w)),20);
               w = sort([w,w_extra]);
            end
         end

% add in frequency from mouse movement
%         if infmat(3,3) ~= 0,
%            w = sort([w, wkw]);
%            infmat(3,3) = find(w == wkw);
%         end

         if infmat(9,1)==1,
            nlo=1;
            if length(nom),
               nlo=squeeze(freqresp(nom,w)).';
            end
%%            clo=qcntbode(cont,w,T).*exp(-i*w*delay);
            clo=qcntbode(cont,w,T);
            lo=nlo.*clo;

         else
            infmat(1,1:2)=[w(1),w(length(w))];
            infmat(2,1:2)=infmat(1,1:2);
            infmat(27,1:2)=infmat(1,1:2);
            lo=qcntbode(cont,w,T).*exp(-i*w*delay);

         end

% update location of red dot
%         if infmat(3,3) ~= 0,
%            set(infmat(6,1),'xdata',qfixfase(lo,axs,infmat(3,3)),...
%                            'ydata',20*log10(abs(lo(infmat(3,3)))),...
%                            'vis','on');
%         end

         set(infmat(6,1),'vis','off');
         infmat(3,3) = 0;
         lomat=[w;lo;ones(1,length(w))];
         if isempty(get(bthan(20),'userdata')),
            set(bthan(1),'userdata',lomat);
         else
            set(bthan(20),'userdata',lomat);
         end
         set(bthan(16),'userdata',infmat);
         if infmat(9,1)==1, qnicplt(f);
         elseif infmat(9,1)==3,
            set(infmat(24,1:2),'xlim',infmat(1,1:2));
            mgphplot(f);
         end
         if strcmp(UIOperation,'Done'),
            set(f2,'vis','off');
         else
            cur_firw = num2str(w(1),3);
            cur_lasw = num2str(w(length(w)),3);
            cur_lenw = int2str(length(w));
            set(infmat(12,2),'string',cur_firw);
            set(infmat(12,3),'string',cur_lasw);
            set(infmat(12,4),'string',cur_lenw);
            set(h(3),'enable','on');
         end

      else
         errordlg('Maximum frequency cannot be greater than pi/T','Message','on');
         set(infmat(12,3),'string',num2str(pi/T));
      end

   case({'Cancel', 'UnDo'}), % UnDo/Cancel

      h = get(infmat(12,2),'userdata');
      if strcmp(UIOperation,'Cancel'),
        lomat = get(h(1),'userdata');
      else
        lomat = get(h(2),'userdata');
      end
      w = lomat(1,:);
      las_fir = log10(w(1));
      las_las = log10(w(length(w)));
      las_len = length(w);
      fir=log10(str2num(get(infmat(12,2),'string')));
      las=log10(str2num(get(infmat(12,3),'string')));
      len=str2num(get(infmat(12,4),'string'));
      if any(abs([fir,las,len]-[las_fir,las_las,las_len])>0.01),
         if isempty(get(bthan(20),'userdata')),
            set(bthan(1),'userdata',lomat);
         else
            set(bthan(20),'userdata',lomat);
         end
         if infmat(9,1)==1, qnicplt(f);
         elseif infmat(9,1)==3,
            infmat(1,1:2)=[w(1),w(length(w))];
            infmat(2,1:2)=infmat(1,1:2);
            infmat(27,1:2)=infmat(1,1:2);
            set(infmat(24,1:2),'xlim',infmat(1,1:2));
            mgphplot(f);
         end
         set(bthan(16),'userdata',infmat);
      end
      if strcmp(UIOperation, 'Cancel'),
         set(f2,'vis','off');
      else
         set(infmat(12,2),'string',num2str(w(1),3));
         set(infmat(12,3),'string',num2str(w(length(w)),3));
         set(infmat(12,4),'string',int2str(length(w)));
         set(h(3),'enable','off');
      end

end
