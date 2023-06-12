function pfshpdef(prob,w,W,P,R,G,H,F0)
% PFSHPDEF Setup PFSHAPE. (Utility Function)
%          PFSHPDEF sets up the CAD environment for PFSHAPE.
%          Variables that are passed in as [] are set to their defaults.

% Author: Craig Borghesani
% Date: 9/6/93
% Revision: 2/16/96 1:11 PM V1.1 updates
% Copyright (c) 2003, Terasoft, Inc.
%       $Revision: 1.8 $

% load user defaults
defs = qftdefs;

T = F0.Ts;

if ~length(w),
   if T > 0,
      maxw = log10(pi/T);
   else
      maxw = defs(4,2);
   end
   w = logspace(defs(4,1),maxw,defs(4,3));
else
   w=w(:)';
   if T > 0,
      if w(end) > pi/T,
         w = w(find(w <= pi/T));
         w = [w,pi/T];
         loc_w = find([0,diff(w)==0]);
         w(loc_w) = [];
      end
   end
end

lw=length(w);
infmat=zeros(35,4);
if ~length(prob),
   error('Must define a problem');
end

if ~any(prob==[1 2 3 4 5 6 7 8 9]),
   error('Invalid problem type');
end

if ~length(P), P=1; end
if ~length(R), R=0; end
if ~length(G), G=1; end
if ~length(H), H=1; end
if prob==7, R=0; end
%zw=qsubset(wbd,w);
zero=find(w==0);
if length(zero),
   disp('Note: w=0 value being replaced with eps. Please do the same in workspace.');
   w(zero)=ones(1,length(zero))*eps;
end

%if length(nf)==1 & length(df)==1,
%   uF=nf/df;
%   if ~finite(uF),
%      error('Incorrect format: denominator polynomial cannot be zero');

%   elseif uF==0,
%      error('Incorrect format: numerator polynomial cannot be zero');

%   end
%end
cont=cntpars(F0);

%[rp,cp]=size(P);
%[rg,cg]=size(G);
%[rh,ch]=size(H);
%[rf,cf]=size(F);
%cm=max([cp cg ch cf]); v=ones(1,cm);

%%%%% V5 change for boolean indexing modification
%if cp == 1 & ~isa(P,'lti'), P=P(:,v); end
%if cg == 1 & ~isa(G,'lti'), G=G(:,v); end
%if ch == 1 & ~isa(H,'lti'), H=H(:,v); end

[jk,cl]=chksiso(prob,w,W,P,R,G,H);

if ~isempty(W),
   if isa(W(1),'lti'),
      if prob ~= 7,
         if length(w) == 1 & prod(size(W)) > 1,
            W = squeeze(freqresp(W, w));
         else
            W = squeeze(freqresp(W, w)).';
         end;%if
         if any(any(isnan(W))),
            error('Frequency vector for weight inconsistent with w.');
         end
         W = abs(W);
      else
         W1 = squeeze(freqresp(W(1), w)).';
         W2 = squeeze(freqresp(W(2), w)).';
         if any(any(isnan(W1))) | any(any(isnan(W2))),
            error('Frequency vector for weight inconsistent with w.');
         end
         W = [abs(W1);abs(W2)];
      end
   end %if isa(W,'lti')

   [rW,cW]=size(W);
   lw=length(w);
   if rW~=0,
      if rW==1,
         if cW==1,
            W=W*ones(2,lw);
         elseif cW==lw,
            W=W([1;1],:);
         end

      elseif rW==2 & cW==1,
         W=W(:,ones(2,lw));

      elseif rW>2,
         error('Incorrect weight vector format (max of 2 rows)');

      elseif cW~=lw & cW~=1,
         error('Weight vector incorrectly formed (length ~= freq vector)');

      end
   %   W=W(:,zw);

   end
   %w=w(zw); cl=cl(:,zw);
end

cls=size(cl);
if cls(1)>1,
   clmax=max(abs(cl));
   clmin=min(abs(cl));

else
   clmax=abs(cl);
   clmin=abs(cl);

end

if cls(2)==1,
   error('FSHAPE cannot be used with only one frequency');
end

sctlt={'Continuous-time Filter Shaping',['Discrete-time Filter Shaping (Ts = ',num2str(T),' sec)']};
proc_str='';
CADLength = length(findobj(allchild(0),'tag','CAD Window')) + 1;
infmat(25,2)=CADLength;
if CADLength > 1,
   proc_str=[' (',int2str(CADLength),')'];
end
dis=1+(T > 0);

% define various shades of grey
grey = get(0,'defaultuicontrolbackground');
ltgrey = [0.5,0.5,0.5]*1.5;
dkgrey = [0.5,0.5,0.5]*0.5;

% compute center of screen
ScreenSize     = get(0,'screensize');
FigureWidth    = ScreenSize(3) - 100;
FigureHeight   = ScreenSize(4) - 100;
FigureLeft     = (ScreenSize(3) - FigureWidth)/2;
FigureBottom   = (ScreenSize(4) - (FigureHeight+20))/2;
FigurePosition = [FigureLeft, FigureBottom, FigureWidth, FigureHeight];
dis=1 + (T > 0);
f = figure('position',FigurePosition,...
           'name',[sctlt{dis},proc_str],...
           'menubar','none',...
           'numbertitle','off',...
           'visible','off',...
           'tag','CAD Window',...
           'windowbuttonmotionfcn','qmouse(''Floating'')',...
           'defaultuicontrolbackgroundcolor',grey);

if datenum(version('-date')) < datenum('September 15, 2014'),
    xstr=['xitcade(',int2str(f),',0);'];
else
    xstr=['xitcade(',int2str(f.Number),',0);'];
end
set(f,'closerequestfcn',xstr);
set(f,'keypressfcn','qkeypress;');

% determine frame sizes depending upon figure size
BorderWidth = 5;
SectionBorderWidth = 3;
MainFrameWidth = 240;
MainFramePosition = [FigureWidth - (MainFrameWidth + 5), ...
                     5, ...
                     MainFrameWidth, ...
                     FigureHeight-10];
SubFramePosition = [MainFramePosition(1) + BorderWidth, ...
                    MainFramePosition(2) + BorderWidth, ...
                    MainFramePosition(3) - 2*BorderWidth, ...
                    FigureHeight - 2*BorderWidth];

SubFrame1Height = 80;
SubFrame2Height = (SubFramePosition(4) - SubFrame1Height) * 0.67 - 4*BorderWidth;
SubFrame3Height = (SubFramePosition(4) - SubFrame1Height) * 0.33 - 4*BorderWidth;

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

SubFrame3Position = [SubFramePosition(1), ...
                     SubFrame2Position(2) - (SubFrame3Height + 10), ...
                     SubFramePosition(3), ...
                     SubFrame3Height];
SubFrame3UILeft   = SubFrame3Position(1) + BorderWidth;
SubFrame3UIBottom = sum(SubFrame3Position([2,4])) - 10;

FullWidth         = SubFramePosition(3) - 2*BorderWidth;
HalfWidth         = FullWidth/2 - BorderWidth/2;
ThirdWidth        = FullWidth/3 - (BorderWidth*2)/3;

StatusBarPos = [5 + BorderWidth, ...
                10, ...
                FigureWidth - 2 * (5 + BorderWidth), ...
                20];

ft=qcntbode(cont,w,T);

if prob==7,
   clmxmn=[clmax; clmin];
else
   clmxmn=[clmax; clmax];
end

lo=clmxmn.*abs(ft([1;1],:));
a=[]; b=[];
if length(W),
   a=20*log10(W(1,:))';
   b=20*log10(W(2,:))';
   axs=[min(w),max(w),min([a(:);...
        b(:);20*log10(lo(:))])-5,max([a(:);b(:);20*log10(lo(:))])+5];
else
   axs=[min(w),max(w),min(20*log10(lo(:)))-5,max(20*log10(lo(:)))+5];
end

ax = axes('box','on',...
          'xgrid','on',...
          'ygrid','on',...
          'units','pixel',...
          'gridlinestyle',':',...
          'pos',[60, 70, FigureWidth - MainFrameWidth - 100, FigureHeight - 100],...
          'nextplot','add',...
          'xlim',axs(1:2),'ylim',axs(3:4),'xscale','log');

if length(a),
   semilogx(w,a,'r--',w,b,'r--');
end

xlabel('Frequency (rad/sec)');
ylabel('Magnitude (dB)');

scrnp=get(0,'screenp');
if scrnp > 80, fnt_size = 8;
else fnt_size = 10; end

comp_type = computer;
if length(comp_type) < 3,
   comp_type = 'NUN';
end

infmat(1,:)=axs;
infmat(9,1)=2;
infmat(23,1)=strcmp('VAX',comp_type(1:3));
infmat(24,1)=ax;
infmat(26,:)=axs;
infmat(27,:)=axs;

bt=qbutnset([],3+(T > 0),1,[],f);
set(f,'userdata',bt);
infmat(8,1)=bt(9);
set(bt(9),'enable','off');

infmat(6,1)=line('xdata',-180,'ydata',0,'marker','o','color','g',...
                 'markersize',8,'vis','off',...
                 'markerfacecolor','r');

if prob~=7,
   vl(1) = line('xdata',0,'ydata',0,'linestyle','-','color','b');
   vl2(1) = line('xdata',0,'ydata',0,'linestyle','--','vis','off','color','b');
else
   vl(1) = line('xdata',0,'ydata',0,'linestyle','-','color','b');
   vl(2) = line('xdata',0,'ydata',0,'linestyle','-','color','g');
   vl2(1) = line('xdata',0,'ydata',0,'linestyle','--','vis','off','color','b');
   vl2(2) = line('xdata',0,'ydata',0,'linestyle','--','color','g','vis','off');
end

lomat=[w;lo;ones(1,length(lo))];

%%% Control panel for controller design
QFTToolData.Frame = uicontrol('style','push',...
                              'pos',MainFramePosition,...
                              'enable','off');

%% Pointer Information (subframe 1)
QFTToolData.Mouse(1) = uicontrol('style','frame',...
                                 'pos',SubFrame1Position,...
                                 'fore','k');
QFTToolData.Mouse(2) = uicontrol('style','text',...
                                 'pos',[SubFrame1UILeft, SubFrame1UIBottom, 100, 17],...
                                 'string',' Pointer Info',...
                                 'horiz','left');

SubFrame1UIBottom = SubFrame1UIBottom - 15;
QFTToolData.Mouse(3) = uicontrol('style','text',...
                                 'pos',[SubFrame1UILeft, SubFrame1UIBottom, FullWidth, 15],...
                                 'string','',...
                                 'horiz','left');

SubFrame1UIBottom = SubFrame1UIBottom - 15;
QFTToolData.Mouse(4) = uicontrol('style','text',...
                                 'pos',[SubFrame1UILeft, SubFrame1UIBottom, FullWidth, 15],...
                                 'string','',...
                                 'horiz','left');

SubFrame1UIBottom = SubFrame1UIBottom - 15;
QFTToolData.Mouse(5) = uicontrol('style','text',...
                                 'pos',[SubFrame1UILeft, SubFrame1UIBottom, FullWidth, 15],...
                                 'string','',...
                                 'horiz','left');

SubFrame1UIBottom = SubFrame1UIBottom - 20;
QFTToolData.Mouse(6) = uicontrol('style','radio',...
                                 'pos',[SubFrame1UILeft, SubFrame1UIBottom, HalfWidth, 20],...
                                 'string','Rad/Sec',...
                                 'value',1,...
                                 'horiz','left',...
                                 'callback','qunitchng',...
                                 'enable','off');
QFTToolData.Mouse(7) = uicontrol('style','radio',...
                                 'pos',[SubFrame1UILeft+HalfWidth+5, SubFrame1UIBottom, HalfWidth, 20],...
                                 'string','Hertz',...
                                 'value',0,...
                                 'horiz','left',...
                                 'callback','qunitchng',...
                                 'enable','off');

set(QFTToolData.Mouse(6),'userdata',QFTToolData.Mouse(7));
set(QFTToolData.Mouse(7),'userdata',QFTToolData.Mouse(6));

infmat(30,1:3)=QFTToolData.Mouse(3:5);

%% Controller design (subframe 2)
QFTToolData.Elements(1) = uicontrol('style','frame',...
                                  'pos',SubFrame2Position,...
                                  'fore','k');
QFTToolData.Elements(2) = uicontrol('style','text',...
                                  'pos',[SubFrame2UILeft, SubFrame2UIBottom, 100, 17],...
                                  'string',' Filter Elements',...
                                  'horiz','left');

SubFrame2UIBottom = SubFrame2UIBottom - 15;
if T == 0,
   ele21String = 'Select elements to Add.';
else
   ele21String = ['Select elements to Add. (Ts = ', num2str(T), 'sec)'];
end
QFTToolData.Elements(21) = uicontrol('style','text',...
                                     'pos',[SubFrame2UILeft, SubFrame2UIBottom, FullWidth, 15],...
                                     'string',ele21String,...
                                     'horiz','left');

if T == 0,
   ElementPopupString = {'Gain',...
                         'Real Pole', ...
                         'Real Zero', ...
                         'Complex Pole', ...
                         'Complex Zero', ...
                         'Notch',...
                         'Integrator/Differentiator'};

else
   ElementPopupString = {'Gain',...
                         'Real Pole', ...
                         'Real Zero', ...
                         'Complex Pole', ...
                         'Complex Zero', ...
                         'Notch',...
                         'Integrator/Differentiator', ...
                         'Delay/Predictor'};

end

SubFrame2UIBottom = SubFrame2UIBottom - 20;
QFTToolData.Elements(3) = uicontrol('style','popup',...
                                  'pos',[SubFrame2UILeft, SubFrame2UIBottom, FullWidth, 20],...
                                  'string',ElementPopupString, ...
                                  'callback','qelmtpopup(''New'')',...
                                  'back','w',...
                                  'horiz','left',...
                                  'value',1);

infmat(16,4) = QFTToolData.Elements(3);

SubFrame2UIBottom = SubFrame2UIBottom - 17;
QFTToolData.Elements(4) = uicontrol('style','text',...
                                  'pos',[SubFrame2UILeft, SubFrame2UIBottom, ThirdWidth, 15],...
                                  'string','pole',...
                                  'horiz','left');
QFTToolData.Elements(5) = uicontrol('style','text',...
                                  'pos',[SubFrame2UILeft+ThirdWidth+5, SubFrame2UIBottom, ThirdWidth, 15],...
                                  'string','',...
                                  'horiz','left');
QFTToolData.Elements(6) = uicontrol('style','text',...
                                  'pos',[SubFrame2UILeft+ThirdWidth*2+5*2, SubFrame2UIBottom, ThirdWidth, 15],...
                                  'string','',...
                                  'horiz','left');
infmat(19,1:3) = QFTToolData.Elements(4:6);

SubFrame2UIBottom = SubFrame2UIBottom - 20;
QFTToolData.Elements(7) = uicontrol('style','edit',...
                                  'pos',[SubFrame2UILeft, SubFrame2UIBottom, ThirdWidth-10, 20],...
                                  'string','',...
                                  'back','w',...
                                  'horiz','right');
QFTToolData.Elements(10) = uicontrol('style','slider',...
                                  'pos',[SubFrame2UILeft+ThirdWidth-10, SubFrame2UIBottom, 10, 20]);

QFTToolData.Elements(8) = uicontrol('style','edit',...
                                  'pos',[SubFrame2UILeft+ThirdWidth+5, SubFrame2UIBottom, ThirdWidth-10, 20],...
                                  'string','',...
                                  'back','w',...
                                  'horiz','right');
QFTToolData.Elements(11) = uicontrol('style','slider',...
                                  'pos',[SubFrame2UILeft+ThirdWidth*2+5-10, SubFrame2UIBottom, 10, 20]);

QFTToolData.Elements(9) = uicontrol('style','edit',...
                                  'pos',[SubFrame2UILeft+ThirdWidth*2+5*2, SubFrame2UIBottom, ThirdWidth-10, 20],...
                                  'string','',...
                                  'back','w',...
                                  'horiz','right');
QFTToolData.Elements(12) = uicontrol('style','slider',...
                                  'pos',[SubFrame2UILeft+ThirdWidth*3+5*2-10, SubFrame2UIBottom, 10, 20]);

infmat(16,1:3) = QFTToolData.Elements(7:9);
infmat(17,1:3) = QFTToolData.Elements(10:12);

SubFrame2UIBottom = SubFrame2UIBottom - 25;
QFTToolData.Elements(13) = uicontrol('style','push',...
                                  'pos',[SubFrame2UILeft, SubFrame2UIBottom, FullWidth, 20],...
                                  'string','Add using Input Fields');
infmat(13,1) = QFTToolData.Elements(13);

%SubFrame2UIBottom = SubFrame2UIBottom - 25;
QFTToolData.Elements(16) = uicontrol('style','push',...
                                  'pos',[SubFrame2UILeft, SubFrame2UIBottom, FullWidth, 20],...
                                  'string','Add using Mouse','vis','off');
infmat(13,4) = QFTToolData.Elements(16);

SubFrame2UIBottom = SubFrame2UIBottom - 25;
QFTToolData.Elements(14) = uicontrol('style','push',...
                                  'pos',[SubFrame2UILeft, SubFrame2UIBottom, HalfWidth, 20],...
                                  'string','Apply',...
                                  'callback','qelmts(-1,0,''Done'')');
infmat(13,2) = QFTToolData.Elements(14);

QFTToolData.Elements(15) = uicontrol('style','push',...
                                  'pos',[SubFrame2UILeft+HalfWidth+5, SubFrame2UIBottom, HalfWidth, 20],...
                                  'string','Cancel',...
                                  'callback','qelmts(-1,0,''Cancel'')');
infmat(13,3) = QFTToolData.Elements(15);

SubFrame2UIBottom = SubFrame2UIBottom - 20;
QFTToolData.Elements(23) = uicontrol('style','text',...
                                     'pos',[SubFrame2UILeft, SubFrame2UIBottom, FullWidth, 15],...
                                     'string','Select elements to Tune.',...
                                     'horiz','left');

% setup uicontextmenu for listbox
UICont = uicontextmenu('callback','qlbxcmenu(''Check'')');
uimenu(UICont,'label','Select All','callback','qlbxcmenu(''Select All'')');
UICont2 = uimenu(UICont,'label','Sort','separator','on','vis','off');
uimenu(UICont2,'label','By Type','callback','qelesort(''ByType'')','checked','on','tag','sortByType','vis','off');
uimenu(UICont2,'label','All','callback','qelesort(''All'')','tag','sortAll','vis','off');

ListboxBottom = SubFrame2Position(2) + 30;
ListboxHeight = (SubFrame2UIBottom) - ListboxBottom;
QFTToolData.Elements(17) = uicontrol('style','listbox',...
                                  'pos',[SubFrame2UILeft, ListboxBottom, FullWidth, ListboxHeight],...
                                  'back','w',...
                                  'string','',...
                                  'callback','qelmtlistbox(''Edit'')',...
                                  'max',2,...
                                  'uicontextmenu',UICont);

SubFrame2UIBottom = SubFrame2Position(2) + 5;
QFTToolData.Elements(18) = uicontrol('style','push',...
                                  'pos',[SubFrame2UILeft, SubFrame2UIBottom, HalfWidth, 20],...
                                  'string','Delete',...
                                  'callback','qdelmts');
QFTToolData.Elements(19) = uicontrol('style','push',...
                                  'pos',[SubFrame2UILeft+HalfWidth+5, SubFrame2UIBottom, HalfWidth, 20],...
                                  'string','Reduction',...
                                  'callback','qrelmts(''Setup'')');
QFTToolData.Elements(20) = uicontrol('style','push',...
                                  'pos',[SubFrame2UILeft+ThirdWidth*2+5*2, SubFrame2UIBottom, ThirdWidth, 20],...
                                  'string','Convert',...
                                  'callback','qcelmts',...
                                  'vis','off');

%% Bound (subframe 3)
QFTToolData.Bounds(1) = uicontrol('style','frame',...
                                  'pos',SubFrame3Position,...
                                  'fore','k',...
                                  'enable','off');
QFTToolData.Bounds(2) = uicontrol('style','text',...
                                  'pos',[SubFrame3UILeft, SubFrame3UIBottom, 100, 17],...
                                  'string',' Bounds',...
                                  'horiz','left',...
                                  'enable','off');

ListboxBottom = SubFrame3Position(2) + 30;
ListboxHeight = (SubFrame3UIBottom) - ListboxBottom;
QFTToolData.Bounds(3) = uicontrol('style','listbox',...
                                  'pos',[SubFrame3UILeft, ListboxBottom, FullWidth, ListboxHeight],...
                                  'back','w',...
                                  'max',2,...
                                  'enable','off');

SubFrame3UIBottom = SubFrame3Position(2) + 5;
QFTToolData.Bounds(4) = uicontrol('style','push',...
                                  'pos',[SubFrame3UILeft, SubFrame3UIBottom, ThirdWidth, 20],...
                                  'string','On/Off',...
                                  'callback','bndonoff(''On/Off'')',...
                                  'enable','off');
QFTToolData.Bounds(5) = uicontrol('style','push',...
                                  'pos',[SubFrame3UILeft+ThirdWidth+5, SubFrame3UIBottom, ThirdWidth, 20],...
                                  'string','All On',...
                                  'callback','bndonoff(''All On'')',...
                                  'enable','off');
QFTToolData.Bounds(6) = uicontrol('style','push',...
                                  'pos',[SubFrame3UILeft+ThirdWidth*2+5*2, SubFrame3UIBottom, ThirdWidth, 20],...
                                  'string','All Off',...
                                  'callback','bndonoff(''All Off'')',...
                                  'enable','off');

%setup hint bar
QFTToolData.HintBar = uicontrol(f,'style','edit','pos',[0,0,FigureWidth-MainFrameWidth-10,20],...
                     'horizontalalignment','left','enable','inactive');

setappdata(f,'QFTToolData',QFTToolData);

set(bt(1),'userdata',lomat);
set(bt(2),'userdata',cl);
set(bt(3),'userdata',cont);
set(bt(4),'userdata',cont);
set(bt(6),'userdata',w);
set(bt(8),'userdata',lo);
set(bt(10),'userdata',vl);
set(bt(13),'userdata',T);
set(bt(16),'userdata',infmat);
set(bt(17),'userdata',clmxmn);
set(bt(18),'userdata',[a;b]);
set(bt(19),'userdata',[]);
set(bt(20),'userdata',[]);
set(bt(21),'userdata',vl2);
set(bt(36),'userdata',QFTToolData.HintBar);

[ControllerString,ListboxInfo] = cntstr(f,cont);
set(QFTToolData.Elements(17), 'string', ControllerString,...
                              'userdata', ListboxInfo);

qmagplt(f);
qelmtpopup;
drawnow;

set(QFTToolData.HintBar,'string','Ready');
set(f,'vis','on', 'resizefcn', 'qresize',...
            'handlevisibility','callback');

