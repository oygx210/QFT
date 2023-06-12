function lpshpdef(w,bdb,P0,C0,phase)
% LPSHPDEF Set up LPSHAPE/DLPSHAPE environments. (Utility Function)
%          LPSHPDEF sets up the CAD environment for LPSHAPE/DLPSHAPE.
%          Variables that are passed in as [] are set to their defaults.

% Author: Craig Borghesani
% Date: 9/6/93
% Revised: 2/16/96 1:06 PM V1.1 updates
% Copyright (c) 2003, Terasoft, Inc.

% load user defaults
defs = qftdefs;

delay = P0.ioDelay;
T = P0.Ts;

if isempty(w) && ~isa(P0, 'frd'),
   if T > 0,
      maxw = log10(pi/T);
   else
      maxw = defs(4,2);
   end
   w = logspace(defs(4,1),maxw,defs(4,3));
elseif ~isa(P0, 'frd'),
   w=w(:)';
   if T > 0,
      if w(length(w)) > pi/T,
         w = w(find(w <= pi/T));
         w = [w,pi/T];
         loc_w = find([0,diff(w)==0]);
         w(loc_w) = [];
      end
   end
else
    w = P0.Frequency(:)';
end

if ~length(phase),
   phase = defs(1,1):defs(1,2):defs(1,3);
end

lw=length(w);
%[rm0,cm0]=size(uL0);
%[rp0,cp0]=size(vL0);

%if rm0>1 & rp0==0,
%   uL0=conj(uL0');
%   [rm0,cm0]=size(uL0);

%elseif rm0>1 | rp0>1,
%   error('Num/Den information must be in row format');

%end

%if cm0~=lw & rp0==0 & rm0~=0,
%   error('Frequency vector does not match complex data');

%elseif cm0==lw,
%%%%%%% V1.1 change: removal of this message
%% disp('Program is assuming data in complex format'); pause(1);
%end

%if (~length(numC0) & length(denC0)) | (~length(denC0) & length(numC0)),
%   error('Initial controller cannot be in complex format');
%end

%if ~length(uL0),
%   uL0=1;
%end

%if ~length(vL0),
%   vL0=1;
%end

%if ~length(numC0),
%   numC0=1;
%end

%if ~length(denC0),
%   denC0=1;
%end

if length(bdb),
   [bdaxs,wbs]=qfindinf(phase,bdb,1);
   w=sort([w wbs]);
end
%nom=[];

%if cm0~=lw,
%   if cp0 ~= 0 & cm0 ~= 0,
%      nom=[zeros(1,cp0-cm0) uL0;zeros(1,cm0-cp0) vL0]
%   else
%      nom=[uL0; vL0];
%   end

%   uL0=qcpqft(nom(1,:),nom(2,:),w,T);
%end

con_ar = cntpars(C0);
uL0=[squeeze(freqresp(P0, w)).'].*qcntbode(con_ar,w,T);

if length(bdb),
   [rb,cb]=size(bdb);
   bdbw=bdb(:,find(bdb(rb-1,:)==wbs(1)));
   maxbd=qfindinf(phase,bdbw,1);
   magbd=10^(maxbd(4)/20);
   z=find(w>=wbs(1)); loc=z(1);
%   k1=magbd/abs(uL0(loc));
%   if con_ar(1,1) == 1,
%      con_ar(1,1) = k1;
%   end
else
   wbs=[]; k1=[];
   bdaxs=[1/eps eps 1/eps eps];
end

%if length(numC0)==1 & length(denC0)==1,
%   m=numC0/denC0;
%   if ~finite(m),
%      error('Incorrect format: denominator polynomial cannot be zero');

%   elseif m==0,
%      error('Incorrect format: numerator polynomial cannot be zero');
%   end

%%%%%%% V1.1 change: removal of auto gain setting
%% if length(k1) & m==1, m=k1; end

%   con_ar(1,1:4)=[m NaN NaN 0];
%   if nargin==8,
%      con_ar(2,1:4)=[0 NaN NaN 0.7];
%   else
%      con_ar(2,1:4)=[0 0 NaN 0.5];
%      con_ar(3,1:4)=[0 NaN NaN 0.6];
%   end

%else
%   con_ar=cntpars(numC0(:)',denC0(:)',T);

%end

if any(isnan(uL0)),
   error('Frequency vector and frequency data for nominal plant must be the same.');
end

axs=[-360,0,min([bdaxs(3),min(20*log10(abs(uL0)))])-5,...
            max([bdaxs(4),max(20*log10(abs(uL0)))])+5];

if length(bdb),
   [coora,coorb]=wherebnd(bdb);
end

% determine figure title
infmat=zeros(35,4);
chil=get(0,'children');
sctlt={'Continuous-time Loop Shaping',['Discrete-time Loop Shaping (Ts = ',num2str(T),' sec)']};
proc_str=[];
len=[28,26];
CADLength = length(findobj(allchild(0),'tag','CAD Window')) + 1;
infmat(25,2)=CADLength;
if CADLength > 1,
   proc_str=[' (',int2str(CADLength),')'];
end

% setup CAD figure window

% define various shades of grey
grey = get(0,'defaultuicontrolbackground');
ltgrey = [0.5,0.5,0.5]*1.5;
dkgrey = [0.5,0.5,0.5]*0.5;

% make it so the toolbar stays with the figure
%DefaultFigureToolbarStatus = get(0,'defaultfiguretoolbar');
%set(0,'defaultfiguretoolbar','figure');

% compute center of screen
ScreenSize     = get(0,'screensize');
FigureWidth    = ScreenSize(3) - 100;
FigureHeight   = ScreenSize(4) - 100;
FigureLeft     = (ScreenSize(3) - FigureWidth)/2;
FigureBottom   = (ScreenSize(4) - (FigureHeight+20))/2;
FigurePosition = [FigureLeft, FigureBottom, FigureWidth, FigureHeight];
dis=1+(T > 0);
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
%set(0,'defaultfiguretoolbar',DefaultFigureToolbarStatus);

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

XTickValues = qaxesadjust(axs(1:2));

% setup axis and labels
a=axes('box','on',...
       'xgrid','on',...
       'ygrid','on',...
       'units','pixel',...
       'pos',[60, 70, FigureWidth - MainFrameWidth - 100, FigureHeight - 100],...
       'nextplot','add',...
       'xlim',axs(1:2),...
       'xtick', XTickValues, ...
       'ylim',axs(3:4));

xlabstr='Open-Loop Phase (deg)'; ylabstr='Open-Loop Gain (dB)';
xlabel(xlabstr);
ylabel(ylabstr);

b=[];
if length(bdb),
   b=qplotbd(phase,bdb,coora,coorb,axs);
end
[rb,cb]=size(b);

scrnp=get(0,'screenp');
if scrnp > 80,
   fnt_size = 8;
else
   fnt_size = 10;
end

comp_type = computer;
if length(comp_type) < 3,
   comp_type = 'NUN';
end

infmat(1,:)=axs;
infmat(2,:)=bdaxs;
infmat(9,1)=1;
infmat(10,1)=delay;
infmat(23,1)=strcmp('VAX',comp_type(1:3));
infmat(24,1)=a;
infmat(26,:)=axs;
infmat(27,1:2)=axs(1:2);
infmat(28,2:3)=axs(1:2);

[bt,bnd_bt]=qbutnset(b,1+(T > 0),delay,P0,f);
set(f,'userdata',bt);
infmat(8,1)=bt(9);
set(bt(9),'enable','off');

infmat(6,1)=line(-180,0,'marker','o','color','g',...
                 'markersize',8,'vis','off',...
                 'markerfacecolor','r');

% setup loop response and frequency markers
vo1=[];
vo2 = [];
v1 = line(0,0,'linestyle','-');
v2 = line(0,0,'linestyle',':','vis','off');
cvec=['r';'g';'b';'c';'m'];
clr=[cvec;cvec;cvec;cvec];
clr=[clr;clr;clr;clr;clr];
clr=[clr;clr;clr;clr;clr];
for j=1:length(wbs),
   vo1(j)=line(0,0,'color',clr(j),...
               'marker','o','markersize',6);
   vo2(j)=line(0,0,'color',clr(j),...
               'marker','o','markersize',6,'vis','off');
end

lomat=[w;uL0;ones(1,length(uL0))];

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
                                 'callback','qunitchng;');
QFTToolData.Mouse(7) = uicontrol('style','radio',...
                                 'pos',[SubFrame1UILeft+HalfWidth+5, SubFrame1UIBottom, HalfWidth, 20],...
                                 'string','Hertz',...
                                 'value',0,...
                                 'horiz','left',...
                                 'callback','qunitchng;');
set(QFTToolData.Mouse(6),'userdata',QFTToolData.Mouse(7));
set(QFTToolData.Mouse(7),'userdata',QFTToolData.Mouse(6));

infmat(30,1:3)=QFTToolData.Mouse(3:5);

%% Controller design (subframe 2)
QFTToolData.Elements(1) = uicontrol('style','frame',...
                                  'pos',SubFrame2Position,...
                                  'fore','k');
QFTToolData.Elements(2) = uicontrol('style','text',...
                                  'pos',[SubFrame2UILeft, SubFrame2UIBottom, 150, 17],...
                                  'string',' Controller Elements',...
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
                         'Lead or Lag', ...
                         'Notch',...
                         'Super 2nd', ...
                         'Complex Lead/Lag', ...
                         'Integrator/Differentiator'};

else
   ElementPopupString = {'Gain',...
                         'Real Pole', ...
                         'Real Zero', ...
                         'Complex Pole', ...
                         'Complex Zero', ...
                         'Lead or Lag', ...
                         'Notch',...
                         'Super 2nd', ...
                         'Complex Lead/Lag', ...
                         'Integrator/Differentiator', ...
                         'Delay/Predictor'};

end

SubFrame2UIBottom = SubFrame2UIBottom - 20;
QFTToolData.Elements(3) = uicontrol('style','popup',...
                                  'pos',[SubFrame2UILeft, SubFrame2UIBottom, FullWidth, 20],...
                                  'string', ElementPopupString, ...
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
                                  'string','Add using Mouse',...
                                  'vis','off');
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
                                  'enable','off',...
                                  'callback','bndonoff(''On/Off'')');

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

if ~isempty(bnd_bt),
   set(QFTToolData.Bounds(3), 'string',bnd_bt(:,1),...
                              'userdata',bnd_bt(:,2));
   set(QFTToolData.Bounds(1:6),'enable','on');
end

setappdata(f,'QFTToolData',QFTToolData);

set(bt(1),'userdata',lomat);
set(bt(2),'userdata',P0);
set(bt(3),'userdata',con_ar);
set(bt(4),'userdata',con_ar);
set(bt(6),'userdata',w);
set(bt(8),'userdata',uL0);
set(bt(10),'userdata',v1);
set(bt(11),'userdata',wbs);
set(bt(12),'userdata',bnd_bt);
set(bt(13),'userdata',T);
set(bt(16),'userdata',infmat);
set(bt(17),'userdata',vo1);
set(bt(19),'userdata',[]);
set(bt(20),'userdata',[]);
set(bt(21),'userdata',v2);
set(bt(22),'userdata',vo2);
set(bt(31),'userdata',{[], []});
set(bt(32),'userdata',bdb);
set(bt(33),'userdata',phase);
set(bt(36),'userdata',QFTToolData.HintBar);

[ControllerString,ListboxInfo] = cntstr(f,con_ar);
set(QFTToolData.Elements(17), 'string', ControllerString,...
                              'userdata', ListboxInfo);

qnicplt(f);
qelmtpopup;
drawnow;

set(QFTToolData.HintBar,'string','Ready');
set(f,'vis','on', 'resizefcn','qresize',...
      'handlevisibility','callback');

