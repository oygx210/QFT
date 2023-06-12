function wht=cntdisp(f,cont,Name)
% CNTDISP Controller display. (Utility Function)
%         CNTDISP displays the present controller.  It also sets up the
%         window for Edit, Iterate, Delete, Conversion and Reduction.

% Author: Craig Borghesani
% 10/10/93
% Copyright (c) 2003, Terasoft, Inc.

bthan=get(f,'userdata');
infmat=get(bthan(16),'userdata');
T=get(bthan(13),'userdata');
hint_bar = get(bthan(36),'userdata');

flag3 = infmat(9,1);
delay = infmat(10,1);
proc_str = [];
proc_num = int2str(infmat(25,2));
if infmat(25,2)>1, proc_str=[' (',int2str(infmat(25,2)),')']; end

ScreenSize     = get(0,'screensize');
FigureWidth    = 300;
FigureHeight   = 250;
FigureLeft     = (ScreenSize(3) - FigureWidth)/2;
FigureBottom   = (ScreenSize(4) - (FigureHeight+20))/2;
FigurePosition = [FigureLeft, FigureBottom, FigureWidth, FigureHeight];
f2 = figure('name',[Name, ' ', proc_str],...
            'numbertitle','off',...
            'position',FigurePosition,...
            'menubar','none',...
            'vis','off',...
            'userdata',f,...
            'tag',['qft3',proc_num]);

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

FullWidth         = SubFramePosition(3) - 2*BorderWidth;

MainFrameLeft = MainFramePosition(1) + 5;
MainFrameRight = FigureWidth - 2*BorderWidth;

SubFrame1Height = MainFramePosition(4) - 20;

SubFrame1Position = [SubFramePosition(1), ...
                     MainFramePosition(4) - (SubFrame1Height + 5), ...
                     SubFramePosition(3), ...
                     SubFrame1Height];
SubFrame1UILeft   = SubFrame1Position(1) + BorderWidth;
SubFrame1UIBottom = sum(SubFrame1Position([2,4])) - 10;

uicontrol('style','push',...
          'pos',MainFramePosition,...
          'enable','inactive');
uicontrol('style','frame',...
          'pos',SubFrame1Position,...
          'fore','k');
uicontrol('style','text',...
          'pos',[SubFrame1UILeft, SubFrame1UIBottom, 100, 17],...
          'string',[' ', Name, ' Elements'],...
          'horiz','left');

ListboxBottom = 20;
ListboxHeight = SubFrame1UIBottom - ListboxBottom;

[NameString] = cntstr(f,cont);

if ~isempty(delay) & delay ~= 0,
   NameString = [{['Delay: ', num2str(delay)]}, NameString];
end;%if

if ~isempty(T) & T ~= 0,
   NameString = [NameString, {['Sampling Time: ', num2str(T), ' sec']}];
end;%if

uicontrol('style','listbox',...
          'pos',[SubFrame1UILeft, ListboxBottom, FullWidth, ListboxHeight],...
          'back','w',...
          'string',NameString);

set(f2,'vis','on',...
       'handlevisibility','callback');

