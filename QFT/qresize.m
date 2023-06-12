function qresize
%QRESIZE Readjust uicontrols depending upon figure size. (Utility Function)

% Author: Craig Borghesani <cborg@terasoft.com>
% Date: 2/1/02 9:37PM
% Copyright (c) 2002, Terasoft, Inc.

f = gcf;
QFTToolData = getappdata(f,'QFTToolData');

FigurePosition = get(f,'position');
FigureLeft = FigurePosition(1);
FigureBottom = FigurePosition(2);
FigureWidth = FigurePosition(3);
FigureHeight = FigurePosition(4);

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

% perform test of listbox size
SubFrame2UIBottomTest = SubFrame2UIBottom - 15;
SubFrame2UIBottomTest = SubFrame2UIBottomTest - 20;
SubFrame2UIBottomTest = SubFrame2UIBottomTest - 17;
SubFrame2UIBottomTest = SubFrame2UIBottomTest - 20;
SubFrame2UIBottomTest = SubFrame2UIBottomTest - 25;
SubFrame2UIBottomTest = SubFrame2UIBottomTest - 25;
SubFrame2UIBottomTest = SubFrame2UIBottomTest - 20;
ListboxBottomTest = SubFrame2Position(2) + 30;
ListboxHeightTest = (SubFrame2UIBottomTest) - ListboxBottomTest;
if ListboxHeightTest < 20, % need to reset figure size
   FigureHeight = 420;
   set(f,'position',[FigureLeft, FigureBottom, FigureWidth, FigureHeight]);

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
end

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

ax = findobj(f,'type','axes');
set(ax,'position',[60, 70, FigureWidth - MainFrameWidth - 100, FigureHeight - 100]);

%%% Control panel for controller design
set(QFTToolData.Frame, 'pos', MainFramePosition);

%% Pointer Information (subframe 1)
set(QFTToolData.Mouse(1),'pos', SubFrame1Position);
set(QFTToolData.Mouse(2),'pos',[SubFrame1UILeft, SubFrame1UIBottom, 100, 17]);

SubFrame1UIBottom = SubFrame1UIBottom - 15;
set(QFTToolData.Mouse(3),'pos',[SubFrame1UILeft, SubFrame1UIBottom, FullWidth, 15]);

SubFrame1UIBottom = SubFrame1UIBottom - 15;
set(QFTToolData.Mouse(4),'pos',[SubFrame1UILeft, SubFrame1UIBottom, FullWidth, 15]);

SubFrame1UIBottom = SubFrame1UIBottom - 15;
set(QFTToolData.Mouse(5),'pos',[SubFrame1UILeft, SubFrame1UIBottom, FullWidth, 15]);

SubFrame1UIBottom = SubFrame1UIBottom - 20;
set(QFTToolData.Mouse(6),'pos',[SubFrame1UILeft, SubFrame1UIBottom, HalfWidth, 20]);
set(QFTToolData.Mouse(7),'pos',[SubFrame1UILeft+HalfWidth+5, SubFrame1UIBottom, HalfWidth, 20]);

%% Controller design (subframe 2)
set(QFTToolData.Elements(1),'pos',SubFrame2Position);
set(QFTToolData.Elements(2),'pos',[SubFrame2UILeft, SubFrame2UIBottom, 100, 17]);

SubFrame2UIBottom = SubFrame2UIBottom - 15;
set(QFTToolData.Elements(21),'pos',[SubFrame2UILeft, SubFrame2UIBottom, FullWidth, 15]);

SubFrame2UIBottom = SubFrame2UIBottom - 20;
set(QFTToolData.Elements(3),'pos',[SubFrame2UILeft, SubFrame2UIBottom, FullWidth, 20]);

SubFrame2UIBottom = SubFrame2UIBottom - 17;
set(QFTToolData.Elements(4),'pos',[SubFrame2UILeft, SubFrame2UIBottom, ThirdWidth, 15]);
set(QFTToolData.Elements(5),'pos',[SubFrame2UILeft+ThirdWidth+5, SubFrame2UIBottom, ThirdWidth, 15]);
set(QFTToolData.Elements(6),'pos',[SubFrame2UILeft+ThirdWidth*2+5*2, SubFrame2UIBottom, ThirdWidth, 15]);

SubFrame2UIBottom = SubFrame2UIBottom - 20;
set(QFTToolData.Elements(7),'pos',[SubFrame2UILeft, SubFrame2UIBottom, ThirdWidth-10, 20]);
set(QFTToolData.Elements(10),'pos',[SubFrame2UILeft+ThirdWidth-10, SubFrame2UIBottom, 10, 20]);

set(QFTToolData.Elements(8),'pos',[SubFrame2UILeft+ThirdWidth+5, SubFrame2UIBottom, ThirdWidth-10, 20]);
set(QFTToolData.Elements(11),'pos',[SubFrame2UILeft+ThirdWidth*2+5-10, SubFrame2UIBottom, 10, 20]);

set(QFTToolData.Elements(9),'pos',[SubFrame2UILeft+ThirdWidth*2+5*2, SubFrame2UIBottom, ThirdWidth-10, 20]);
set(QFTToolData.Elements(12),'pos',[SubFrame2UILeft+ThirdWidth*3+5*2-10, SubFrame2UIBottom, 10, 20]);

SubFrame2UIBottom = SubFrame2UIBottom - 25;
set(QFTToolData.Elements(13),'pos',[SubFrame2UILeft, SubFrame2UIBottom, FullWidth, 20]);

%SubFrame2UIBottom = SubFrame2UIBottom - 25;
set(QFTToolData.Elements(16),'pos',[SubFrame2UILeft, SubFrame2UIBottom, FullWidth, 20]);

SubFrame2UIBottom = SubFrame2UIBottom - 25;
set(QFTToolData.Elements(14),'pos',[SubFrame2UILeft, SubFrame2UIBottom, HalfWidth, 20]);

set(QFTToolData.Elements(15),'pos',[SubFrame2UILeft+HalfWidth+5, SubFrame2UIBottom, HalfWidth, 20]);

SubFrame2UIBottom = SubFrame2UIBottom - 20;
set(QFTToolData.Elements(23),'pos',[SubFrame2UILeft, SubFrame2UIBottom, FullWidth, 15]);

ListboxBottom = SubFrame2Position(2) + 30;
ListboxHeight = (SubFrame2UIBottom) - ListboxBottom;
set(QFTToolData.Elements(17),'pos',[SubFrame2UILeft, ListboxBottom, FullWidth, ListboxHeight]);

SubFrame2UIBottom = SubFrame2Position(2) + 5;
set(QFTToolData.Elements(18),'pos',[SubFrame2UILeft, SubFrame2UIBottom, HalfWidth, 20]);
set(QFTToolData.Elements(19),'pos',[SubFrame2UILeft+HalfWidth+5, SubFrame2UIBottom, HalfWidth, 20]);
set(QFTToolData.Elements(20),'pos',[SubFrame2UILeft+ThirdWidth*2+5*2, SubFrame2UIBottom, ThirdWidth, 20]);

%% Bound (subframe 3)
set(QFTToolData.Bounds(1),'pos',SubFrame3Position);
set(QFTToolData.Bounds(2),'pos',[SubFrame3UILeft, SubFrame3UIBottom, 100, 17]);

ListboxBottom = SubFrame3Position(2) + 30;
ListboxHeight = (SubFrame3UIBottom) - ListboxBottom;
set(QFTToolData.Bounds(3),'pos',[SubFrame3UILeft, ListboxBottom, FullWidth, ListboxHeight]);

SubFrame3UIBottom = SubFrame3Position(2) + 5;
set(QFTToolData.Bounds(4),'pos',[SubFrame3UILeft, SubFrame3UIBottom, ThirdWidth, 20]);
set(QFTToolData.Bounds(5),'pos',[SubFrame3UILeft+ThirdWidth+5, SubFrame3UIBottom, ThirdWidth, 20]);
set(QFTToolData.Bounds(6),'pos',[SubFrame3UILeft+ThirdWidth*2+5*2, SubFrame3UIBottom, ThirdWidth, 20]);

set(QFTToolData.HintBar,'pos',[0,0,FigureWidth-MainFrameWidth-10,20]);

