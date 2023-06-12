function qzoomaxs(flag)
% QZOOMAXS Zoom tool for frequency weight IDE. (Utility Function)
%          QZOOMAXS is the zoom tool for changing the axis settings of the
%          frequency weight IDE.
%

% Author: Craig Borghesani
% 5/12/93
% Copyright (c) 2003, Terasoft, Inc.
%       $Revision: 1.5 $

a=gca;
f2=gcf;
f=get(f2,'userdata');
bthan=get(f,'userdata');
han = get(bthan(26),'userdata');
lims = get(han(11),'userdata');

if flag==0, % Initialize Zoom
 set(han(7),'userdata',get(f2,'windowbuttondownfcn'));
 set(f2,'pointer','cross');
 set(f2,'windowbuttondownfcn','qzoomaxs(1)','interruptible','On',...
        'windowbuttonupfcn','1;','units','pixels');
elseif flag==1,
 pt1=get(a,'currentpoint');
 rbbox([get(f2,'currentpoint'),0,0],get(f2,'currentpoint'));
 drawnow;
 pt2=get(a,'currentpoint');
 x = [pt1(1,1),pt2(1,1)];
 xlimp = [min(x),max(x)];
 y = [pt1(1,2),pt2(1,2)];
 ylimp = [min(y),max(y)];
 set(a,'xlim',xlimp,'ylim',ylimp);
 old_button_down = get(han(7),'userdata');
 set(f2,'pointer','arrow','windowbuttondownfcn',old_button_down,...
        'windowbuttonupfcn','','interruptible','Off','units','norm');
elseif flag==3,
 set(a,'xlim',lims(1:2),'ylim',lims(3:4));
end
