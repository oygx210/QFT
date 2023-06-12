function qscrelmt(flag)
% QSCRELMT Setup mouse window function. (Utility Function)
%          QSCRELMT initializes the specific mouse button down functions
%          for the adding of elements from the screen: K, 1, 2, L/L, NTC

% Author: Craig Borghesani
% 9/6/93
% Copyright (c) 2003, Terasoft, Inc.

f=gcf;
bthan=get(f,'userdata');
infmat=get(bthan(16),'userdata');
hint_bar = get(bthan(36),'userdata');

set(infmat(8,1),'enable','on');
%set(bthan(27),'userdata',lomat);
%set(bthan(28),'userdata',cont);

T=get(bthan(13),'userdata');

if infmat(9,1)==3,
   set(infmat(6,1:2),'vis','off','userdata',[]);
else
   set(infmat(6,1),'vis','off','userdata',[]);
end

set(bthan(31),'userdata',{[], []});

% reseting button-state functions
set(f,'windowbuttonupfcn','');
if length(get(f,'windowbuttonmotionfcn'))>1,
   set(f,'windowbuttonmotionfcn','qmouse(''Floating'')');
end

% reset complex and frequency values
infmat(4,3)=0;
infmat(3,3)=0;
set(bthan(16),'userdata',infmat);

if any(infmat(9,1)==[1,2,3]),

 opr=['mogain(0,0)  ';
      'mofirst(0,0) ';
      'mosecond(0,0)';
      'moldlg(0,0)  ';
      'montch(0,0)  ';
      'mo2ovr2(0)   ';
      'mocpld(0,0)'];

 set(f,'windowbuttondownfcn',opr(flag,:))
 set(hint_bar,'string','Select a point on the frequency response by pressing and holding down the mouse button');

end
