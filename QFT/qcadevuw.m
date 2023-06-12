function qcadevuw(flag)
% QCADEVUW QView pulldown menu in IDEs. (Utility Function)
%          QCADEVUW manages the QView functions that implement Full, In, Out,
%          Last, Move, Zoom, Axis, and Elements.

% Author: Craig Borghesani
% 9/11/93
% Copyright (c) 2003, Terasoft, Inc.


f=gcf;
bthan=get(f,'userdata');
infmat=get(bthan(16),'userdata');
hint_bar = get(bthan(36),'userdata');

if flag==1,
 set(hint_bar,'string','Setting axis to FULL limits');
elseif flag==2,
 set(hint_bar,'string','Zooming IN by 20% of vertical axis');
elseif flag==3,
 set(hint_bar,'string','Zooming OUT by 20% of vertical axis');
elseif flag==4,
 set(hint_bar,'string','Setting axis to LAST (previous) limits');
elseif flag==11,
 set(hint_bar,'string','Select two points to define movement of axis limits');
elseif flag==13,
 set(hint_bar,'string','Press and drag mouse to define new axis limits window');
end
drawnow;

if infmat(9,1)==1,
 qnicvuw(flag);
elseif infmat(9,1)==2,
 qmagvuw(flag);
elseif infmat(9,1)==3,
 bodevuw(flag);
end
