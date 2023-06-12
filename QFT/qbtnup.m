function qbtnup
% QBTNUP Button up handler. (Utility Function)
%        QBTNUP handles the button up function of the mouse and uses the
%        this instance to store x-y data regarding line objects within the
%        frequency weight IDE.

% Author: Craig Borghesani
% 5/12/93
% Copyright (c) 2003, Terasoft, Inc.


f2=gcf;
f=get(f2,'userdata');
bthan=get(f,'userdata');
a=gca;
han = get(bthan(26),'userdata');

cur_obj = get(f2,'currentobject');
objects = get(han(3),'userdata');
sel_mat = get(han(6),'userdata');
obj_len = length(objects);
sz=size(sel_mat);

sloc = 2:3:sz(1);
sloc2 = sel_mat(sloc,1)';
for obj_loc = sloc2,
 mxdata = get(objects(obj_loc),'xdata');
 mydata = get(objects(obj_loc),'ydata');
 set(objects(obj_loc-1),'userdata',[mxdata(1),mydata(1)]);
 set(objects(obj_loc+1),'userdata',[mxdata(2),mydata(2)]);
 set(objects(obj_loc),'buttondownfcn','qbtnpres');
 set(objects(obj_loc),'color',get(objects(obj_loc),'userdata'));
end
set(f2,'windowbuttonupfcn','','windowbuttonmotionfcn','','pointer','arrow');
set(han(4:5),'enable','on');
set(han(6:10),'enable','off');
set(han(6),'userdata',[]);
