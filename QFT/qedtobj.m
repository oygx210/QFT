function qedtobj(flag)
% QEDTOBJ Edit frequency weight line object. (Utility Function)
%         QEDTOBJ edits line objects within the Frequency Weight IDE.  It
%         handles all the special cases concerning the altering of line
%         properties.

% Author: Craig Borghesani
% 5/12/93
% Copyright (c) 2003, Terasoft, Inc.


f2=gcf;
f=get(f2,'userdata');
bthan=get(f,'userdata');
a=gca;
han = get(bthan(26),'userdata');
lims = get(han(11),'userdata');

red = [1,0,0];
blue = [0,0,1];
cyan = [0,1,1];

cur_obj = get(f2,'currentobject');
objects = get(han(3),'userdata');
obj_len = length(objects);

if flag==0,
 set(f2,'windowbuttonupfcn','qedtobj(1)');
end

x_pts = eps;
if obj_len,
 ct=1;
 for k=1:2:obj_len,
  x_pts(ct) = get(objects(k),'xdata');
  ct=ct+1;
 end
end

pt1 = get(a,'currentpoint');

if all(pt1(1,1:2)>lims([1,3])) & all(pt1(1,1:2)<lims([2,4])),
 xlt = find(pt1(1,1) < x_pts);
 if length(xlt),
  objx_loc = 2*xlt(1)-1;
  if objx_loc > 1,
   if strcmp(get(objects(objx_loc-1),'vis'),'on'),
    pt0 = get(objects(objx_loc-2),'userdata');
    set(objects(objx_loc-1),'xdata',[pt0(1),pt1(1,1)],'ydata',[pt0(2),pt1(1,2)]);
    clr = get(objects(objx_loc-1),'userdata');
    pt0 = get(objects(objx_loc),'userdata');
    objects = [objects,0,0];
    objects((objx_loc+2):(obj_len+2)) = objects(objx_loc:obj_len);
    objects(objx_loc) = line('xdata',pt1(1,1),'ydata',pt1(1,2),...
                             'marker','o','color','g','markersize',4,...
                             'userdata',pt1(1,1:2),'erase','xor');
    objects(objx_loc+1) = line('xdata',[pt1(1,1),pt0(1)],'ydata',[pt1(1,2),pt0(2)],...
                               'linestyle','-','color',clr,...
                               'userdata',clr,'erase','xor',...
                               'buttondownfcn','qbtnpres');
    set(han(3),'userdata',objects);
   else
    errordlg('No line to EDIT','Message','on');
   end
  else
   errordlg('No line to EDIT','Message','on');
  end
 else
  errordlg('No line to EDIT','Message','on');
 end
end
