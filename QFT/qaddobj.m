function qaddobj(flag)
% QADDOBJ Add frequency weight line object. (Utility Function)
%         QADDOBJ adds line objects within the Frequency Weight IDE.  It
%         handles all the special cases concerning placement of the lines
%         depending upon point selection from the user.

% Author: Craig Borghesani
% 5/12/93
% Copyright (c) 2003, Terasoft, Inc.


f2=gcf;
f=get(f2,'userdata');
a=gca;
bthan = get(f,'userdata');
han = get(bthan(26),'userdata');
lims = get(han(11),'userdata');

red = [1,0,0];
blue = [0,0,1];
cyan = [0,1,1];

cur_obj = get(f,'currentobject');
objects = get(han(3),'userdata');
obj_len = length(objects);

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
 if flag==0, % begin adding a new line
  set(han(7),'userdata',1);
  set(han([3,12]),'enable','off');
  set(han([5,11]),'enable','on');
  if all(x_pts < pt1(1,1)), % adding new line away from any other lines
   if obj_len,
    obj_len = obj_len + 1;
    pt0 = get(objects(obj_len-1),'userdata');
    objects(obj_len) = line('xdata',[pt0(1,1),pt1(1,1)],'ydata',[pt0(1,2),pt1(1,2)],...
                            'linestyle','-','color','r','erase','xor',...
                            'userdata',red,...
                            'buttondownfcn','qbtnpres','vis','off');
    objects(obj_len+1) = line('xdata',pt1(1,1),'ydata',pt1(1,2),...
                              'marker','o','color','g','markersize',4,...
                              'userdata',pt1(1,1:2),'erase','xor');
   else
    obj_len = obj_len + 1;
    objects(obj_len) = line('xdata',pt1(1,1),'ydata',pt1(1,2),...
                            'marker','o','color','g','markersize',4,...
                            'userdata',pt1(1,1:2),'erase','xor');
   end
   set(a,'userdata',[]);

% adding a new line and the beginning starts over/under another
  elseif any(x_pts < pt1(1,1)),
   xlt = find(pt1(1,1) < x_pts);
   objx_loc = 2*xlt(1)-1;
   if objx_loc > 1 & objx_loc < obj_len,  % over/under a present line
    if strcmp(get(objects(objx_loc-1),'vis'),'off'),
     pt0 = get(objects(objx_loc-2),'userdata');
     set(objects(objx_loc-1),'xdata',[pt0(1),pt1(1,1)],'ydata',[pt0(2),pt1(1,2)]);
     pt0 = get(objects(objx_loc),'userdata');
     objects = [objects,0,0];
     objects((objx_loc+2):(obj_len+2)) = objects(objx_loc:obj_len);
     objects(objx_loc) = line('xdata',pt1(1,1),'ydata',pt1(1,2),...
                             'marker','o','color','g','markersize',4,...
                             'userdata',pt1(1,1:2),'erase','xor');
     objects(objx_loc+1) = line('xdata',[pt1(1,1),pt0(1)],'ydata',[pt1(1,2),pt0(2)],...
                             'linestyle','-','color','r',...
                             'userdata',red,'erase','xor',...
                             'buttondownfcn','qbtnpres','vis','off');
     set(a,'userdata',objx_loc);
    else
     if strcmp(get(objects(objx_loc+1),'vis'),'off'),
      set(objects(objx_loc),'vis','off');
     end
     xdata = get(objects(objx_loc-1),'xdata');
     ydata = get(objects(objx_loc-1),'ydata');
     set(objects(objx_loc-1),'xdata',[xdata(1),pt1(1,1)]);
     pt0 = get(objects(objx_loc),'userdata');
     objects = [objects,0,0,0,0];
     objects((objx_loc+4):(obj_len+4))=objects(objx_loc:obj_len);
     objects(objx_loc) = line('xdata',pt1(1,1),'ydata',ydata(2),...
                             'marker','o','color','g','markersize',4,...
                             'userdata',[pt1(1,1),ydata(2)],'erase','xor');
     objects(objx_loc+1) = line('xdata',[pt1(1,1),pt1(1,1)],'ydata',[ydata(2),pt1(1,2)],...
                             'linestyle','-','color',red,...
                             'userdata',red,'erase','xor',...
                             'buttondownfcn','qbtnpres','vis','off');
     objects(objx_loc+2) = line('xdata',pt1(1,1),'ydata',pt1(1,2),...
                             'marker','o','color','g','markersize',4,...
                             'userdata',pt1(1,1:2),'erase','xor');
     objects(objx_loc+3) = line('xdata',[pt1(1,1),pt0(1,1)],'ydata',[pt1(1,2),pt0(1,2)],...
                             'linestyle','-','color',red,...
                             'userdata',red,'erase','xor',...
                             'buttondownfcn','qbtnpres','vis','off');
     set(a,'userdata',objx_loc+2);
    end
   else
    xdata = get(objects(objx_loc-1),'xdata');
    ydata = get(objects(objx_loc-1),'ydata');
    set(objects(objx_loc-1),'xdata',[xdata(1),pt1(1,1)]);
    set(objects(objx_loc),'xdata',pt1(1,1),'ydata',ydata(2));
    set(objects(objx_loc),'userdata',[pt1(1,1),ydata(2)]);
    obj1 = line('xdata',[pt1(1,1),pt1(1,1)],'ydata',[ydata(1),pt1(1,2)],...
                'linestyle','-','color',red,...
                'userdata',red,'erase','xor',...
                'buttondownfcn','qbtnpres','vis','off');
    obj2 = line('xdata',pt1(1,1),'ydata',pt1(1,2),...
                'marker','o','color','g','markersize',4,...
                'userdata',pt1(1,1:2),'erase','xor');
    objects = [objects,obj1,obj2];
    set(a,'userdata',[]);
   end
  else % adding a line and the beginning starts before all lines
   pt0 = get(objects(1),'userdata');
   obj1 = line('xdata',pt1(1,1),'ydata',pt1(1,2),...
               'marker','o','color','g','markersize',4,...
               'userdata',pt1(1,1:2),'erase','xor');
   obj2 = line('xdata',[pt1(1,1),pt0(1,1)],'ydata',[pt1(1,2),pt0(1,2)],...
               'linestyle','-','color','r',...
               'userdata',red,'erase','xor',...
               'buttondownfcn','qbtnpres','vis','off');
   objects = [obj1,obj2,objects];
   set(a,'userdata',1);
  end
  set(f2,'windowbuttonupfcn','qaddobj(1)');
  set(han(3),'userdata',objects);
 elseif flag==1,
  set(han(7),'userdata',0);
  set(han([3,12]),'enable','on');
  lastx = get(a,'userdata');
  if any(x_pts > pt1(1,1)) & length(lastx),
   xlt = find(pt1(1,1) < x_pts);
   obj_loc = 2*xlt(1)-1;
   if (obj_loc > lastx) & (obj_loc <= obj_len),
    if (obj_loc-lastx)>2,
     if obj_loc < obj_len,
      set(objects((lastx+2):(obj_loc-3)),'vis','off');
      if strcmp(get(objects(obj_loc+1),'vis'),'off'),
       set(objects(obj_loc),'vis','off');
      end
      pt0 = get(objects(obj_loc),'userdata');
      set(objects(obj_loc-1),'xdata',[pt1(1,1),pt0(1)],'ydata',[pt1(1,2),pt0(2)],...
                             'vis','off');
      set(objects(obj_loc-2),'xdata',pt1(1,1),'ydata',pt1(1,2));
      if strcmp(get(objects(obj_loc-2),'vis'),'off'),
       set(objects(obj_loc-2),'vis','on');
      end
      set(objects(obj_loc-2),'userdata',pt1(1,1:2));
      pt0 = get(objects(lastx),'userdata');
      set(objects(lastx+1),'xdata',[pt0(1),pt1(1,1)],...
                             'ydata',[pt0(2),pt1(1,2)],'vis','on');
      set(objects(lastx+1),'color',red);
      set(objects(lastx+1),'userdata',red);
      objects((lastx+2):(obj_loc-3))=[];
      set(a,'userdata',lastx+2);
     else
      set(objects((lastx+3):obj_loc),'vis','off');
      pt0 = get(objects(lastx),'userdata');
      set(objects(lastx+1),'xdata',[pt0(1),pt1(1,1)],...
                             'ydata',[pt0(2),pt1(1,2)],'vis','on');
      set(objects(lastx+1),'color',red);
      set(objects(lastx+1),'userdata',red);
      if strcmp(get(objects(lastx+2),'vis'),'off'),
       set(objects(lastx+2),'xdata',pt1(1,1),'ydata',pt1(1,2),'vis','on');
      else
       set(objects(lastx+2),'xdata',pt1(1,1),'ydata',pt1(1,2));
      end
      set(objects(lastx+2),'userdata',pt1(1,1:2));
      objects((lastx+3):obj_loc)=[];
      set(a,'userdata',[]);
     end
    else
     objects = [objects,0,0];
     objects((obj_loc+2):(obj_len+2))=objects(obj_loc:obj_len);
     pt0 = get(objects(lastx),'userdata');
     set(objects(lastx+1),'xdata',[pt0(1),pt1(1,1)],'ydata',[pt0(2),pt1(1,2)],...
                          'vis','on');
     set(objects(lastx+1),'color',red);
     set(objects(lastx+1),'userdata',red);
     pt0 = get(objects(obj_loc),'userdata');
     objects(obj_loc) = line('xdata',pt1(1,1),'ydata',pt1(1,2),...
                             'marker','o','color','g','markersize',4,...
                             'userdata',[pt1(1,1:2)],'erase','xor');
     objects(obj_loc+1) = line('xdata',[pt1(1,1),pt0(1)],'ydata',[pt1(1,2),pt0(2)],...
                             'linestyle','-','color','r',...
                             'userdata',red,'erase','xor',...
                             'buttondownfcn','qbtnpres','vis','off');
     set(a,'userdata',lastx+2);
    end
   else  % at the last segment of a line
    pt0 = get(objects(lastx),'userdata');
    set(objects(obj_loc-1:obj_loc),'vis','off');
    set(objects(obj_loc-2),'xdata',pt1(1,1),'ydata',pt1(1,2));
    set(objects(obj_loc-2),'userdata',pt1(1,1:2));
    set(objects(lastx+1),'xdata',[pt0(1),pt1(1,1)],'ydata',[pt0(2),pt1(1,2)],...
                         'vis','on');
    set(objects(lastx+1),'color',red);
    set(objects(lastx+1),'userdata',red);
    objects((obj_loc-1):obj_loc)=[];
    set(a,'userdata',[]);
   end
  elseif all(x_pts < pt1(1,1)),  % add to the end of a line
   if length(lastx),
    pt0 = get(objects(lastx),'userdata');
    set(objects((lastx+2):(obj_len-1)),'vis','off');
    set(objects(obj_len),'xdata',pt1(1,1),'ydata',pt1(1,2));
    set(objects(obj_len),'userdata',pt1(1,1:2));
    set(objects(lastx+1),'xdata',[pt0(1),pt1(1,1)],...
                         'ydata',[pt0(2),pt1(1,2)],'vis','on');
    objects((lastx+2):(obj_len-1))=[];
    set(a,'userdata',[]);
   else
    obj_len = obj_len + 1;
    pt0 = get(objects(obj_len-1),'userdata');
    objects(obj_len) = line('xdata',[pt0(1),pt1(1,1)],'ydata',[pt0(2),pt1(1,2)],...
                            'linestyle','-','color','r',...
                            'userdata',red,'erase','xor',...
                            'buttondownfcn','qbtnpres');
    objects(obj_len+1) = line('xdata',pt1(1,1),'ydata',pt1(1,2),...
                            'marker','o','color','g','markersize',4,...
                            'userdata',pt1(1,1:2),'erase','xor');
   end
  else
   errordlg('Invalid operation','Message','on');
  end
  set(han(3),'userdata',objects);
 end
end
