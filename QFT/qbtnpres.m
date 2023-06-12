function qbtnpres
% QBTNPRES Button press handler. (Utility Function)
%          QBTNPRES handles the presses of the mouse button and selection
%          of line objects within the frequency weight IDE.

% Author: Craig Borghesani
% 5/12/93
% Copyright (c) 2003, Terasoft, Inc.



f2=gcf;
f=get(f2,'userdata');
bthan=get(f,'userdata');
a=gca;
han = get(bthan(26),'userdata');
on_top = get(han(7),'userdata');
if ~on_top,
 objects = get(han(3),'userdata');
 sel_mat = get(han(6),'userdata');

 obj_len = length(objects);
 cur_obj = get(f2,'currentobject');
 obj_loc = find(cur_obj==objects);
 obj_vec = 1:length(objects);
 set(han([4:10]),'enable','off');

 set(f2,'windowbuttonupfcn','','windowbuttondownfcn','',...
        'windowbuttonmotionfcn','','pointer','arrow');

 sloc=[];
 if length(sel_mat),
  sloc = find(obj_loc==sel_mat(:,1));
 end

 if ~length(sloc),
  set(han(10),'enable','on');
  set(objects(obj_loc),'color','y');
  if (obj_loc-1) > 1,
   xdata = get(objects(obj_loc-2),'xdata');
   ydata = get(objects(obj_loc-2),'ydata');
   sel_mat = [sel_mat;0,objects(obj_loc-2),xdata,ydata];
  else
   sel_mat = [sel_mat;0,0,0,0,0,0];
  end

  xdata = get(objects(obj_loc),'xdata');
  ydata = get(objects(obj_loc),'ydata');
  sel_mat = [sel_mat;obj_loc,objects(obj_loc),xdata,ydata];

  if (obj_loc+1) < obj_len,
   xdata = get(objects(obj_loc+2),'xdata');
   ydata = get(objects(obj_loc+2),'ydata');
   sel_mat = [sel_mat;0,objects(obj_loc+2),xdata,ydata];
  else
   sel_mat = [sel_mat;0,0,0,0,0,0];
  end
 elseif length(sel_mat),
  set(han(10),'enable','on');
  set(objects(obj_loc),'color',get(objects(obj_loc),'userdata'));
  sel_mat((sloc-1):(sloc+1),:)=[];
  if ~length(sel_mat),
   set(han([4,5]),'enable','on');
   set(han(10),'enable','off');
  end
 end
 sz = size(sel_mat);

 if sz(1) == 3,
  set(han([6,7]),'enable','on');
  obj_locs = sel_mat(2,1);
  if ((obj_locs-1) > 1) & ((obj_locs+1) < obj_len),
   if strcmp(get(objects(obj_locs-2),'vis'),'on') & ...
      strcmp(get(objects(obj_locs+2),'vis'),'on'),
    set(han(8),'enable','on');
   end
  end
 elseif sz(1)~=0,
  set(han(7),'enable','on');
  obj_locs = sel_mat(2:3:sz(1),1);
  min_obj = min(obj_locs);
  max_obj = max(obj_locs);
  if all(diff(sort(obj_locs))==2),
   set(han(6),'enable','on');
   if all(obj_locs ~= 2) & all(obj_locs ~= (obj_len-1)),
    if strcmp(get(objects(min_obj-2),'vis'),'on') & ...
         strcmp(get(objects(max_obj+2),'vis'),'on'),
     set(han(8),'enable','on');
    end
   end
  elseif sz(1)==6,
   if diff(sort(obj_locs)) > 2,
    if diff(sort(obj_locs)) > 4,
     set(han(9),'enable','on');
    elseif strcmp(get(objects(min_obj+2),'vis'),'off'),
     set(han(9),'enable','on');
    end
   end
  end
 end
 set(han(6),'userdata',sel_mat);
else
 errordlg('Cannot be on top of line object','Message','on');
end
