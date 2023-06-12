function qdelobj(flag)
% QDELOBJ Delete frequency weight line object. (Utility Function)
%         QDELOBJ deletes line objects within the Frequency Weight IDE.  It
%         handles all the special cases concerning removal of the lines
%         selected by the user.

% Author: Craig Borghesani
% 5/12/93
% Copyright (c) 2003, Terasoft, Inc.


f2=gcf;
f=get(f2,'userdata');
bthan=get(f,'userdata');
a=gca;
han = get(bthan(26),'userdata');
objects = get(han(3),'userdata');
obj_len = length(objects);
sel_mat = get(han(6),'userdata');
sz = size(sel_mat);
tvec=[]; tl=0;

red = [1,0,0];
blue = [0,0,1];
cyan = [0,1,1];

if flag==0, % Delete
 obj_locs = sort(sel_mat(2:3:sz(1),1)');
 for obj_loc=obj_locs,
  if length(tvec),
   tl = tl + length(tvec);
   obj_loc = obj_loc-tl;
   obj_len = obj_len-length(tvec);
  end
  if ((obj_loc+1) == obj_len),
   tvec = obj_loc:(obj_loc+1);
   if obj_len > 3,
    if strcmp(get(objects(obj_loc-2),'vis'),'off'),
     tvec = (obj_loc-2):(obj_loc+1);
    end
    set(objects(tvec),'vis','off');
    objects(tvec)=[];
   else
    set(objects,'vis','off');
    objects=[];
   end
  elseif ((obj_loc-1) == 1),
   tvec = (obj_loc-1):obj_loc;
   if strcmp(get(objects(obj_loc+2),'vis'),'off'),
    tvec = (obj_loc-1):(obj_loc+2);
   end
   set(objects(tvec),'vis','off');
   objects(tvec)=[];
  else
   vec1 = [-1,3,2];
   tvec = obj_loc:(obj_loc+1);
   if strcmp(get(objects(obj_loc-2),'vis'),'off') & strcmp(get(objects(obj_loc+2),'vis'),'off'),
    vec1 = [-3,3,2];
    tvec = (obj_loc-2):(obj_loc+1);
   elseif strcmp(get(objects(obj_loc-2),'vis'),'off'),
    vec1 = [-3,1,-2];
    tvec = (obj_loc-1):obj_loc;
   end
   set(objects(tvec),'vis','off');
   xdata1 = get(objects(obj_loc+vec1(1)),'xdata');
   ydata1 = get(objects(obj_loc+vec1(1)),'ydata');
   xdata2 = get(objects(obj_loc+vec1(2)),'xdata');
   ydata2 = get(objects(obj_loc+vec1(2)),'ydata');
   set(objects(obj_loc+vec1(3)),'xdata',[xdata1,xdata2],'ydata',[ydata1,ydata2]);
   objects(tvec)=[];
  end
 end
elseif flag==1, % Break
 sel_mat(sel_mat(:,1)==0,1) = sel_mat(sel_mat(:,1)==0,1)+1000;
 [jk,min_obj] = min(sel_mat(:,1));
 sel_mat(sel_mat(:,1)==1000,1) = sel_mat(sel_mat(:,1)==1000,1)-1000;
 [jk,max_obj] = max(sel_mat(:,1));
 obj_vec = sel_mat(min_obj,1):sel_mat(max_obj,1);
 set(objects(obj_vec),'vis','off');
 if length(obj_vec) > 1,
  set(objects(obj_vec(1)),'xdata',[sel_mat(min_obj,3),sel_mat(max_obj,4)],...
                          'ydata',[sel_mat(min_obj,5),sel_mat(max_obj,6)],...
                          'vis','off');
  objects((sel_mat(min_obj,1)+1):sel_mat(max_obj,1))=[];
 end
elseif flag==2, % Connect
 sel_mat(sel_mat(:,1)==0,1) = sel_mat(sel_mat(:,1)==0,1)+1000;
 [jk,min_obj] = min(sel_mat(:,1));
 sel_mat(sel_mat(:,1)==1000,1) = sel_mat(sel_mat(:,1)==1000,1)-1000;
 [jk,max_obj] = max(sel_mat(:,1));
 obj_vec = sel_mat(min_obj,1):sel_mat(max_obj,1);
 clr = red;
 if all(get(sel_mat(min_obj,2),'userdata')==get(sel_mat(max_obj,2),'userdata')),
  clr = get(sel_mat(min_obj,2),'userdata');
 end
 if length(obj_vec) > 5,
  set(objects(obj_vec(4):(sel_mat(max_obj,1)-2)),'vis','off');
  set(objects(obj_vec(3)),'xdata',[sel_mat(min_obj,4),sel_mat(max_obj,3)],...
                          'ydata',[sel_mat(min_obj,6),sel_mat(max_obj,5)]);
  if strcmp(get(objects(obj_vec(3)),'vis'),'off'),
   set(objects(obj_vec(3)),'vis','on');
   set(objects(obj_vec(3)),'userdata',clr);
  end
  objects(obj_vec(4):(sel_mat(max_obj,1)-2))=[];
 else
  set(objects(obj_vec(3)),'vis','on');
  set(objects(obj_vec(3)),'userdata',clr);
 end
elseif flag==3, % Clear
 delete(objects);
 objects = [];
elseif flag==4, % Select-None
 if length(sel_mat),
  obj_locs = sel_mat(2:3:sz(1),1);
  set(objects(obj_locs),'color','r');
 end
elseif flag==5, % Select-All
 sel_mat=[];
 for obj_loc = 2:2:obj_len,
  if strcmp(get(objects(obj_loc),'vis'),'on'),
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
  end
 end
 set(han(6),'userdata',sel_mat);
elseif flag==6, % Input
 obj_locs = sel_mat(2:3:sz(1),1);
 set(objects(obj_locs),'userdata',blue);
elseif flag==7, % Output
 obj_locs = sel_mat(2:3:sz(1),1);
 set(objects(obj_locs),'userdata',cyan);
elseif flag==8, % Both
 obj_locs = sel_mat(2:3:sz(1),1);
 set(objects(obj_locs),'userdata',red);
end
if flag~=5,
 set(han(6),'userdata',[]);
 set(han(3),'userdata',objects);
 set(han(6:10),'enable','off');
 set(han(4:5),'enable','on');
 obj_len = length(objects);
 for k = 2:2:obj_len,
  set(objects(k),'color',get(objects(k),'userdata'));
 end
 if ~obj_len,
  set(han([5,11]),'enable','off');
 end
else
 set(han([4:10]),'enable','off');
 set(han([7,10]),'enable','on');
 sz=size(sel_mat);
 obj_locs = sel_mat(2:3:sz(1),1);
 if all(diff(sort(obj_locs))==2) | sz(1)==3,
  set(han(6),'enable','on');
 end
end
set(f2,'pointer','arrow','windowbuttonupfcn','','windowbuttondownfcn','',...
       'windowbuttonmotionfcn','');
