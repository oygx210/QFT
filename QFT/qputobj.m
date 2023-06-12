function qputobj(flag)
% QPUTOBJ Read frequency weight files. (Utility Function)
%         QPUTOBJ reads and plots frequency weight files from within the
%         freuquency weight IDE.

% Author: Craig Borghesani
% 10/10/93
% Copyright (c) 2003, Terasoft, Inc.


a=gca;
f2=gcf;
f=get(f2,'userdata');
bthan=get(f,'userdata');

infmat = get(bthan(16),'userdata');
han = get(bthan(26),'userdata');
objects = get(han(3),'userdata');
lims = get(han(11),'userdata');

proc_str=[];
if infmat(25,2)>1, proc_str=['(',int2str(infmat(25,2)),')']; end

if flag==0,
 [fname,pth]=uigetfile('*.fwf',['Open Frequency Weight Function ',proc_str]);
 if isstr(fname),
  delete(objects);
  eval(['load ',pth,fname,' -mat']);
  obj_len = length(obj_data(:,1));
  xmin = min(obj_data(:,1));
  xmax = max(obj_data(:,3));
  set(a,'xlim',[min(lims(1),xmin),max(lims(2),xmax)]);
  lims(1:2) = [min(lims(1),xmin),max(lims(2),xmax)];
  ct=1;
  for k = 1:obj_len,
   str='on';
   if ~obj_data(k,5), str = 'off'; end
   pt0(1) = obj_data(k,1);
   pt0(2) = 10^((log10(lims(4))-log10(lims(3)))*obj_data(k,2)+log10(lims(3)));
   if rem(k,2),
    objects(k) = line('xdata',pt0(1),'ydata',pt0(2),...
                      'marker','o','color','g','markersize',4,...
                      'userdata',pt0,'erase','xor');
   else
    pt1(1) = obj_data(k,3);
    pt1(2) = 10^((log10(lims(4))-log10(lims(3)))*obj_data(k,4)+log10(lims(3)));
    clr = obj_data(k,6:8);
    objects(k) = line('xdata',[pt0(1),pt1(1)],'ydata',[pt0(2),pt1(2)],...
                      'linestyle','-','color',clr,'erase','xor',...
                      'userdata',clr,...
                      'buttondownfcn','qbtnpres',...
                      'vis',str);
   end
  end
  set(han(3),'userdata',objects);
  set(han(11),'userdata',lims);
 end
elseif flag==1,
 [fname,pth]=uiputfile('*.fwf',['Save Frequency Weight Function ',proc_str]);
 if isstr(fname),
  obj_len = length(objects);
  for k = 1:obj_len,
   vis = 1;
   if strcmp(get(objects(k),'vis'),'off'), vis = 0; end
   if rem(k,2),
    data = get(objects(k),'userdata');
    x = data(1);
    yper = (log10(data(2))-log10(lims(3)))/(log10(lims(4))-log10(lims(3)));
    obj_data(k,:) = [x,yper,0,0,vis,0,0,0];
   else
    data = get(objects(k+1),'userdata');
    clr = get(objects(k),'userdata');
    x2 = data(1);
    y2per = (log10(data(2))-log10(lims(3)))/(log10(lims(4))-log10(lims(3)));
    obj_data(k,:) = [x,yper,x2,y2per,vis,clr];
   end
  end
  eval(['save ',pth,fname,' obj_data;']);
 end
end
