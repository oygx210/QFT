function qdelelmt(flag)
% QDELELMT Set up delete interface. (Utility Function)
%          QDELELMT sets up the interface to delete elements chosen by
%          the user.

% Author: Craig Borghesani
% 9/5/93
% Copyright (c) 2003, Terasoft, Inc.


if flag,
 f2=gcf;
 f=get(f2,'userdata');
else
 qclswin(0);
 f=gcf;
end
bthan=get(f,'userdata');
infmat=get(bthan(16),'userdata');

if flag==0,
 T=get(bthan(13),'userdata');
 cont=get(bthan(3),'userdata');

 fig_color=[0.5,0.5,0.5];
 proc_str=[];
 if infmat(25,2)>1, proc_str=[' (',int2str(infmat(25,2)),')']; end
 % check if there is anything to delete
 if T == 0 & length(cont(:,1))==2 & cont(2,1)==0,
   errordlg('Nothing to Delete','Message','on');
 elseif T > 0 & length(cont(:,1))==3,
  if sum([cont(2,1:2) cont(3,1)])==0,
   errordlg('Nothing to Delete','Message','on');
  else
   set(bthan([1,8,29]),'enable','off');
   cntdisp(f,cont,7);
  end
 else
% setup uicontrols of elements in controller matrix
  set(bthan([1,8,29]),'enable','off');
  cntdisp(f,cont,7);
 end
elseif flag==1, % Delete
 flag2=infmat(9,1);
 T=get(bthan(13),'userdata');
 lomat=get(bthan(1),'userdata');
 cont = get(bthan(3),'userdata');
 t_ele = get(bthan(23),'userdata');
 set(bthan(27),'userdata',lomat);
 set(bthan(28),'userdata',cont);
 cur_sel = gco;
 go_for_it=1;

 data = get(cur_sel,'userdata');
 delete = get(cur_sel,'value');
 cur_ele = data(1:4);

 if cur_ele(4)==0.7,
  val=str2num(get(data(5),'string'));
  if length(val) & abs(val)<=abs(cur_ele(1)),
   if cur_ele(1)>0,
    del_ele=[val,cur_ele(2:4)];
   else
    del_ele=[-val,cur_ele(2:4)];
   end
  else
   errordlg('Cannot DELETE more than you have.','Message','on');
   set(data(5),'string',int2str(abs(cur_ele(1))));
   set(cur_sel,'value',1);
   go_for_it=0;
  end
 elseif cur_ele(4)==0.6,
  val=str2num(get(data(5),'string'));
  if length(val) & abs(val)<=abs(cur_ele(1)),
   if cur_ele(1)>0,
    del_ele=[val,cur_ele(2:4)];
   else
    del_ele=[-val,cur_ele(2:4)];
   end
  else
   errordlg('Cannot DELETE more than you have.','Message','on');
   set(data(5),'string',int2str(abs(cur_ele(1))));
   set(cur_sel,'value',1);
   go_for_it=0;
  end
 elseif cur_ele(4)==0.5,
  val1=str2num(get(data(5),'string'));
  val2=str2num(get(data(6),'string'));
  if length([val1,val2])==2 & abs(val1)<=abs(cur_ele(1)) ...
                            & abs(val2)<=abs(cur_ele(2)),

   del_ele=[val1,val2,cur_ele(3:4)];
  else
   errordlg('Cannot DELETE more than you have.','Message','on');
   set(data(5),'string',int2str(cur_ele(1)));
   set(data(6),'string',int2str(cur_ele(2)));
   set(cur_sel,'value',1);
   go_for_it=0;
  end
 else
  del_ele = cur_ele;
 end
 if go_for_it,
  if any(infmat(9,1)==[1 3]), % SHAPE/DSHAPE/BODPLOT/DBODPLOT
   q=1; s=0;
  elseif infmat(9,1)==2, % FSHAPE/DFSHAPE
   q=[1;1]; s=1;
  end
  set(infmat(21,1),'enable','on');
  set(infmat(21,2),'callback','qdelelmt(3)');
  lomat2 = get(bthan(20),'userdata');
  if ~length(lomat2), lomat2 = lomat; end
  del_cp = qcntbode(del_ele,lomat2(1,:),T);
  if delete==0,
   lomat2(2:2+s,:) = lomat2(2:2+s,:)./del_cp(q,:);
   if length(data)==5, set(data(5),'enable','off');
   elseif length(data)==6, set(data(5:6),'enable','off'); end
  else
   if length(data)==6,
    set(data(5),'string',num2str(cur_ele(1)),'enable','on');
    set(data(6),'string',num2str(cur_ele(2)),'enable','on');
   elseif length(data)==5,
    set(data(5),'string',num2str(abs(cur_ele(1))),'enable','on');
   end
   lomat2(2:2+s,:) = lomat2(2:2+s,:).*del_cp(q,:);
  end

  done_flag = 1;
  for tct = 1:length(t_ele),
   if strcmp(get(t_ele(tct),'style'),'checkbox'),
    if ~get(t_ele(tct),'value') & strcmp(get(t_ele(tct),'vis'),'on'),
     done_flag = 0;
     break;
    end
   end
  end

  if done_flag,
   set(infmat(21,1),'enable','off');
   v2=get(bthan(21),'userdata');
   vo2=get(bthan(22),'userdata');
   if infmat(9,1)==1,
    set([v2(:);vo2(:)],'vis','off');
   else
    set(v2,'vis','off');
   end
   set(bthan(20),'userdata',[]);
  else
   set(bthan(20),'userdata',lomat2);
   if flag2==1, qnicplt(f);
   elseif flag2==2, qmagplt(f);
   elseif flag2==3, mgphplot(f);
   end
  end
 end
elseif flag == 2, % Done
 cont = get(bthan(3),'userdata');
 T = get(bthan(13),'userdata');
 if T > 0, cont_new = cont(1:3,:);
 else cont_new = cont(1:2,:); end
 objects = get(bthan(23),'userdata');
 for k = 1:length(objects),
  if strcmp(get(objects(k),'vis'),'on'),
   data = get(objects(k),'userdata');
   kept = get(objects(k),'value');
   if ~kept & length(data)==6, % discrete integrators/differentiators
    loc = find(data(4)==cont_new(:,4));
    cont_new(loc,1:2) = cont_new(loc,1:2) - ...
            [str2num(get(data(5),'string')),str2num(get(data(6),'string'))];
   elseif ~kept & length(data)==5, % integ/diff/pred/delay
    loc = find(data(4)==cont_new(:,4));
    if cont_new(loc,1) > 0, % integ/pred
     cont_new(loc,1) = cont_new(loc,1) - str2num(get(data(5),'string'));
    else % diff/delay
     cont_new(loc,1) = cont_new(loc,1) + str2num(get(data(5),'string'));
    end
   elseif kept & length(data)==4,
    cont_new = [cont_new;data];
   end
  end
 end
 lomat = get(bthan(20),'userdata');
 if length(lomat),
  set(bthan(1),'userdata',lomat);
 end
 set(bthan(3),'userdata',cont_new);
 set(bthan([1,8,29]),'enable','on');
 set(bthan([19,20]),'userdata',[]);
 v=get(bthan(10),'userdata');
 v2=get(bthan(21),'userdata');
 if infmat(9,1)==1,
  vo2=get(bthan(22),'userdata');
  vo=get(bthan(17),'userdata');
  set([v,vo],'vis','off');
  set(bthan(22),'userdata',vo);
  set(bthan(17),'userdata',vo2);
 else
  set(v,'vis','off');
 end
 set(v2,'linestyle','-');
 set(v,'linestyle','--');
 set(bthan(10),'userdata',v2);
 set(bthan(21),'userdata',v);
 set(infmat(8,1),'enable','on');
 set(f2,'vis','off');

elseif any(flag == [3,4]), % Cancel
 set(bthan([1,8,29]),'enable','on');
 if flag==3, % Cancel
  v2=get(bthan(21),'userdata');
  vo2=get(bthan(22),'userdata');
  if infmat(9,1)==1,
   set([v2(:);vo2(:)],'vis','off');
  else
   set(v2,'vis','off');
  end
 end
 set(bthan([19,20]),'userdata',[]);
 set(f2,'vis','off');
end
