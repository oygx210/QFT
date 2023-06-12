function qmoveobj
% QMOVEOBJ Moves frequency weight line objects. (Utility Function)
%          QMOVEOBJ moves line objects within the Frequency Weight IDE.  It
%          handles all the special cases concerning movement of the lines
%          the user has selected.

% Author: Craig Borghesani
% 5/12/93
% Copyright (c) 2003, Terasoft, Inc.


f2=gcf;
f=get(f2,'userdata');
bthan=get(f,'userdata');
a=gca;
han = get(bthan(26),'userdata');
lims = get(han(11),'userdata');

cur_obj = get(f2,'currentobject');
objects = get(han(3),'userdata');
sel_mat = get(han(6),'userdata');
sz = size(sel_mat);
obj_len = length(objects);

pt1 = get(a,'currentpoint');

if all(pt1(1,1:2)>lims([1,3])) & all(pt1(1,1:2)<lims([2,4])),
 pt0 = get(a,'userdata');

 sel_mat(sel_mat(:,1)==0,1) = sel_mat(sel_mat(:,1)==0,1)+1000;
 [jk,sloc_mn] = min(sel_mat(:,1));
 sel_mat(sel_mat(:,1)==1000,1) = sel_mat(sel_mat(:,1)==1000,1)-1000;
 [jk,sloc_mx] = max(sel_mat(:,1));

 if ~any(cur_obj==[f2,a]),
  lcir = sel_mat(sloc_mn,[3,5]);
  rcir = sel_mat(sloc_mx,[4,6]);
  delx = log10(pt1(1,1))-log10(pt0(1,1));
  dely = log10(pt1(1,2))-log10(pt0(1,2));
  if (sel_mat(sloc_mn,1)-2) >= 2,
   lcir2 = sel_mat(sloc_mn-1,[3,5]);
   if (log10(lcir(1))+delx) < log10(lcir2(1)),
    delx=log10(lcir2(1))-log10(lcir(1));
   end
  end
  if (sel_mat(sloc_mx,1)+2) <= obj_len,
   rcir2 = sel_mat(sloc_mx+1,[4,6]);
   if (log10(rcir(1))+delx) > log10(rcir2(1)),
    delx=log10(rcir2(1))-log10(rcir(1));
   end
  end
  sloc = 2:3:sz(1);
  sloc2 = sel_mat(2:3:sz(1),1)';
  ct = 1;
  for lloc = sloc2,
   xdata = log10(sel_mat(sloc(ct),3:4))+delx;
   ydata = log10(sel_mat(sloc(ct),5:6))+dely;
   set(objects(lloc-1),'xdata',10^xdata(1),'ydata',10^ydata(1));
   set(objects(lloc),'xdata',10.^xdata,'ydata',10.^ydata);
   set(objects(lloc+1),'xdata',10^xdata(2),'ydata',10^ydata(2));
   ct=ct+1;
  end
  if (sel_mat(sloc_mn,1)-2) >= 2,
   set(sel_mat(sloc_mn-1,2),'xdata',[lcir2(1),10^(log10(lcir(1))+delx)],...
                            'ydata',[lcir2(2),10^(log10(lcir(2))+dely)]);
  end
  if (sel_mat(sloc_mx,1)+2) <= obj_len,
   set(sel_mat(sloc_mx+1,2),'xdata',[10^(log10(rcir(1))+delx),rcir2(1)],...
                            'ydata',[10^(log10(rcir(2))+dely),rcir2(2)]);
  end
  set(f2,'windowbuttonupfcn','qbtnup');
 end
end
