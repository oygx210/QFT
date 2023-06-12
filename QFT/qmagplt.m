function qmagplt(f)
% QMAGPLT Plot magnitude frequency response. (Utility Function)
%         QMAGPLT plots the frequency response within the axis limits for
%         the CAD environment PFSHAPE.

% Author: Craig Borghesani
% 9/6/93
% Copyright (c) 2003, Terasoft, Inc.


bthan=get(f,'userdata');
infmat=get(bthan(16),'userdata');

axs=infmat(1,:);
vax=infmat(23,1);
lomat=get(bthan(1),'userdata');
lomat2=get(bthan(20),'userdata');

if ~length(lomat2),
 v=get(bthan(10),'userdata');

% need to create two seperate frequency vectors
 w1=lomat(1,:);
 w2=w1;

 cl=get(bthan(2),'userdata');

 db1=20*log10(abs(lomat(2,:)));
 if ~vax,
  rmv_pt1=find(db1<axs(3) | db1>axs(4) | w1<axs(1) | w1>axs(2));
  db1(rmv_pt1)=db1(rmv_pt1)+NaN;
  w1(rmv_pt1)=w1(rmv_pt1)+NaN;
 end

% Each response has its own frequency vector because one of the
% responses might exceed the axis limits before the other.  If only
% one frequency vector is used, then the other plot may begin
% disappearing prior to its actually exceeding its limits because of
% the other plot exceeding its limits

 set(v(1),'xdata',w1,'ydata',db1);
 set(v(1),'vis','on');

 if length(v)==2,
  db2=20*log10(abs(lomat(3,:)));
  if ~vax,
   rmv_pt2=find(db2<axs(3) | db2>axs(4) | w1<axs(1) | w1>axs(2));
   db2(rmv_pt2)=db2(rmv_pt2)+NaN;
   w2(rmv_pt2)=w2(rmv_pt2)+NaN;
  end
  set(v(2),'xdata',w2,'ydata',db2);
  set(v(2),'vis','on');
 end

else
 v=get(bthan(21),'userdata');
% need to create two seperate frequency vectors
 w1=lomat2(1,:);
 w2=w1;

 cl=get(bthan(2),'userdata');

 db1=20*log10(abs(lomat2(2,:)));
 if ~vax,
  rmv_pt1=find(db1<axs(3) | db1>axs(4) | w1<axs(1) | w1>axs(2));
  db1(rmv_pt1)=db1(rmv_pt1)+NaN;
  w1(rmv_pt1)=w1(rmv_pt1)+NaN;
 end

% Each response has its own frequency vector because one of the
% responses might exceed the axis limits before the other.  If only
% one frequency vector is used, then the other plot may begin
% disappearing prior to its actually exceeding its limits because of
% the other plot exceeding its limits

 set(v(1),'xdata',w1,'ydata',db1);
 set(v(1),'vis','on');

 if length(v)==2,
  db2=20*log10(abs(lomat2(3,:)));
  if ~vax,
   rmv_pt2=find(db2<axs(3) | db2>axs(4) | w1<axs(1) | w1>axs(2));
   db2(rmv_pt2)=db2(rmv_pt2)+NaN;
   w2(rmv_pt2)=w2(rmv_pt2)+NaN;
  end
  set(v(2),'xdata',w2,'ydata',db2);
  set(v(2),'vis','on');
 end
end
