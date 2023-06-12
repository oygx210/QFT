function b=qplotbd2(bph,bdb,coora,coorb,axs,CurrentAxes)
% QPLOTBD Bound plotting driver. (Utility Function)
%         QPLOTBD is the driver to PLOTBNDS and plots the bounds in the
%         figure set up by PLOTBNDS.

% Author: Craig Borghesani
% 9/6/93
% Copyright (c) 2003, Terasoft, Inc.
%       $Revision: 1.3 $

[r,c]=size(bdb);
cvec=['r';'g';'b';'c';'m'];
clr=[cvec;cvec;cvec;cvec];
clr=[clr;clr;clr;clr;clr];
clr=[clr;clr;clr;clr;clr];
w=bdb(r-1,:); w=sort(w); w(find(diff(w)==0))=[];
myeps=1e-16;
ers = 'norm';
CurrentAxesNextPlot = get(CurrentAxes,'nextplot');
set(CurrentAxes,'nextplot','add');

dbmyeps=20*log10(myeps); db1myeps=20*log10(1/myeps);
offset=(r-2)/2;
for k=1:length(coora(:,1)),
 abdb=0; abph=0; bldb=0; blph=0; ldb=0; lph=0; rdb=0; rph=0; t=0;
 astr='-w';bstr=astr;sstr=astr;
 z=coora(k,1); flag1=1; flag2=1; once=1;
 if bdb(r,z)==13, flag=0; else flag=1; end
 az=coora(k,2):coora(k,3); bz=coorb(k,2):coorb(k,3);
 laz=length(az); lbz=length(bz);
 if az(1)~=0,
  astr=['-',clr(find(bdb(r-1,z)==w))];
  abdb=bdb(az,z); abph=bph(az);
  if once & flag,
   nstr=int2str(bdb(r,z));
%   t=text('pos',[abph(1)-8,abdb(1)],'string',nstr,...
%          'horizontalalignment','right');
%   set(t,'clipping','on','erase',ers);
   once=0;
  end
 end
 if bz(1)~=0,
  bstr=['--',clr(find(bdb(r-1,z)==w))];
  bldb=bdb(bz+offset,z); blph=bph(bz);
  if once & flag,
   nstr=int2str(bdb(r,z));
%   t=text('pos',[blph(1)-8,bldb(1)],'string',nstr,...
%          'horizontalalignment','right');
%   set(t,'clipping','on','erase',ers);
   once=0;
  end
 end
 if (laz~=offset & az(1)~=0) | (lbz~=offset & bz(1)~=0),
  sstr=['-',clr(find(bdb(r-1,z)==w))];
  if az(1)~=0,
   if ~any(bph(az(1))==axs(1:2)) & bdb(az(1)+offset,z)~=db1myeps,
    rph=[bph(az(1)); bph(az(1))];
    rdb=[bdb(az(1),z); bdb(az(1)+offset,z)]; flag1=0;
   end
   if ~any(bph(az(laz))==axs(1:2)) & bdb(az(laz)+offset,z)~=db1myeps,
    lph=[bph(az(laz)); bph(az(laz))];
    ldb=[bdb(az(laz),z); bdb(az(laz)+offset,z)]; flag2=0;
   end
  end
  if bz(1)~=0,
   if flag1 & (~any(bph(bz(1))==axs(1:2))) & bdb(bz(1),z)~=dbmyeps,
    rph=[bph(bz(1)); bph(bz(1))];
    rdb=[bdb(bz(1),z); bdb(bz(1)+offset,z)];
   end
   if flag2 & (~any(bph(bz(lbz))==axs(1:2))) & bdb(bz(lbz),z)~=dbmyeps,
    lph=[bph(bz(lbz)); bph(bz(lbz))];
    ldb=[bdb(bz(lbz),z); bdb(bz(lbz)+offset,z)];
   end
  end
 end
 b(1:4,k)=plot(abph,abdb,astr,blph,bldb,bstr,lph,ldb,sstr,rph,rdb,sstr,'parent',CurrentAxes);
 set(b(1:4,k),'erase',ers);
 b(1:6,k)=[b(1:4,k);t;w(find(bdb(r-1,z)==w))];
end
set(CurrentAxes,'nextplot',CurrentAxesNextPlot);
