function copybnds(f)
% COPYBNDS Copy bounds. (Utility Function)
%          COPYBNDS copies the bounds from -360deg to 0deg to other ranges if
%          the user specifies new axis limits.

% Author: Craig Borghesani
% Date: 10/10/93
% Revised: 2/17/96 2:53 PM V1.1 update
% Copyright (c) 2003, Terasoft, Inc.


bthan=get(f,'userdata');
infmat=get(bthan(16),'userdata');
bdb=get(bthan(32),'userdata');
phase=get(bthan(33),'userdata');
bhan=get(bthan(12),'userdata');  % handle to all bound buttons
axs=infmat(1,:);
lastphmin=infmat(28,2); lastphmax=infmat(28,3);
phmin=axs(1); phmax=axs(2);
%figure(f);
set(0,'currentfigure',f);

QFTToolData = getappdata(f,'QFTToolData');
BoundsListbox = QFTToolData.Bounds(3);
BoundsListboxInfo = get(BoundsListbox,'userdata');

if length(bdb),
 [coora,coorb]=wherebnd(bdb);
 if phmin<lastphmin,
  newbd=ceil(abs(phmin)/360)-1;
  if newbd==0, newbd=1; end
%  taxs=axs;
%tphase=phase;
  for k=newbd:-1:1,
   tphase=phase-360*k;
   taxs(1:2)=axs(1:2)-360*k;
   newbhan=qplotbd(tphase,bdb,coora,coorb,taxs);
   ow=newbhan(6,:); ow=sort(ow); ow(find(diff(ow)==0))=[];
   for h=1:length(ow),
    bnddata = BoundsListboxInfo{h};
%    bnddata=get(bhan(h),'userdata');
    bndvis = get(bnddata(1),'vis');
    newdata=newbhan(1:5,find(ow(h)==newbhan(6,:)));
    bnddata=[bnddata;newdata];
    BoundsListboxInfo{h} = bnddata;
%    set(bhan(h),'userdata',bnddata);
    set(bnddata,'vis',bndvis);
   end

% this is the On/Off button
%   bnddata=get(bhan(length(ow)+1),'userdata');
%   set(bhan(length(ow)+1),'userdata',[bnddata;newbhan(1:5,:)]);

  end
 end
 if phmax>lastphmax,
  newbd=ceil(abs(phmax)/360);
  if newbd==0, newbd=1; end
%  taxs=axs;
%tphase=phase;
  for k=newbd:-1:1,
   tphase=phase+360*k;
   taxs(1:2)=axs(1:2)+360*k;
   newbhan=qplotbd(tphase,bdb,coora,coorb,taxs);
   ow=newbhan(6,:); ow=sort(ow); ow(find(diff(ow)==0))=[];
   for h=1:length(ow),
    bnddata = BoundsListboxInfo{h};
%    bnddata=get(bhan(h),'userdata');
    bndvis = get(bnddata(1),'vis');
    newdata=newbhan(1:5,find(ow(h)==newbhan(6,:)));
    bnddata=[bnddata;newdata];
    BoundsListboxInfo{h} = bnddata;
%    set(bhan(h),'userdata',bnddata);
    set(bnddata,'vis',bndvis);
   end

% this is the on/off button
%   bnddata=get(bhan(length(ow)+1),'userdata');
%   set(bhan(length(ow)+1),'userdata',[bnddata;newbhan(1:5,:)]);
  end
 end
end

if phmin<lastphmin,
 lastphmin=floor(phmin/360)*360;
 infmat(28,2)=lastphmin;
end
if phmax>lastphmax,
 lastphmax=ceil(phmax/360)*360;
 infmat(28,3)=lastphmax;
end
set(bthan(16),'userdata',infmat);
set(BoundsListbox, 'userdata', BoundsListboxInfo);

