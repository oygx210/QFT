function [fphs]=qfixfase(lo,axs,wloc)
% QFIXFASE Extends phase range outside [-360,0] so that it
%          becomes continuous across Nichols chart phase axis and
%          guarantees that it wraps correctly depending upon
%          axis limits

% Author: Yossi Chait, Craig Borghesani
% 5/9/94
% Copyright (c) 2003, Terasoft, Inc.


phs=qatan4(lo)*180/pi;
ip=length(phs);
fphs=phs;
if any(fphs~=0),

 d=diff(fphs);
 id=find(abs(d) > 180)+1; % add 1 to correct for location after diff
 for j=1:length(id),
  fphs(id(j):ip)=fphs(id(j):ip)-sign(diff(phs([id(j)-1,id(j)])))*360;
 end

%ax1=floor(axs(1)/360)*360; ax2=ceil(axs(2)/360)*360;
%range=ax2-ax1; shift=ceil(range/360)*360;
 ax1=axs(1); ax2=axs(2);
 if ax2<0,
  shift1=ceil(abs(ax2)/360)*360;
 else
  shift1=floor(abs(ax2)/360)*360;
 end
 fphs=fphs+(sign(ax2)*shift1);
 shift2=floor((ax2-ax1)/360)*360;
 if shift2==0, shift2=360; end
 brkl=find(fphs<ax1); brkr=find(fphs>ax2);
 while length(brkl), fphs(brkl)=fphs(brkl)+shift2; brkl=find(fphs<ax1); end
 while length(brkr), fphs(brkr)=fphs(brkr)-shift2; brkr=find(fphs>ax2); end

 if nargin==3,
  fphs=fphs(wloc);
 end
elseif nargin==3,
 fphs=fphs(wloc);
end
