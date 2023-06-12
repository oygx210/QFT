function kw = qfindfrq(x,y,freq_data,axs,flag);
% QFINDFRQ Find x-y location of pointer. (Utility Function)
%          QFINDFRQ finds the nearest frequency to the point that was
%          selected by the mouse on the plot.

% Author: Craig Borghesani
% 9/5/93
% Copyright (c) 2003, Terasoft, Inc.


kw=[];
xlim = get(axs,'xlim');
ylim = get(axs,'ylim');
if nargin==4, % SHAPE/FSHAPE/BODPLOT(magnitude plot)
 if strcmp(get(axs,'xscale'),'log'),
  x_pts = log(freq_data(1,:));
  xlim = log(xlim);
  x = log(x);
 else
  x_pts = qfixfase(freq_data(2,:),xlim);
 end
 y_pts = 20*log10(abs(freq_data(2,:)));
else % BODPLOT(phase plot)
 x_pts = log(freq_data(1,:));
 xlim = log(xlim);
 x = log(x);
 y_pts = qfixfase(freq_data(2,:),ylim);
end

% normalize everything
normx_pts = (x_pts - xlim(1))/(xlim(2)-xlim(1));
normy_pts = (y_pts - ylim(1))/(ylim(2)-ylim(1));
normx = (x - xlim(1))/(xlim(2)-xlim(1));
normy = (y - ylim(1))/(ylim(2)-ylim(1));

difx=normx_pts-normx;
dify=normy_pts-normy;
dist=sqrt((difx).^2+(dify).^2);
[mindist,k]=min(dist);
if mindist<0.02,
 kw=k;
end
