function qwatecad(flag)
% QWATECAD Frequency weight IDE. (Utility Function)
%          QWATECAD creates the frequency weight IDE within LPSHAPE/DLPSHAPE
%          environments.

% Author: Craig Borghesani
% 10/10/93
% Copyright (c) 2003, Terasoft, Inc.

f2=gcf;
f=get(f2,'userdata');

bthan = get(f,'userdata');
infmat = get(bthan(16),'userdata');

h = get(bthan(25),'userdata');
h2 = get(bthan(26),'userdata');
wfunc = get(h2(3),'userdata');

proc_str=[];
if infmat(25,2)>1, proc_str=[' (',int2str(infmat(25,2)),')']; end

if flag==0,
   cont=get(bthan(19),'userdata');
   lomat=get(bthan(20),'userdata');
   if isempty(cont),
      cont = get(bthan(3),'userdata');
      lomat=get(bthan(1),'userdata');
      set(bthan(19),'userdata',cont);
      set(bthan(20),'userdata',lomat);
   end
   w = lomat(1,:);
   lo = abs(lomat(2,:));
   set(f2,'name',['Frequency Weighting',proc_str]);
   t = get(h2(12),'userdata');
   set(h([4,5,6,7,20]),'enable','off');
   set(h(8),'callback','qrelmts(''Use Weights'')');
   set(h(9),'callback','qwatecad(1)');
   set(h2([3,12]),'enable','on');
   xlims = [min(w),max(w)];
   ylims = [min(lo)/10,max(lo)*10];
   pos = get(h(3),'userdata');
   set(h(1),'xscale','log','yscale','log',...
            'xlim',xlims,'ylim',ylims,'xtickmode','auto','ytickmode','auto');
   set(h2(11),'userdata',[xlims,ylims]);
   set(h2(6),'userdata',[]);
   set(h(2),'xdata',w,'ydata',lo);
   set(h(3),'vis','off');
   set(wfunc(1:2:length(wfunc)),'color','g');
   for loc = 2:2:length(wfunc),
      set(wfunc(loc),'color',get(wfunc(loc),'userdata'));
   end
   set(h2(4:11),'enable','off');
   if length(wfunc),
      set(h2([4,5,11]),'enable','on');
   else
      set(h2(4),'enable','on');
   end
% if length(t),
%  ct=1;
%  for val = ylims(1):((ylims(2)-ylims(1))/4):ylims(2),
%   tval = (val - ylims(1))/(ylims(2) - ylims(1));
%   set(t(ct),'pos',[xlims(2)+xlims(2)*0.07,val],'string',sprintf('%1.4f',tval));
%   ct = ct + 1;
%  end
% else
%  ct=1;
%  for val = ylims(1):((ylims(2)-ylims(1))/4):ylims(2),
%   tval = (val - ylims(1))/(ylims(2) - ylims(1));
%   t(ct) = text('pos',[xlims(2)+xlims(2)*0.07,val],'string',sprintf('%1.4f',tval));
%   ct = ct + 1;
%  end
%  set(h2(12),'userdata',t);
% end
elseif flag==1,
   set(f2,'pointer','arrow','windowbuttondownfcn','',...
          'windowbuttonupfcn','','windowbuttonmotionfcn','',...
          'name',['Hankel Singular Values',proc_str]);
   set(h2([3,12]),'enable','off');
   set(wfunc,'color',[0,0,0]);
   hsv = get(h(2),'userdata');
   set(h(2),'ydata',hsv,'xdata',1:length(hsv));
   set(h(3),'ydata',hsv,'xdata',1:length(hsv),'vis','on');
   pos=get(h(3),'userdata');
   xticks = 0:length(hsv)+1;
   set(h(1),'pos',pos,'xtick',xticks,'xlim',[0,length(hsv)+1],...
            'ylim',[10 .^[floor(log10(min(hsv))),ceil(log10(max(hsv)))]],...
            'yscale','log','xscale','linear');
   set(h([4,5,6,7]),'enable','on');
   set(h(8),'callback','qrelmts(''Done'')');
   set(h(9),'callback','qrelmts(''Cancel'')');
end
