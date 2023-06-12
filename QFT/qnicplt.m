function qnicplt(f)
% QNICPLT Plot nichols response. (Utility Function)
%         QNICPLT plots the magnitude and phase response for the CAD
%         environments LPSHAPE/DLPSHAPE

% Author: Craig Borghesani
% 8/31/93
% Copyright (c) 2003, Terasoft, Inc.
%       $Revision: 1.3 $

bthan=get(f,'userdata');
infmat=get(bthan(16),'userdata');
wbs=get(bthan(11),'userdata');
axs=infmat(1,:);
vax = infmat(23,1);
elementMarker = infmat(6,1);
kw = infmat(3,3);

% compute stabilizing beta
if infmat(25,4) & strcmp(get(infmat(25,4),'vis'),'on'),
   edt = get(bthan(39),'userdata');
   b = str2num(get(edt(1),'string'));
   data = stabbeta(b, lomat(1,:), lomat(2,:));
   set(edt(2),'string',sprintf('%4.10f',data(1)));
   set(edt(3),'string',sprintf('%4.10f',data(2)));
end

for klo = 1:2,
   if klo == 1,
      lomat = get(bthan(1),'userdata');
      v=get(bthan(10),'userdata');
      vo=get(bthan(17),'userdata');

   else
      lomat=get(bthan(20),'userdata');
      v=get(bthan(21),'userdata');
      vo=get(bthan(22),'userdata');

   end

   if length(lomat),
      w=lomat(1,:); lo=lomat(2,:);
      db=20*log10(abs(lo));
      ph=qfixfase(lo,axs);
      brk=find(abs(diff(ph))>170);
      t=1;
      pht=[];
      dbt=[];
      wt=[];
      for k=brk,
         pht=[pht,ph(t:k),NaN];
         dbt=[dbt,db(t:k),NaN];
         wt=[wt,w(t:k),NaN];
         t=k+1;
      end
      pht=[pht,ph(t:length(ph))];
      dbt=[dbt,db(t:length(db))];
      wt=[wt,w(t:length(w))];

      if ~vax,
         rmv_pt=find(dbt>axs(4) | dbt<axs(3) | pht>axs(2) | pht<axs(1));
         pht(rmv_pt)=pht(rmv_pt)+NaN;
         dbt(rmv_pt)=dbt(rmv_pt)+NaN;
      end
      
      set(v, 'xdata', pht, 'ydata', dbt, 'vis', 'on');
      
      if klo == 2 & kw ~= 0,
         set(elementMarker, 'xdata', pht(kw), 'ydata', dbt(kw));
      end

      if length(wbs),
         for j=1:length(wbs),
            q=find(wt>=wbs(j) & ~isnan(wt));
            q2=find(wt<=wbs(j) & ~isnan(wt));
            if ~length(q) | ~length(q2),
               loc(j)=NaN;
            else
               loc(j)=q(1);
            end
         end

         for j=1:length(wbs),
            if ~isnan(loc(j)),
               set(vo(j),'xdata',pht(loc(j)),'ydata',dbt(loc(j)));
               set(vo(j),'vis','on');
            end
         end
      end
   end
end
