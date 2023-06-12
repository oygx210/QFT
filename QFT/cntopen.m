function cntopen(flag)
% CNTOPEN Open controller file. (Utility Function)
%         CNTOPEN allows the user to open a saved controller matrix.

% Author: Craig Borghesani
% 9/6/93
% Copyright (c) 2003, Terasoft, Inc.


f=gcf;
bthan=get(f,'userdata');
infmat=get(bthan(16),'userdata');
QFTToolData = getappdata(f,'QFTToolData');

old_cont=get(bthan(3),'userdata');
old_lomat=get(bthan(1),'userdata');
hint_bar = get(bthan(36),'userdata');

set(bthan(27),'userdata',old_lomat);
set(bthan(28),'userdata',old_cont);
set(infmat(8,1),'enable','on');
lomat=get(bthan(1),'userdata');

T=get(bthan(13),'userdata');
last_name = get(bthan(15),'userdata');
flag2=infmat(9,1);

if flag==1,
 proc_str=[];
 if infmat(25,2)>1, proc_str=['(',int2str(infmat(25,2)),')']; end

 str=str2mat('*.shp;*.dsh','*.fsh;*.dfs');
 file_name = str(flag2,:);

 str2=str2mat(['Open Controller File ',proc_str],...
              ['Open Filter File ',proc_str],...
              ['Open Loop File ',proc_str]);

 [fname,pth]=uigetfile(file_name,str2(infmat(9,1),:));

else

 pth = '';
 if flag2==1, fname='shape.shp';
 elseif flag2==3, fname='fshape.fsh';
 end

 if ~exist(fname),
  errordlg(['Quick Open could not find <',fname,'>'],'Message','on');
  fname = 0;
 end

end

q=1;
if infmat(9,1)==2, q=[1;1]; end

if isstr(fname),
 set(bthan(15),'userdata',fname);
 qclswin(0);
 eval(['load ''',pth,fname,''' -mat']);
 wl=lomat(1,:); nom=get(bthan(2),'userdata');
 delay=infmat(10,1);
 if isempty(T), T = 0; end
 T2=T;
 T=get(bthan(13),'userdata');
 go_for_it=1;
 if T > 0 & T2 > 0,
  if T~=T2,
   go_for_it=2;
  end
 elseif T2 == 0 & T > 0,
  go_for_it=3;
 elseif T == 0 & T2 > 0,
  go_for_it=4;
 end

 if go_for_it==1,
    if T == 0,
      if isnan(cont_r(2,2)),
         if cont_r(2,1) >= 0,
            cont_r(2,2) = 0;
         elseif cont_r(2,1) < 0,
            cont_r(2,2) = abs(cont(2,1));
            cont_r(2,1) = 0;
         end
      end
    else
      if isnan(cont_r(3,2)),
         if cont_r(3,1) > 0,
            cont_r(3,2) = 0;
         elseif cont_r(3,1) < 0,
            cont_r(3,2) = abs(cont(3,1));
            cont_r(3,1) = 0;
         end
      end

    end
  cp=qcntbode(cont_r,wl,T);
  if length(nom) & infmat(9,1)==1,
   ncp=squeeze(freqresp(nom,wl)).';
  elseif infmat(9,1)==1,
   conL0=get(bthan(4),'userdata');
   L0=get(bthan(8),'userdata');
   cp0=qcntbode(conL0,wl,T); ncp=L0./cp0;
  elseif infmat(9,1)==2,
   ncp=get(bthan(17),'userdata');
  elseif infmat(9,1)==3,
   if length(nom),
    ncp=ones(1,length(wl));
   else
    ncp=get(bthan(8),'userdata');
   end
  end
  lomat(2:2+(infmat(9,1)==2),:)=cp(q,:).*ncp;
  set(bthan(3),'userdata',cont_r);
  set(bthan(1),'userdata',lomat);
  if infmat(9,1)==1, qnicplt(f);
  elseif infmat(9,1)==2, qmagplt(f);
  elseif infmat(9,1)==3, mgphplot(f);
  end
  set(hint_bar,'string',['Opening <',pth,fname,'>...']);

  [ControllerString,ListboxInfo] = cntstr(f,cont_r);
  set(QFTToolData.Elements(17), 'string', ControllerString,...
                                'userdata', ListboxInfo);

 elseif go_for_it==2,
  if infmat(9,1)==1,
   errordlg(['Sampling time mismatch. Controller (',num2str(T2),') Environment (',num2str(T),')'],'Message','on');
  elseif infmat(9,1)==2,
   errordlg(['Sampling time mismatch. Pre-filter (',num2str(T2),') Environment (',num2str(T),')'],'Message','on');
  end
 elseif go_for_it==3,
  if infmat(9,1)==1,
   errordlg('Continuous-time controller cannot be used in DLPSHAPE','Message','on');
  elseif infmat(9,1)==2,
   errordlg('Continuous-time pre-filter cannot be used in DPFSHAPE','Message','on');
  end
 elseif go_for_it==4,
  if infmat(9,1)==1,
   errordlg('Discrete-time controller cannot be used in LPSHAPE','Message','on');
  elseif infmat(9,1)==2,
   errordlg('Discrete-time pre-filter cannot be used in PFSHAPE','Message','on');
  end
 end
end
