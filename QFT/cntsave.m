function cntsave(flag)
% CNTSAVE Controller save. (Utility Function)
%         CNTSAVE allows the user to save the controller matrix to a file.

% Author: Craig Borghesani
% 9/6/93
% Copyright (c) 2003, Terasoft, Inc.


if flag==3,
 f2 = gcf;
 f = get(f2,'userdata');
else
 f=gcf;
end

bthan=get(f,'userdata');
infmat=get(bthan(16),'userdata');
T=get(bthan(13),'userdata');
lomat = get(bthan(1),'userdata');
hint_bar = get(bthan(36),'userdata');
last_name = get(bthan(15),'userdata');
proc_num = int2str(infmat(25,2));

wl = lomat(1,:);
cont_r=get(bthan(3),'userdata');
cont2 = get(bthan(31),'userdata');
if ~isempty(cont2{2}),
 cont_r(1,1) = cont_r(1,1)*cont2{2}(1,1);
 if T > 0,
  cont_r(3,1) = cont_r(3,1)+cont2{2}(3,1);
  cont2{2}(1:3,:) = [];
 else
  cont2{2}(1:2,:) = [];
 end
 cont_r = [cont_r;cont2{2}];
end

flag2=infmat(9,1)*2-(T==0);

if any(flag==[1,3]),

 proc_str=[];
 if infmat(25,2)>1, proc_str=['(',int2str(infmat(25,2)),')']; end

 str=str2mat('*.shp','*.dsh','*.fsh','*.dfs','*.bod','*.dbo');
 if ~length(last_name),
  file_name = str(flag2,:);
 else
  file_name = last_name;
 end

 str2=['Save Elements File ',proc_str];

 [fname,pth]=uiputfile(file_name,str2);

 if isstr(fname),
% add the extension if it is not there
  if isempty(findstr('.', fname)),
   fname = [fname, str(flag2, 2:5)];
  end
  set(bthan(15),'userdata',fname);
  eval(['save ''',pth,fname,''' cont_r wl T -mat']);
  set(hint_bar,'string',['Saving <',pth,fname,'>...']);
 end

 if flag==3 & fname,
  close(findobj('tag',['qft1',proc_num]));
  close(findobj('tag',['qft2',proc_num]));
  close(findobj('tag',['qft3',proc_num]));
  close(findobj('tag',['qft4',proc_num]));
  close(findobj('tag',['qft5',proc_num]));
  close(findobj('tag',['qft6',proc_num]));
  close(findobj('tag',['qft7',proc_num]));
  close(f);
  leftover = get(0,'chil');
  for k=1:length(leftover),
   if strcmp(get(leftover(k),'name'),'Message'),
    close(leftover(k));
   end
  end
 end

elseif flag==2,

 if flag2==1, fname='shape.shp';
 elseif flag2==2, fname='dshape.dsh';
 elseif flag2==3, fname='fshape.fsh';
 elseif flag2==4, fname='dfshape.dfs';
 elseif flag2==5, fname='bodplot.bod';
 elseif flag2==6, fname='dbodplot.dbo';
 end

 eval(['save ',fname,' cont_r wl T;']);
 set(hint_bar,'string',['Saving <',fname,'>...']);

end
