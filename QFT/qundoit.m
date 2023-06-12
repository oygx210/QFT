function qundoit(flag)
% QUNDOIT Un-do. (Utility Function)
%         QUNDOIT undoes the last frequency response change.

% Author: Craig Borghesani
% 10/10/93
% Copyright (c) 2003, Terasoft, Inc.
%       $Revision: 1.5 $

T = [];

if nargin==0,
 f=gcf;
 bthan=get(f,'userdata');
 infmat=get(bthan(16),'userdata');
 set(infmat(8,1),'enable','off');

 old_lomat=get(bthan(27),'userdata');
 old_cont=get(bthan(28),'userdata');

 qclswin(0);
 if length(old_lomat),
  set(bthan(1),'userdata',old_lomat);
  if infmat(9,1)==1, qnicplt(f);
  elseif infmat(9,1)==2, qmagplt(f);
  elseif infmat(9,1)==3, mgphplot(f);
  end
 end
 set(bthan(3),'userdata',old_cont);

elseif any(flag==[1,2]), % UnDo for Add and Edit
 f2=gcf;
 f=get(f2,'userdata');
 bthan=get(f,'userdata');
 infmat=get(bthan(16),'userdata');
 T = get(bthan(13),'userdata');
 cont = get(bthan(19),'userdata');
 lomat = get(bthan(20),'userdata');
 w = lomat(1,:);
 selctn=get(bthan(30),'userdata');

 q=1; r=2;
 if infmat(9,1)==2, q=[1;1]; r=2:3; end
 if flag==1, % Undo for Add
  edt1 = 16; txt1 = 19;
  loc = length(cont(:,1));
  if selctn(1,1),
   loc = selctn(1,1);
   val = selctn(1,2);
   ele = cont(loc,4);
   if any(ele==[0.6,0.7]),
    old_ele=[val,NaN,NaN,ele];
    cont(loc,1)=-val+cont(loc,1);
   elseif ele == 0.5, % discrete integrators
    old_ele=[0,0,NaN,ele];
    old_ele(1,1+(val<0))=old_ele(1,1+(val<0))+abs(val);
    cont(2,1+(val<0))=cont(2,1+(val<0))-abs(val);
    if (cont(2,1)-cont(2,2))==0, cont(2,1:2)=[0 0]; end
   end
  else
   old_ele=cont(loc,:);
   ele = old_ele(4);
   cont(loc,:)=[];
  end
  old_cp = qcntbode(old_ele,w,T);
  set(infmat(13,4),'callback','qelmts(7,0)');

 else % Undo for Edit
  edt1 = 17; txt1 = 20;
  loc = selctn(1,1);
  butn_sel = selctn(1,2);
  cur_str = get(butn_sel,'string');
  brac = find(cur_str=='[');
  new_str = cur_str(1:brac);
  cur_ele = cont(loc,:);
  old_ele = [selctn(1,3:5),cur_ele(4)];
  ele = cur_ele(4);
  old_cp = qcntbode(cur_ele,w,T)./qcntbode(old_ele,w,T);
  cont(loc,:) = old_ele;
 end

 lomat(r,:) = lomat(r,:)./old_cp(q,:);

 selctn(1,:)=[];
 set(bthan(19),'userdata',cont);
 set(bthan(20),'userdata',lomat);
 set(bthan(30),'userdata',selctn);
 set(infmat(txt1,1:3),'string','');
 set(infmat(edt1,1:3),'string','','enable','off');
 drawnow;

 if ele==0,
  set(infmat(edt1,1),'string',num2str(old_ele(1)),'enable','on');
  set(infmat(txt1,1),'string','gain=');
  if flag==2,
   v1=num2str(old_ele(1),4);
   set(butn_sel,'string',[new_str,v1,']']);
  end

 elseif any(ele==[1 2]),
  str=['pole=';'zero='];
  set(infmat(edt1,1),'string',num2str(old_ele(1)),'enable','on');
  set(infmat(txt1,1),'string',str(ele,:));
  if flag==2,
   v1=num2str(old_ele(1),4);
   set(butn_sel,'string',[new_str,v1,']']);
  end

 elseif any(ele==[3 4]),
  set(infmat(edt1,2),'string',num2str(old_ele(1)),'enable','on');
  set(infmat(edt1,1),'string',num2str(old_ele(2)),'enable','on');
  if ele==3,
   set(infmat(txt1,2),'string','zeta(pole)=');
   set(infmat(txt1,1),'string','wn(pole)=');
  else
   set(infmat(txt1,2),'string','zeta(zero)=');
   set(infmat(txt1,1),'string','wn(zero)=');
  end
  if flag==2,
   v1=num2str(old_ele(1),4);
   v2=num2str(old_ele(2),4);
   set(butn_sel,'string',[new_str,v1,', ',v2,']']);
  end

 elseif any(ele==[0.6,0.7]),  % add cont. integrator/differentiator
  set(infmat(edt1,1),'string',num2str(old_ele(1)),'enable','on');
  set(infmat(txt1,1),'string','n=');

 elseif ele==0.5,  % add discrete integrator/differentiator
  set(infmat(edt1,1),'string',num2str(old_ele(1,1+(val<0))),'enable','on');
  set(infmat(txt1,1),'string','n=');

 elseif ele==5,
  set(infmat(edt1,2),'string',num2str(old_ele(1)),'enable','on');
  set(infmat(edt1,1),'string',num2str(old_ele(2)),'enable','on');
  set(infmat(txt1,2),'string','phase=');
  set(infmat(txt1,1),'string','w=');
  if flag==2,
   v1=num2str(old_ele(1),4);
   v2=num2str(old_ele(2),4);
   set(butn_sel,'string',[new_str,v1,', ',v2,']']);
  end

 elseif ele==6,
  set(infmat(edt1,3),'string',num2str(old_ele(1)),'enable','on');
  set(infmat(edt1,2),'string',num2str(old_ele(2)),'enable','on');
  set(infmat(edt1,1),'string',num2str(old_ele(3)),'enable','on');
  set(infmat(txt1,3),'string','zeta1(zero)=');
  set(infmat(txt1,2),'string','zeta2(pole)=');
  set(infmat(txt1,1),'string','wn=');
  if flag==2,
   v1=num2str(old_ele(1),4);
   v2=num2str(old_ele(2),4);
   v3=num2str(old_ele(3),4);
   set(butn_sel,'string',[new_str,v1,', ',v2,', ',v3,']']);
  end
 end

 if ~length(selctn),
  v2=get(bthan(21),'userdata');
  vo2=get(bthan(22),'userdata');
  if flag==1,
   set(infmat(13,3),'enable','off');
  else
   set(infmat(21,3),'enable','off');
  end
  if infmat(9,1)==1,
   set([v2(:);vo2(:)],'vis','off');
  else
   set(v2(:),'vis','off');
  end
 else
  if infmat(9,1)==1, qnicplt(f);
  elseif infmat(9,1)==2, qmagplt(f);
  elseif infmat(9,1)==3, mgphplot(f);
  end
 end
 if flag == 1,
  set(infmat(13,1),'callback',['qelmts(',num2str(ele),',0)']);
  set(infmat(13,4),'callback','qelmts(7,0)');
  cntdisp(f,cont,0);
 else
  other_rads = get(butn_sel,'userdata');
  set(butn_sel,'value',1); set(other_rads,'value',0);
  set(infmat(21,1),'userdata',butn_sel);
  set(infmat(21,[1,4]),'callback',['qelmts(',num2str(ele),',',int2str(loc),');']);
 end

elseif flag == 3, % Undo for Iterate

 f2=gcf;
 f=get(f2,'userdata');
 bthan=get(f,'userdata');
 infmat=get(bthan(16),'userdata');
 cont_old = get(bthan(3),'userdata');
 cont_new = get(bthan(19),'userdata');
 lomat_new = get(bthan(20),'userdata');
 w = lomat_new(1,:);
 data = get(infmat(21,3),'userdata');
 flag2 = data(1);
 butn_sel = data(2);
 cur_str = get(butn_sel,'string');
 brac = find(cur_str=='[');
 new_str = cur_str(1:brac);
 q=1; r=2;
 if infmat(9,1)==2, q=[1;1]; r=2:3; end

 old_cp = qcntbode(cont_new(flag2,:),w,T)./qcntbode(cont_old(flag2,:),w,T);
 lomat_new(r,:)=lomat_new(r,:)./old_cp(q,:);
 cont_new(flag2,:)=cont_old(flag2,:);

 if cont_new(flag2,4)==0,
  set(infmat(17,1),'value',cont_new(flag2,1));
  if cont_new(flag2,1)>0,
    set(infmat(17,1),'min',10^((20*log10(cont_new(flag2,1))-5)/20),...
                     'max',10^((20*log10(cont_new(flag2,1))+5)/20));
  else
    set(infmat(17,1),'max',sign(cont_new(flag2,1))*10^((20*log10(abs(cont_new(flag2,1)))-5)/20),...
            'min',sign(cont_new(flag2,1))*10^((20*log10(abs(cont_new(flag2,1)))+5)/20));
  end

  v1=num2str(cont_new(flag2,1),4);
  set(butn_sel,'string',[new_str,v1,']']);

 elseif any(cont_new(flag2,4)==[1 2]),
  set(infmat(17,1),'min',real(cont_new(flag2,1))*0.5,'max',real(cont_new(flag2,1))*1.5,...
        'value',real(cont_new(flag2,1)));

  v1=num2str(cont_new(flag2,1),4);
  set(butn_sel,'string',[new_str,v1,']']);

 elseif any(cont_new(flag2,4)==[3 4]),
  val=cont_new(flag2,1); val2=cont_new(flag2,2);
  set(infmat(17,2),'min',val/2,'max',val*1.5,'value',val);
  set(infmat(17,1),'min',val2/2,'max',val2*1.5,'value',val2);

  v1=num2str(val,4);
  v2=num2str(val2,4);
  set(butn_sel,'string',[new_str,v1,', ',v2,']']);

 elseif cont_new(flag2,4)==5,
  set(infmat(17,2),'min',-87,'max',87,'value',cont_new(flag2,1));
  set(infmat(17,1),'min',cont_new(flag2,2)/2,'max',cont_new(flag2,2)*1.5,...
         'value',cont_new(flag2,2));

  v1=num2str(cont_new(flag2,1),4);
  v2=num2str(cont_new(flag2,2),4);
  set(butn_sel,'string',[new_str,v1,', ',v2,']']);

 elseif cont_new(flag2,4)==6,
  val=cont_new(flag2,1); val2=cont_new(flag2,2); val3=cont_new(flag2,3);
  set(infmat(17,3),'min',val/2,'max',val*1.5,'value',val);
  set(infmat(17,2),'min',val2/2,'max',val2*1.5,'value',val2);
  set(infmat(17,1),'min',val3/2,'max',val3*1.5,'value',val3);

  v1=num2str(val,4);
  v2=num2str(val2,4);
  v3=num2str(val3,4);
  set(butn_sel,'string',[new_str,v1,', ',v2,', ',v3,']']);

 end
 set(infmat(21,3),'enable','off');

 set(bthan(20),'userdata',lomat_new);
 set(bthan(19),'userdata',cont_new);
 if infmat(9,1)==1, qnicplt(f);
 elseif infmat(9,1)==2, qmagplt(f);
 elseif infmat(9,1)==3, mgphplot(f);
 end

end
