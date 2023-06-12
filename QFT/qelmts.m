function qelmts(ElementType, ElementLocation, ElementOperation)
% QELMTS Compute and store individual terms. (Utility Function)
%        QELMTS computes the various frequency responses of the element
%        that is selected from within the Add or Edit window.

% Author: Craig Borghesani
% 9/5/93
% Copyright (c) 2003, Terasoft, Inc.
%       $Revision: 1.5 $

f=gcf;
bthan=get(f,'userdata');
infmat=get(bthan(16),'userdata');
QFTToolData = getappdata(f,'QFTToolData');

butn = get(f,'currentobject');

ElementEnvironment=infmat(9,1);
cont=get(bthan(19),'userdata');
lomat=get(bthan(20),'userdata');
if isempty(cont),
	cont = get(bthan(3),'userdata');
	lomat=get(bthan(1),'userdata');
	set(bthan(19),'userdata',cont);
	set(bthan(20),'userdata',lomat);
end
T=get(bthan(13),'userdata');
selctn=get(bthan(30),'userdata');
ElementsListbox = QFTToolData.Elements(17);
elementMarker = infmat(6,1);
kw = infmat(3,3);

% depending upon whether the mode is add or edit determines the location
% of the data that was entered
loc=16; loc2=13;

% setting up for whether FSHAPE/DFSHAPE is being used
if any(ElementEnvironment==[1 3]), % SHAPE/DSHAPE/BODPLOT/DBODPLOT
 q=1; s=0;
elseif ElementEnvironment==2, % FSHAPE/DFSHAPE
 q=[1;1]; s=1;
end

go_for_it=0;

if ElementType == 0,

   oldval=cont(ElementLocation,1);

   switch ElementOperation,
      case('Edit'),
         val=str2num(get(infmat(loc,1),'string'));
         if val > 0,
            set(infmat(17,1),'min',10^((20*log10(val)-5)/20),...
                             'max',10^((20*log10(val)+5)/20),...
                             'value',val,...
                             'sliderstep',[0.01, 0.1]);
         else
            set(infmat(17,1),'min',-10^((20*log10(abs(val))+5)/20),...
                             'max',-10^((20*log10(abs(val))-5)/20),...
                             'value',val,...
                             'sliderstep',[0.01, 0.1]);
         end;%if

      case('Iterate'),
         val = get(infmat(17,1),'value');
         slider_min=get(infmat(17,1),'min');
         slider_max=get(infmat(17,1),'max');
         if val <= slider_min | val >= slider_max,
            if val > 0,
               set(infmat(17,1),'min',10^((20*log10(val)-5)/20),...
                                'max',10^((20*log10(val)+5)/20),...
                                'value',val,...
                                'sliderstep',[0.01, 0.1]);
            else
               set(infmat(17,1),'min',-10^((20*log10(abs(val))+5)/20),...
                                'max',-10^((20*log10(abs(val))-5)/20),...
                                'value',val,...
                                'sliderstep',[0.01, 0.1]);
            end;%if
         end
%         val = sign(cont(ElementLocation,1))*val;
         set(infmat(16,1),'string',num2str(val));

   end

   if length(val),
      lomat(2:2+s,:)=lomat(2:2+s,:).*(val/oldval);
      cont(1,1)=val;
      go_for_it=1;
   else
      errordlg('Gain value needs to be a number','Message','on');
   end

elseif any(ElementType==[1 2]), % first order

   str=['Pole';'Zero'];

   switch ElementOperation,
      case('Add'),
         val=str2num(get(infmat(loc,1),'string'));

      case('Edit'),
         oldval=cont(ElementLocation,1);
         val=str2num(get(infmat(loc,1),'string'));
         if val > 0,
            set(infmat(17,1),'min',real(val)*0.5,...
                              'max',real(val)*1.5,...
                              'value',val,...
                              'sliderstep',[0.01, 0.1]);
         else
            set(infmat(17,1),'min',real(val)*1.5,...
                              'max',real(val)*0.5,...
                              'value',val,...
                              'sliderstep',[0.01, 0.1]);

         end;%if

      case('Iterate'),
         oldval=cont(ElementLocation,1);
         val=get(infmat(17,1),'value');
         slider_min=get(infmat(17,1),'min');
         slider_max=get(infmat(17,1),'max');
         if real(val) <= slider_min | real(val) >= slider_max,
            if val > 0,
               set(infmat(17,1),'min',real(val)*0.5,...
                                'max',real(val)*1.5,...
                                'sliderstep',[0.01, 0.1]);
            else
               set(infmat(17,1),'min',real(val)*1.5,...
                                'max',real(val)*0.5,...
                                'sliderstep',[0.01, 0.1]);

            end;%if
         end
         set(infmat(16,1),'string',num2str(val));

   end

   if length(val) & val~=0,
      if imag(val)~=0,
         sign_imag = sign(imag(val));
         val = real(val) + sign_imag*pi/T*i;
      end
      rtnv=rlroot(val,lomat(1,:),[(ElementType==2)-(ElementType==1) T]);
      if strcmp(ElementOperation,'Add'),
         rt = ones(size(rtnv));
         cont=[cont;val,NaN,NaN,ElementType];
         ElementLocation = size(cont,1);
      else
         rt=rlroot(oldval,lomat(1,:),[(ElementType==2)-(ElementType==1) T]);
         cont(ElementLocation,1)=val;
      end
      lomat(2:2+s,:)=lomat(2:2+s,:).*(rtnv(q,:)./rt(q,:));
      go_for_it=1;

   elseif length(val),
      errordlg([str(ElementType,:),' cannot be = 0'],'Message','on');

   else
      errordlg('First Order value must be a number','Message','on');

   end

elseif any(ElementType==[3 4]), % second order

   delay = infmat(10,1);
   switch ElementOperation,
      case('Add'),
         val = str2num(get(infmat(loc,1),'string'));
         val2 = str2num(get(infmat(loc,2),'string'));

      case('Edit'),
         val = str2num(get(infmat(loc,1),'string'));
         val2 = str2num(get(infmat(loc,2),'string'));
         oldval = cont(ElementLocation, 1);
         oldval2 = cont(ElementLocation, 2);

         set(infmat(17,1),'min',val*0.5,...
                           'max',val*1.5,...
                           'value',val,...
                           'sliderstep',[0.001, 0.1]);
         set(infmat(17,2),'min',val2*0.5,...
                           'max',val2*1.5,...
                           'value',val2,...
                           'sliderstep',[0.01, 0.1]);

      case('Iterate'),
         oldval = cont(ElementLocation,1);
         oldval2 = cont(ElementLocation,2);
         val = get(infmat(17,1),'value');
         val2 = get(infmat(17,2),'value');

         slider_min=get(infmat(17,1),'min');
         slider_max=get(infmat(17,1),'max');
         if val <= slider_min | val >= slider_max,
            set(infmat(17,1),'min',val*0.5,...
                             'max',val*1.5,...
                             'sliderstep',[0.001, 0.1]);
         end

         slider_min=get(infmat(17,2),'min');
         slider_max=get(infmat(17,2),'max');
         if val2 <= slider_min | val2 >= slider_max,
            set(infmat(17,2),'min',val2*0.5,...
                             'max',val2*1.5,...
                             'sliderstep',[0.01, 0.1]);
         end
         set(infmat(16,1),'string',num2str(val));
         set(infmat(16,2),'string',num2str(val2));

   end

   if length([val,val2]) == 2,
      if val2~=0,
         if strcmp(ElementOperation, 'Add'),
            nom=get(bthan(2),'userdata');
            cont=[cont;val,val2,NaN,ElementType];
            ElementLocation = size(cont,1);
            if length(nom) & abs(val)<0.7 & ElementEnvironment~=2,
               wlt=qfrqenh(val2,val,lomat(1,:),T);
               if ElementEnvironment==1,
                  nlo=squeeze(freqresp(nom,wlt)).';
                  clo=qcntbode(cont,wlt,T).*exp(-i*wlt*delay);
                  lo=nlo.*clo;
               else
                  lo=qcntbode(cont,wlt,T).*exp(-i*wlt*delay);
               end
               lomat=[wlt;lo;ones(1,length(lo))];

            else
               rt=cproot(val,val2,lomat(1,:),[(ElementType==4)-(ElementType==3) T]);
               lomat(2:2+s,:)=lomat(2:2+s,:).*rt(q,:);
            end

         else
            rtnv=cproot(val,val2,lomat(1,:),[(ElementType==4)-(ElementType==3) T]);
            rt=cproot(oldval,oldval2,lomat(1,:),[(ElementType==4)-(ElementType==3) T]);
            lomat(2:2+s,:)=lomat(2:2+s,:).*(rtnv(q,:)./rt(q,:));
            cont(ElementLocation,1:2)=[val,val2];
         end
         go_for_it=1;

      else
         errordlg('wn cannot be = 0','Message','on');

      end

   else
      errordlg('Second Order values must be numbers','Message','on');

   end

% continuous integrators/differentiators or
% discrete predictors/delays
elseif any(ElementType==[0.6,0.7]),

   switch ElementOperation,
      case('Add'),
         val = str2num(get(infmat(loc,1),'string'));
         val2 =str2num(get(infmat(loc,2),'string'));

      case('Edit'),
         val = str2num(get(infmat(loc,1),'string'));
         val2 = str2num(get(infmat(loc,2),'string'));
         oldval = cont(ElementLocation, 1);
         oldval2 = cont(ElementLocation, 2);

         set(infmat(17,1),'min',0,...
                           'max',100,...
                           'value',val,...
                           'sliderstep',[1/100, 10/100]);
         set(infmat(17,2),'min',0,...
                           'max',100,...
                           'value',val2,...
                           'sliderstep',[1/100, 10/100]);

      case('Iterate'),
         oldval = cont(ElementLocation,1);
         oldval2 = cont(ElementLocation,2);
         val = get(infmat(17,1),'value');
         val2 = get(infmat(17,2),'value');

%         slider_min=get(infmat(17,1),'min');
%         slider_max=get(infmat(17,1),'max');
%         if val <= slider_min | val >= slider_max,
%            set(infmat(17,1),'min',val*0.5,...
%                             'max',val*1.5,...
%                             'sliderstep',[0.001, 0.1]);
%         end

%         slider_min=get(infmat(17,2),'min');
%         slider_max=get(infmat(17,2),'max');
%         if val2 <= slider_min | val2 >= slider_max,
%            set(infmat(17,2),'min',val2*0.5,...
%                             'max',val2*1.5,...
%                             'sliderstep',[0.01, 0.1]);
%         end
         set(infmat(16,1),'string',num2str(val));
         set(infmat(16,2),'string',num2str(val2));

   end

   if length(val),
      if rem(val,1) == 0,
         if strcmp(ElementOperation, 'Add'),
            rtnv=cintegtr(val,lomat(1,:),T) .* cintegtr(-val2,lomat(1,:),T);
            rt = ones(size(rtnv));
            cont(2+(T > 0),1)=val+cont(2+(T > 0),1);
            cont(2+(T > 0),2)=val2+cont(2+(T > 0),2);
            ElementLocation = 2 + (T > 0);

         else
            rtnv=cintegtr(val,lomat(1,:),T) .* cintegtr(-val2,lomat(1,:),T);
            rt = cintegtr(oldval, lomat(1,:), T) .* cintegtr(-oldval2, lomat(1,:), T);
            cont(ElementLocation,1)=(val-oldval)+cont(ElementLocation,1);
            cont(ElementLocation,2)=(val2-oldval2)+cont(ElementLocation,2);

         end
         lomat(2:2+s,:)=lomat(2:2+s,:).*(rtnv(q,:)./rt(q,:));
         go_for_it=1;

      else
         errordlg('Value must be an integer.','Message','on');

      end

   else
      errordlg('Value must be an integer.','Message','on');

   end

elseif ElementType==0.5, % discrete integrators

   switch ElementOperation,
      case('Add'),
         val = str2num(get(infmat(loc,1),'string'));
         val2 =str2num(get(infmat(loc,2),'string'));

      case('Edit'),
         val = str2num(get(infmat(loc,1),'string'));
         val2 = str2num(get(infmat(loc,2),'string'));
         oldval = cont(ElementLocation, 1);
         oldval2 = cont(ElementLocation, 2);
         if val < 4,
            set(infmat(17,1),'min',0,...
                              'max',3,...
                              'value',val,...
                              'sliderstep',[1/3, 10/3]);
         end
         if val2 < 4,
            set(infmat(17,2),'min',0,...
                              'max',3,...
                              'value',val2,...
                              'sliderstep',[1/3, 10/3]);
         end

      case('Iterate'),
         oldval = cont(ElementLocation,1);
         oldval2 = cont(ElementLocation,2);
         val = get(infmat(17,1),'value');
         val2 = get(infmat(17,2),'value');

%         slider_min=get(infmat(17,1),'min');
%         slider_max=get(infmat(17,1),'max');
%         if val <= slider_min | val >= slider_max,
%            set(infmat(17,1),'min',val*0.5,...
%                             'max',val*1.5,...
%                             'sliderstep',[0.001, 0.1]);
%         end

%         slider_min=get(infmat(17,2),'min');
%         slider_max=get(infmat(17,2),'max');
%         if val2 <= slider_min | val2 >= slider_max,
%            set(infmat(17,2),'min',val2*0.5,...
%                             'max',val2*1.5,...
%                             'sliderstep',[0.01, 0.1]);
%         end
         if val < 4,
            set(infmat(16,1),'string',num2str(val));
         end
         if val2 < 4,
            set(infmat(16,2),'string',num2str(val2));
         end

   end

   if length(val),
      if rem(val, 1) == 0,
         if strcmp(ElementOperation,'Add') & val < 4,
            canv=dintegtr(val,lomat(1,:),T,-1);
            ca=dintegtr(cont(2,1),lomat(1,:),T,-1);
            cont(2,1) = val;
            ElementLocation = 2;
            go_for_ita = 1;

         elseif (strcmp(ElementOperation, 'Edit') | strcmp(ElementOperation, 'Iterate')) & val < 4,
            canv=dintegtr(val, lomat(1,:),T,-1);
            ca=dintegtr(oldval,lomat(1,:),T,-1);
            cont(ElementLocation,1) = val;
            go_for_ita = 1;

         else
            set(infmat(17,1),'value',3);
            go_for_ita = 0;
            errordlg('No more than 3 integrators allowed.',...
                     'Message','on');
         end

      else
         errordlg('Integrator values must be integers.','Message','on');

      end
   else
      canv = ones(size(lomat(1,:)));
      ca=ones(size(canv));
      cont(2,1) = 0;
      ElementLocation = 2;
      go_for_ita = 1;

   end

   if length(val2),
      if rem(val2, 1) == 0,
         if strcmp(ElementOperation, 'Add') & val2 < 4,
            cbnv=dintegtr(val2,lomat(1,:),T,1);
            cb=dintegtr(cont(2,2),lomat(1,:),T,1);
            cont(2,2) = val2;
            ElementLocation = 2;
            go_for_itb = 1;

         elseif (strcmp(ElementOperation, 'Edit') | strcmp(ElementOperation, 'Iterate')) & val2 < 4,
            cbnv=dintegtr(val2, lomat(1,:),T,1);
            cb=dintegtr(oldval2,lomat(1,:),T,1);
            cont(ElementLocation,2) = val2;
            go_for_itb = 1;

         else
            set(infmat(17,2),'value',3);
            go_for_itb = 0;
            errordlg('No more than 3 differentiators allowed.',...
                     'Message','on');

         end

      else
         errordlg('Differentiator values must be integers.','Message','on');

      end

   else
      cbnv = ones(size(lomat(1,:)));
      cb=ones(size(cbnv));
      cont(2,2) = 0;
      ElementLocation = 2;
      go_for_itb = 1;

   end

   go_for_it = go_for_ita * go_for_itb;

   if go_for_it,
      cpnv = canv.*cbnv;
      cp = ca.*cb;
      lomat(2:2+s,:)=lomat(2:2+s,:).*(cpnv(q,:)./cp(q,:));
   end

elseif ElementType==5, % lead/lag

   switch ElementOperation,
      case('Add'),
         val=str2num(get(infmat(loc,1),'string'));
         val2=str2num(get(infmat(loc,2),'string'));

      case('Edit'),
         val=str2num(get(infmat(loc,1),'string'));
         val2=str2num(get(infmat(loc,2),'string'));

         set(infmat(17,1),'min',val*0.5,...
                           'max',val*1.5,...
                           'value',val,...
                           'sliderstep',[0.01, 0.1]);

         if val2 > 0,
            set(infmat(17,2),'min',val2*0.5,...
                              'max',val2*1.5,...
                              'value',val2,...
                              'sliderstep',[0.01, 0.1]);
         else
            set(infmat(17,2),'min',val2*1.5,...
                              'max',val2*0.5,...
                              'value',val2,...
                              'sliderstep',[0.01, 0.1]);
         end;%if

         oldval = cont(ElementLocation, 1);
         oldval2 = cont(ElementLocation, 2);

      case('Iterate'),
         oldval = cont(ElementLocation, 1);
         oldval2 = cont(ElementLocation, 2);
         val=get(infmat(17,1),'value');
         val2=get(infmat(17,2),'value');

         slider_min=get(infmat(17,1),'min');
         slider_max=get(infmat(17,1),'max');
         if val <= slider_min | val >= slider_max,
            set(infmat(17,1),'min',val*0.5,...
                             'max',val*1.5,...
                             'sliderstep',[0.01, 0.1]);
         end

         slider_min=get(infmat(17,2),'min');
         slider_max=get(infmat(17,2),'max');
         if val2 <= slider_min | val2 >= slider_max,
            if val2 > 0,
               set(infmat(17,2),'min',val2*0.5,...
                                'max',val2*1.5,...
                                'sliderstep',[0.01, 0.1]);
            else
               set(infmat(17,2),'min',val2*1.5,...
                                 'max',val2*0.5,...
                                 'sliderstep',[0.01, 0.1]);
            end
         end
         set(infmat(16,1),'string',num2str(val));
         set(infmat(16,2),'string',num2str(val2));

   end

   if length([val,val2]) == 2,
      if val~=0 & abs(val)<88,
         rtnv=ldlgcplx(val,val2,lomat(1,:),T);
         if strcmp(ElementOperation,'Add'),
            rt = ones(size(rtnv));
            cont=[cont;val,val2,NaN,5];
            ElementLocation = size(cont,1);

         else
            rt=ldlgcplx(oldval,oldval2,lomat(1,:),T);
            cont(ElementLocation,1:2)=[val,val2];

         end
         lomat(2,:)=lomat(2,:).*(rtnv./rt);
         go_for_it=1;
      else
         if length(val) & val==0,
            errordlg('w cannot be = 0','Message','on');
         elseif length(val2) & abs(val2)>=88,
            errordlg('Phase change must be < 88 degrees','Message','on');
         end
      end
   else
      errordlg('Phase and Frequency must be numbers','Message','on');
   end

elseif ElementType==6, % Notch

   switch ElementOperation,
      case('Add'),
         val=str2num(get(infmat(loc,1),'string'));
         val2=str2num(get(infmat(loc,2),'string'));
         val3=str2num(get(infmat(loc,3),'string'));

      case('Edit'),
         val=str2num(get(infmat(loc,1),'string'));
         val2=str2num(get(infmat(loc,2),'string'));
         val3=str2num(get(infmat(loc,3),'string'));

         set(infmat(17,1),'min',val*0.5,...
                           'max',val*1.5,...
                           'value',val,...
                           'sliderstep',[0.001, 0.1]);
         set(infmat(17,2),'min',val2*0.5,...
                           'max',val2*1.5,...
                           'value',val2,...
                           'sliderstep',[0.001, 0.1]);
         set(infmat(17,3),'min',val3*0.5,...
                           'max',val3*1.5,...
                           'value',val3,...
                           'sliderstep',[0.01, 0.1]);

         oldval = cont(ElementLocation, 1);
         oldval2 = cont(ElementLocation, 2);
         oldval3 = cont(ElementLocation, 3);

      case('Iterate'),
         oldval = cont(ElementLocation, 1);
         oldval2 = cont(ElementLocation, 2);
         oldval3 = cont(ElementLocation, 3);
         val=get(infmat(17,1),'value');
         val2=get(infmat(17,2),'value');
         val3=get(infmat(17,3),'value');

         slider_min=get(infmat(17,1),'min');
         slider_max=get(infmat(17,1),'max');
         if val <= slider_min | val >= slider_max,
            set(infmat(17,1),'min',val*0.5,...
                             'max',val*1.5,...
                             'sliderstep',[0.001, 0.1]);
         end

         slider_min=get(infmat(17,2),'min');
         slider_max=get(infmat(17,2),'max');
         if val2 <= slider_min | val2 >= slider_max,
            set(infmat(17,2),'min',val2*0.5,...
                             'max',val2*1.5,...
                             'sliderstep',[0.001, 0.1]);
         end

         slider_min=get(infmat(17,3),'min');
         slider_max=get(infmat(17,3),'max');
         if val3 <= slider_min | val3 >= slider_max,
            set(infmat(17,3),'min',val3*0.5,...
                             'max',val3*1.5,...
                             'sliderstep',[0.01, 0.1]);
         end

         set(infmat(16,1),'string',num2str(val));
         set(infmat(16,2),'string',num2str(val2));
         set(infmat(16,3),'string',num2str(val3));

   end

   if length([val,val2,val3])==3,
      ztas=[val,val2];
      freq=val3;
      delay=infmat(10,1);
      zta=ztas(1+(ztas(1)>ztas(2)));
      if freq~=0,
         if strcmp(ElementOperation,'Add'),
            cont=[cont;ztas(1),ztas(2),freq,6];
            ElementLocation = size(cont,1);
            nom=get(bthan(2),'userdata');
            if length(nom) & ElementEnvironment~=2,
               wlt=qfrqenh(freq,zta,lomat(1,:),T);
               if ElementEnvironment==1,
                  nlo=squeeze(freqresp(nom,wlt)).';
%                  nlo=qcpqft(nom(1,:),nom(2,:),wlt,T);
%                  clo=qcntbode(cont,wlt,T).*exp(-i*wlt*delay);
                  clo=qcntbode(cont,wlt,T);
                  lo=nlo.*clo;
               else
                  lo=qcntbode(cont,wlt,T).*exp(-i*wlt*delay);
               end
               lomat=[wlt;lo;ones(1,length(lo))];
            else
               rt=ntchcplx(ztas(1),ztas(2),freq,lomat(1,:),T);
               lomat(2:2+s,:)=lomat(2:2+s,:).*rt(q,:);
            end
         else
            rtnv=ntchcplx(ztas(1),ztas(2),freq,lomat(1,:),T);
            rt=ntchcplx(cont(ElementLocation,1),cont(ElementLocation,2),cont(ElementLocation,3),lomat(1,:),T);
            lomat(2:2+s,:)=lomat(2:2+s,:).*(rtnv(q,:)./rt(q,:));
            cont(ElementLocation,1:3)=[ztas(1),ztas(2),freq];

         end
         go_for_it=1;

      else
         errordlg('wn cannot be = 0','Message','on');

      end

   else
      errordlg('Notch values must be numbers','Message','on');

   end

elseif ElementType==7, % Complex Lead

   switch ElementOperation,
      case('Add'),
         val=str2num(get(infmat(loc,1),'string'));
         val2=str2num(get(infmat(loc,2),'string'));
         val3=str2num(get(infmat(loc,3),'string'));

      case('Edit'),
         val=str2num(get(infmat(loc,1),'string'));
         val2=str2num(get(infmat(loc,2),'string'));
         val3=str2num(get(infmat(loc,3),'string'));

         if val > 0,
            set(infmat(17,1),'min',val*0.5,...
                              'max',val*1.5,...
                              'value',val,...
                              'sliderstep',[0.01, 0.1]);
         else
            set(infmat(17,1),'min',val*1.5,...
                              'max',val*0.5,...
                              'value',val,...
                              'sliderstep',[0.01, 0.1]);
         end;%if

         set(infmat(17,2),'min',val2*0.5,...
                           'max',val2*1.5,...
                           'value',val2,...
                           'sliderstep',[0.01, 0.1]);
         set(infmat(17,3),'min',val3*0.5,...
                           'max',val3*1.5,...
                           'value',val3,...
                           'sliderstep',[0.01, 0.1]);

         oldval = cont(ElementLocation, 1);
         oldval2 = cont(ElementLocation, 2);
         oldval3 = cont(ElementLocation, 3);

      case('Iterate'),
         oldval = cont(ElementLocation, 1);
         oldval2 = cont(ElementLocation, 2);
         oldval3 = cont(ElementLocation, 3);
         val=get(infmat(17,1),'value');
         val2=get(infmat(17,2),'value');
         val3=get(infmat(17,3),'value');

         slider_min=get(infmat(17,1),'min');
         slider_max=get(infmat(17,1),'max');
         if val <= slider_min | val >= slider_max,
            if val > 0,
               set(infmat(17,1),'min',val*0.5,...
                                 'max',val*1.5,...
                                 'sliderstep',[0.01, 0.1]);
            else
               set(infmat(17,1),'min',val*1.5,...
                                 'max',val*0.5,...
                                 'sliderstep',[0.01, 0.1]);
            end;%if
         end

         slider_min=get(infmat(17,2),'min');
         slider_max=get(infmat(17,2),'max');
         if val2 <= slider_min | val2 >= slider_max,
            set(infmat(17,2),'min',val2*0.5,...
                             'max',val2*1.5,...
                             'sliderstep',[0.01, 0.1]);
         end

         slider_min=get(infmat(17,3),'min');
         slider_max=get(infmat(17,3),'max');
         if val3 <= slider_min | val3 >= slider_max,
            set(infmat(17,3),'min',val3*0.5,...
                             'max',val3*1.5,...
                             'sliderstep',[0.01, 0.1]);
         end

         set(infmat(16,1),'string',num2str(val));
         set(infmat(16,2),'string',num2str(val2));
         set(infmat(16,3),'string',num2str(val3));

   end

   if length([val,val2,val3])==3,
      delph = val;
      delmg = val2;
      freq = val3;
      delay = infmat(10,1);
      if freq~=0,
         if strcmp(ElementOperation,'Add'),
            cont=[cont;delph,delmg,freq,7];
            ElementLocation = size(cont,1);
            nom=get(bthan(2),'userdata');
 %           if length(nom) & ElementEnvironment~=2,
 %              wlt=qfrqenh(freq,zta,lomat(1,:),T);
 %              if ElementEnvironment==1,
 %                 nlo=qcpqft(nom(1,:),nom(2,:),wlt,T);
 %                 clo=qcntbode(cont,wlt,T).*exp(-i*wlt*delay);
 %                 lo=nlo.*clo;
 %              else
 %                 lo=qcntbode(cont,wlt,T).*exp(-i*wlt*delay);
 %              end
 %              lomat=[wlt;lo;ones(1,length(lo))];
 %           else
               rt=complexlead(delph,delmg,freq,lomat(1,:),T);
               lomat(2:2+s,:)=lomat(2:2+s,:).*rt(q,:);
 %           end
         else
            rtnv=complexlead(delph,delmg,freq,lomat(1,:),T);
            rt=complexlead(cont(ElementLocation,1),cont(ElementLocation,2),cont(ElementLocation,3),lomat(1,:),T);
            lomat(2:2+s,:)=lomat(2:2+s,:).*(rtnv(q,:)./rt(q,:));
            cont(ElementLocation,1:3)=[delph,delmg,freq];

         end
         go_for_it=1;

      else
         errordlg('w cannot be = 0','Message','on');

      end

   else
      errordlg('Complex Lead values must be numbers','Message','on');

   end

elseif ElementType == -1,
   go_for_it = 1;

end

if go_for_it,
   if strcmp(ElementOperation,'Done'),
      set(bthan(19),'userdata',[]);
      set(bthan(20),'userdata',[]);
      set(bthan(1),'userdata',lomat);
      set(bthan(3),'userdata',cont);

      set(bthan([1,8,29]),'enable','on');
      v=get(bthan(10),'userdata');
      v2=get(bthan(21),'userdata');
      if infmat(9,1)==1,
         vo2=get(bthan(22),'userdata');
         vo=get(bthan(17),'userdata');
         set([v,vo],'xdata',0,'ydata',0);
         set(bthan(22),'userdata',vo);
         set(bthan(17),'userdata',vo2);
         set(v2,'linestyle','-');
         set(v,'linestyle',':');
         set(bthan(10),'userdata',v2);
         set(bthan(21),'userdata',v);
         set(infmat(8,1),'enable','on');
      else
         if ElementEnvironment==1, qnicplt(f);
         elseif ElementEnvironment==2, qmagplt(f);
         elseif ElementEnvironment==3, mgphplot(f);
         end
      end
      [ControllerString,ListboxInfo] = cntstr(f,cont);
      set(ElementsListbox, 'string', ControllerString,...
                           'userdata', ListboxInfo);
%      qelmtpopup('Done');
      qelmtlistbox('Iterate');
      qclswin(0);

   elseif strcmp(ElementOperation,'Cancel'),
      v2=get(bthan(21),'userdata');
      vo2=get(bthan(22),'userdata');
      if infmat(9,1)==1,
         set([v2(:);vo2(:)],'vis','off');
      else
         set(v2,'xdata',0,'ydata',0);
      end
      cont = get(bthan(3),'userdata');
      lomat = get(bthan(1),'userdata');
      set(bthan(19),'userdata',[]);
      set(bthan(20),'userdata',[]);
      [ControllerString,ListboxInfo] = cntstr(f,cont);
      ListboxValue = get(ElementsListbox,'value');
      if ListboxValue > length(ControllerString),
         ListboxValue = length(ControllerString);
      end
      set(ElementsListbox, 'string', ControllerString,...
                           'userdata', ListboxInfo,...
                           'value',ListboxValue);

% go into edit/iterate mode
      qelmtlistbox('Edit');
      qclswin(0);

   else
      set(bthan(19),'userdata',cont);
      set(bthan(20),'userdata',lomat);
      if ElementEnvironment==1, qnicplt(f);
      elseif ElementEnvironment==2, qmagplt(f);
      elseif ElementEnvironment==3, mgphplot(f);
      end
      set(infmat(13,2),'callback','qelmts(-1,0,''Done'')','enable','on');
      set(infmat(13,3),'callback','qelmts(-1,0,''Cancel'')','enable','on');
      [ControllerString,ListboxInfo,ListboxValue] = cntstr(f,cont,ElementLocation);
      set(ElementsListbox, 'string', ControllerString,...
                           'userdata', ListboxInfo,...
                           'value',ListboxValue);

      if kw ~= 0,
         infmat(4,3)=rtnv(kw);
         set(bthan(16),'userdata',infmat);
      end

      if strcmp(ElementOperation, 'Add'),
         qelmtlistbox('Add');
      end
   end

end
