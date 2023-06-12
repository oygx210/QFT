function plottmp(w,P,nom,pos)
% PLOTTMPL Plot frequency response arrays (QFT templates).
%          PLOTTMPL(W,P,NOM) plots the frequency response array
%          defined by the LTI/FRD array P at the frequencies W.
%          NOM designates the location of the nominal plant
%          (for display purposes).
%
%          PLOTTMPL(W,P) uses the default value for NOM (the first array element).
%
%          See also PLOTBNDS, SISOBNDS, GENBNDS, SECTBNDS.

% Author: Craig Borghesani
% Date: 9/6/93
% Revised: 2/17/96 10:00 PM V1.1 updates.
%          8/28/02 1:45PM V2.1 updates.
%          4/10/03 3:29PM V2.5 updates.
% Copyright (c) 2003, Terasoft, Inc.

% bound button initial locations
lf1=0.0781;
ht=0.0417;
tp=1-0.0625;

if ~isstr(w),

% remove wbd
   wbd = [];

   cvec=['r';'g';'b';'c';'m'];
   clr=[cvec;cvec;cvec;cvec];
   clr=[clr;clr;clr;clr;clr];
   clr=[clr;clr;clr;clr;clr];
   fpos=[0.333,0.28,0.6620,0.6604];
   if nargin==2,
    nom=1;
   elseif nargin==4,
    fpos=pos;
    if ~length(nom), nom=1; end
   end

   if ~length(wbd), wbd=w; end

   [jk,jk,jk,mag,ph]=bndsdef(w,[],1,P,[],[],[],[],[],0);
   sbds=qsubset(wbd,w);
   [jk,ind]=sort(wbd); sbds=sbds(ind);

   ph=ph*180/pi;
   mag=20*log10(mag);

   f1 = figure('name','Plot Template','numbertitle','off','units','norm',...
          'pos',fpos,'vis','off','tag','qft');

   a=gca;
   apos=get(a,'position');
   set(a,'box','on','xgrid','on','ygrid','on',...
         'gridlinestyle',':',...
         'nextplot','add','xlim',[-360,0],'ylim',[min(mag(:))-5,max(mag(:))+5]);

   UICont = uicontextmenu;

   c=1; tmpdata2=[]; h3=0;
   for k=sbds,
    h1=plot(ph(:,k),mag(:,k),['o',clr(c)]);
    h2=plot(ph(nom,k),mag(nom,k),'w*');
    h3=text('pos',[ph(nom,k)+5,mag(nom,k)],'string',num2str(w(k)),...
             'horizontalalignment','left');
    wstr=sprintf('%4.4g',w(k));
    loc=find(wstr=='e');
    if length(loc),
     wstr(find(wstr(loc:length(wstr))=='0')+(loc-1))=[];
    else
     wstr(find(wstr==' '))=[];
    end
    st2=[wstr,', ',clr(c)];
    uimenu(UICont,'label',st2,'callback','plottmpl(''Toggle'')',...
                  'userdata', k, 'checked', 'on');
    st=[wstr];
    tmpdata=[h1(:)',h2,h3];
    tmpdata2=[tmpdata2;tmpdata];
    tmphandle(c) = tmpdata(1);
    tmplabel{c} = st;
    tp=tp-ht;
    c=c+1;
   end
   uimenu(UICont,'label','Toggle Off','callback','plottmpl(''Toggle All'')',...
                 'separator','on');
   set(a,'uicontextmenu',UICont,'userdata',tmpdata2);

   legend(tmphandle, tmplabel);
   xlabel('Open-Loop Phase (deg)');
   ylabel('Open-Loop Gain (dB)');
   xlim = get(a,'xlim'); ylim = get(a,'ylim');
   if diff(xlim) >= 360,
      XTickValues = qaxesadjust(xlim);
      set(gca,'xtick',XTickValues);
   end

   hToolbar = findall(f1, 'type', 'uitoolbar');
   hChil = allchild(hToolbar);
   tooltips = get(hChil, 'tooltip');
   zoomIn = strmatch('Zoom In', tooltips, 'exact');
   zoomOut = strmatch('Zoom Out', tooltips, 'exact');
   delete(hChil(zoomOut));
   set(hChil(zoomIn),'clickedcallback','plottmpl(''Zoom Set'')');

   set(f1,'vis','on','userdata',[xlim,ylim]);

elseif strcmp(w, 'Toggle'),

   tmpdata2 = get(gca,'userdata');
   CurrentItem = get(gcbo,'userdata');

   if strcmp(get(tmpdata2(CurrentItem,1),'vis'),'on'),
      set(tmpdata2(CurrentItem,:),'vis','off');
      set(gcbo,'checked','off');
   else
      set(tmpdata2(CurrentItem,:),'vis','on');
      set(gcbo,'checked','on');
   end;%if

elseif strcmp(w, 'Toggle All'),

   CurrentMenu = gcbo;
   CurrentParent = get(CurrentMenu,'parent');
   CurrentChil = flipud(get(CurrentParent,'children'));
   tmpdata2 = get(gca,'userdata');

   if strcmp(get(CurrentMenu,'label'), 'Toggle Off'),
      set(CurrentMenu,'label','Toggle On');
      for k = 1:length(CurrentChil)-1,
         CurrentItem = get(CurrentChil(k),'userdata');
         set(CurrentChil(k),'checked','off');
         set(tmpdata2(k,:),'vis','off');
      end; %for


   else
      set(CurrentMenu,'label','Toggle Off');
      for k = 1:length(CurrentChil)-1,
         CurrentItem = get(CurrentChil(k),'userdata');
         set(CurrentChil(k),'checked','on');
         set(tmpdata2(k,:),'vis','on');
      end; %for

   end

elseif strcmp(w, 'Zoom Set'),

   f = gcf;
   if strcmp(get(f,'windowbuttondownfcn'),'plottmpl(''Zoom'')'),
      set(f,'windowbuttondownfcn','');
      set(f,'pointer','arrow');
      txt = findobj(f,'tag','uicontextholder');
      UICont = get(txt,'uicontextmenu');
      set(gca,'uicontextmenu',UICont);
      delete(txt);

   else
      txt = text(0,0,' ','vis','off','tag','uicontextholder');
      UICont = get(gca,'uicontextmenu');
      set(txt,'uicontextmenu',UICont);
      set(gca,'uicontextmenu',[]);
      set(f,'windowbuttondownfcn','plottmpl(''Zoom'')');
      setptr(f,'glass');

   end;%if

elseif strcmp(w,'Zoom'),

   f = gcf;
   am = gca;

   FigureSelectionType = get(f,'selectiontype');

   CurrentLimits1 = [get(am,'xlim'), get(am,'ylim')];

   switch FigureSelectionType,

      case('normal'),

         CurrentPoint1Axis1 = get(am,'currentpoint');
         fUnits = get(f,'units');
         set(f,'units','pixel');
         rbbox([get(f,'currentpoint') 0 0],...
                  get(f,'currentpoint'));
         set(f,'units',fUnits);
         CurrentPoint2Axis1 = get(am,'currentpoint');

         CurrentDiff = sum(abs(CurrentPoint1Axis1(1,:) - CurrentPoint2Axis1(1,:)));

% store the current limits
         ZoomMemoryData = getappdata(f, 'ZoomMemory');
         if ~isempty(ZoomMemoryData),
            ZoomMemoryData = [ZoomMemoryData; CurrentLimits1];

         else
            ZoomMemoryData = CurrentLimits1;

         end
         setappdata(f,'ZoomMemory',ZoomMemoryData);

         if CurrentDiff > 0.1,
            XLimits1 = [min(CurrentPoint1Axis1(1,1),CurrentPoint2Axis1(1,1)),...
                        max(CurrentPoint1Axis1(1,1),CurrentPoint2Axis1(1,1))];

            YLimits1 = [min(CurrentPoint1Axis1(1,2),CurrentPoint2Axis1(1,2)),...
                        max(CurrentPoint1Axis1(1,2),CurrentPoint2Axis1(1,2))];

         else

% lets zoom in by 50% around the selected point
            DiffXLimits1 = diff(CurrentLimits1(1:2))/2;
            DiffYLimits1 = diff(CurrentLimits1(3:4))/2;

            XLimits1 = CurrentPoint1Axis1(1,1) + [-DiffXLimits1/2, DiffXLimits1/2];
            YLimits1 = CurrentPoint1Axis1(1,2) + [-DiffYLimits1/2, DiffYLimits1/2];

         end

         if diff(XLimits1) >= 360,
            XTickValues = qaxesadjust(XLimits1);
            set(am, 'xlim', XLimits1, 'ylim', YLimits1,...
                     'xtick',XTickValues);

         else
            set(am, 'xlim', XLimits1, 'ylim', YLimits1,...
                     'xtickmode','auto');

         end

      case('alt'),

         ZoomMemoryData = getappdata(f,'ZoomMemory');
         if length(ZoomMemoryData),
            XLimits1 = ZoomMemoryData(end, 1:2);
            YLimits1 = ZoomMemoryData(end, 3:4);
            ZoomMemoryData(end, :) = [];
            setappdata(f, 'ZoomMemory', ZoomMemoryData);
            if diff(XLimits1) >= 360,
               XTickValues = qaxesadjust(XLimits1);
               set(am, 'xlim', XLimits1, 'ylim', YLimits1,...
                        'xtick',XTickValues);

            else
               set(am, 'xlim', XLimits1, 'ylim', YLimits1,...
                        'xtickmode','auto');

            end
         end

   end % FigureSelectionType

end
