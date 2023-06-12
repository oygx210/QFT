function qmouse(UIOperation)

f = gcf;
if isempty(findstr(get(f,'name'),'Shaping')),
   return;
end

bthan = get(f,'userdata');
infmat = get(bthan(16),'userdata');
QFTToolData = getappdata(f,'QFTToolData');

if infmat(9,1)==3,
   axs=infmat(2,3:4);
else
   axs=infmat(1,:);
end
am = infmat(24,1);
ap = infmat(24,2);
txt1 = infmat(30,1);
txt2 = infmat(30,2);
txt3 = infmat(30,3);
kw = infmat(3,3);
v2 = get(bthan(21), 'userdata');
zoomMenuItem = findobj(f,'tag','zoommode');

PopupIndex = get(infmat(16,4),'value');
PopupString = get(infmat(16,4),'string');

% are we over the axes? display current location of mouse.
pointer_loc = get(f,'currentpoint');
axs_pos = get(am,'position');

if all(pointer_loc > [axs_pos(1), axs_pos(2)]) & ...
   all(pointer_loc < [sum(axs_pos([1,3])), sum(axs_pos([2,4]))]),

   switch UIOperation,
      case('Floating'),

         lomat=get(bthan(20),'userdata');
         if isempty(lomat),
            lomat = get(bthan(1),'userdata');
         end
         w=lomat(1,:);
         lo=lomat(2,:);

         FrequencyUnitStr = 'rad/sec';
         FrequencyUnit = 1;
         if infmat(9,1) == 1,
            if ~get(QFTToolData.Mouse(6),'value'),
               FrequencyUnitStr = 'hz';
               FrequencyUnit = 1/(2*pi);
            end
         end

         if any(infmat(9,1)==[1,2]),
            pt=get(am,'currentpoint');
            wloc=qfindfrq(pt(1,1),pt(1,2),lomat,am);
         else
            am_pos = get(am,'pos');
            ap_pos = get(ap,'pos');
            dist = am_pos(2) - sum(ap_pos([2,4]));
            middle = am_pos(2) - dist/2;
            if pointer_loc(2) >= middle, % Magnitude subplot
               pt=get(am,'currentpoint');
               wloc=qfindfrq(pt(1,1),pt(1,2),lomat,am);

            else % Phase subplot
               pt=get(ap,'currentpoint');
               wloc=qfindfrq(pt(1,1),pt(1,2),lomat,ap,1);

            end
         end

         if length(wloc),
            set(f,'windowbuttondownfcn','');
            if wloc == kw,
               butn_dwn = get(infmat(6,1),'userdata');
               set(f,'windowbuttondownfcn',butn_dwn);
               switch PopupString{PopupIndex},
                  case('Gain'),
                     setptr(f,'uddrag');

                  case({'Real Pole', 'Real Zero'}),
                     if infmat(9,1) == 1,
                        setptr(f,'lrdrag');
                     else
                        setptr(f,'uddrag');
                     end

                  case({'Complex Pole', 'Complex Zero'}),
                     if infmat(9,1) == 1,
                        setptr(f,'fleur');
                     else
                        setptr(f,'uddrag');
                     end

                  case('Lead or Lag'),
                     setptr(f,'lrdrag');

                  case('Notch'),
                     setptr(f,'uddrag');

                  case('Super 2nd'),
                     setptr(f,'fleur');

                  case('Complex Lead/Lag'),
                     setptr(f,'lrdrag');

               end

            elseif kw > 0,
               if strcmp(get(zoomMenuItem,'checked'),'on'),
                  set(f,'windowbuttondownfcn','qmouse(''Zoom'')');
                  setptr(f,'glass');
               else
                  set(f,'windowbuttondownfcn','');
                  setptr(f,'arrow');
               end

            elseif kw == 0,
               if infmat(9,1)==1,
                  set(infmat(6,1),'xdata',qfixfase(lo,axs,wloc),...
                                  'ydata',20*log10(abs(lo(1,wloc))),...
                                  'vis','on');
               elseif infmat(9,1)==2,
                  set(infmat(6,1),'xdata',w(wloc),'ydata',20*log10(abs(lo(1,wloc))), ...
                                  'vis','on');
               end

               switch PopupString{PopupIndex},
                  case('Gain'),
                     setptr(f,'uddrag');
                     set(f,'windowbuttondownfcn','mogain(0,0)');

                  case({'Real Pole', 'Real Zero'}),
                     if infmat(9,1) == 1,
                        setptr(f,'lrdrag');
                     else
                        setptr(f,'uddrag');
                     end
                     set(f,'windowbuttondownfcn','mofirst(0,0)');

                  case({'Complex Pole', 'Complex Zero'}),
                     if infmat(9,1) == 1,
                        setptr(f,'fleur');
                     else
                        setptr(f,'uddrag');
                     end
                     set(f,'windowbuttondownfcn','mosecond(0,0)');

                  case('Lead or Lag'),
                     setptr(f,'lrdrag');
                     set(f,'windowbuttondownfcn','moldlg(0,0)');

                  case('Notch'),
                     setptr(f,'uddrag');
                     set(f,'windowbuttondownfcn','montch(0,0)');

                  case('Super 2nd'),
                     setptr(f,'fleur');
                     set(f,'windowbuttondownfcn','mo2ovr2(0)');

                  case('Complex Lead/Lag'),
                     setptr(f,'lrdrag');
                     set(f,'windowbuttondownfcn','mocpld(0,0)');
               end

            end

         else
            if strcmp(get(zoomMenuItem,'checked'),'on'),
               set(f,'windowbuttondownfcn','qmouse(''Zoom'')');
               setptr(f,'glass');
            else
               set(f,'windowbuttondownfcn','');
               setptr(f,'arrow');
            end
            if kw == 0, set(infmat(6,1),'vis','off'); end

         end

         if infmat(9,1)==1,  % SHAPE, DSHAPE
            db_ol=sprintf('%0.2f',pt(1,2));
            ph_ol=sprintf('%0.2f',pt(1,1));
            mag=10^(pt(1,2)/20);
            ph=pi/180*pt(1,1);
            ol=mag*exp(i*ph);
            db_cl=sprintf('%0.2f',20*log10(abs(ol/(1+ol))));
            ph_cl=sprintf('%0.2f',180/pi*qatan4(ol/(1+ol)));
            set(txt1,'string',['Open-Loop: ',ph_ol,' deg, ',db_ol,' dB']);
            set(txt2,'string',['Closed-Loop: ',ph_cl,' deg, ',db_cl,' dB']);

            if length(wloc),
               set(txt3,'string',sprintf('Frequency: %0.3g %s',w(wloc)*FrequencyUnit, FrequencyUnitStr));
            else
               set(txt3,'string',['Frequency: n/a ',FrequencyUnitStr]);
            end

         elseif infmat(9,1)==2, % FSHAPE, DFSHAPE
            db_ol=sprintf('%0.2f',pt(1,2));
            fr_ol=sprintf('%0.3g',pt(1,1)*FrequencyUnit);
            if length(wloc),
               magn = 20*log10(abs(lomat(2,wloc)));
               phn = 180/pi*qatan4(lomat(2,wloc));
               db_ol=sprintf('%0.2f',magn);
               ph_ol=sprintf('%0.2f',phn);
               mag=10^(magn/20);
               ph=pi/180*phn;
               ol=mag*exp(i*ph);
               db_cl=sprintf('%0.2f',20*log10(abs(ol/(1+ol))));
               ph_cl=sprintf('%0.2f',180/pi*qatan4(ol/(1+ol)));
               set(txt1,'string',['Open-Loop: ',ph_ol,' deg, ',db_ol,' dB']);
               set(txt2,'string',['Closed-Loop: ',ph_cl,' deg, ',db_cl,' dB']);
            else
               set(txt1,'string',['Open-Loop: ',db_ol,' dB']);
               set(txt2,'string','');
            end
            set(txt3,'string',['Frequency: ',fr_ol,' ',FrequencyUnitStr]);

         elseif infmat(9,1)==3, % BODPLOT/DBODPLOT

            if pointer_loc(2) >= middle, % Magnitude subplot
               db_ol=sprintf('%0.2f',pt(1,2));
               fr_ol=sprintf('%0.2f',pt(1,1)*FrequencyUnit);
               if length(wloc),
                  magn = 20*log10(abs(lomat(2,wloc)));
                  phn = 180/pi*qatan4(lomat(2,wloc));
                  db_ol=sprintf('%0.2f',magn);
                  ph_ol=sprintf('%0.2f',phn);
                  mag=10^(magn/20);
                  ph=pi/180*phn;
                  ol=mag*exp(i*ph);
                  db_cl=sprintf('%0.2f',20*log10(abs(ol/(1+ol))));
                  ph_cl=sprintf('%0.2f',180/pi*qatan4(ol/(1+ol)));
                  set(txt1,'string',[ph_ol,'deg, ',db_ol,'dB']);
                  set(txt2,'string',[ph_cl,'deg, ',db_cl,'dB']);
               else
                  set(txt1,'string',[db_ol,'dB']);
                  set(txt2,'string','');
               end
               set(txt3,'string',['Frequency: ',fr_ol,' ',FrequencyUnitStr]);

            else % Phase subplot
               ph_ol=sprintf('%0.2f',pt(1,2));
               fr_ol=sprintf('%0.2f',pt(1,1)*FrequencyUnit);
               if length(wloc),
                  magn = 20*log10(abs(lomat(2,wloc)));
                  phn = 180/pi*qatan4(lomat(2,wloc));
                  db_ol=sprintf('%0.2f',magn);
                  ph_ol=sprintf('%0.2f',phn);
                  mag=10^(magn/20);
                  ph=pi/180*phn;
                  ol=mag*exp(i*ph);
                  db_cl=sprintf('%0.2f',20*log10(abs(ol/(1+ol))));
                  ph_cl=sprintf('%0.2f',180/pi*qatan4(ol/(1+ol)));
                  set(txt1,'string',[ph_ol,'deg, ',db_ol,'dB']);
                  set(txt2,'string',[ph_cl,'deg, ',db_cl,'dB']);
               else
                  set(txt1,'string',[ph_ol,176]);
                  set(txt2,'string','');
               end
               set(txt3,'string',[fr_ol,' ',FrequencyUnitStr]);
            end
         end

      case('Zoom'),

         FigureSelectionType = get(f,'selectiontype');

         CurrentLimits1 = [get(am,'xlim'), get(am,'ylim')];

         switch FigureSelectionType,

            case('normal'),

               CurrentPoint1Axis1 = get(am,'currentpoint');
               rbbox([get(f,'currentpoint') 0 0],...
                      get(f,'currentpoint'));
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

               if infmat(9,1)==1,
                  if diff(XLimits1) >= 360,
                     XTickValues = qaxesadjust(XLimits1);
                     set(am, 'xlim', XLimits1, 'ylim', YLimits1,...
                             'xtick',XTickValues);

                  else
                     set(am, 'xlim', XLimits1, 'ylim', YLimits1,...
                             'xtickmode','auto');

                  end
               else
                  set(am, 'xlim', XLimits1, 'ylim', YLimits1);

               end

            case('alt'),

               ZoomMemoryData = getappdata(f,'ZoomMemory');
               if length(ZoomMemoryData),
                  XLimits1 = ZoomMemoryData(end, 1:2);
                  YLimits1 = ZoomMemoryData(end, 3:4);
                  ZoomMemoryData(end, :) = [];
                  setappdata(f, 'ZoomMemory', ZoomMemoryData);
                  if infmat(9,1)==1,
                     if diff(XLimits1) >= 360,
                        XTickValues = qaxesadjust(XLimits1);
                        set(am, 'xlim', XLimits1, 'ylim', YLimits1,...
                                'xtick',XTickValues);

                     else
                        set(am, 'xlim', XLimits1, 'ylim', YLimits1,...
                                'xtickmode','auto');

                     end
                  else
                     set(am, 'xlim', XLimits1, 'ylim', YLimits1);

                  end
               end

         end % FigureSelectionType

   end

else
   setptr(f,'arrow');

   if kw == 0,
      set(infmat(6,1),'vis','off');
   end

end

