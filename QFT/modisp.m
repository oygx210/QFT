function modisp
% MODISP Display legend data. (Utility Function)
%        MODISP takes the current mouse location and determines its
%        plot location within the IDE.  The results are then displayed in
%        the upper-right hand corner of the screen.

% Author: Craig Borghesani
% 9/5/93
% Copyright (c) 2003, Terasoft, Inc.

do_it = 1;
pointer_figure = get(0,'pointerwindow');
temp=get(pointer_figure,'userdata');
size_temp=size(temp);
if all(size_temp==[1,1]),
 f=temp;
 bthan=get(f,'userdata');
 if ~all(size(bthan)==[1,40]),
  do_it = 0;
 end
elseif all(size_temp==[1,40]),
 f=pointer_figure;
 bthan=temp;
else
 do_it = 0;
end

if do_it,
 pointer_loc = get(f,'currentpoint');
 QFTToolData = getappdata(f,'QFTToolData');
 infmat=get(bthan(16),'userdata');
 lomat=get(bthan(20),'userdata');
 if isempty(lomat),
   lomat = get(bthan(1),'userdata');
 end
 w=lomat(1,:); lo=lomat(2,:);
 if infmat(9,1)==3, axs=infmat(2,3:4);
 else axs=infmat(1,:); end
 am = infmat(24,1);
 ap = infmat(24,2);
 txt1 = infmat(30,1);
 txt2 = infmat(30,2);
 txt3 = infmat(30,3);
 butn_dwn = get(infmat(6,1),'userdata');
 kw = infmat(3,3);

 if any(infmat(9,1)==[1,2]),
  pt=get(am,'currentpoint');
  wloc=qfindfrq(pt(1,1),pt(1,2),lomat,am);
 else
  am_pos = get(am,'pos'); ap_pos = get(ap,'pos');
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
  if length(butn_dwn) & wloc==kw,
   set(f,'pointer','crosshair');
   set(f,'windowbuttondownfcn',butn_dwn);
  elseif length(butn_dwn),
   set(f,'pointer','arrow');
   set(f,'windowbuttondownfcn','');
  else
   set(f,'pointer','crosshair');
  end
 else
  set(f,'pointer','arrow');
  if length(butn_dwn),
   set(f,'windowbuttondownfcn','');
  end
 end

 FrequencyUnitStr = 'rad/sec';
 FrequencyUnit = 1;
 if infmat(9,1) == 1,
   if ~get(QFTToolData.Mouse(6),'value'),
    FrequencyUnitStr = 'hz';
    FrequencyUnit = 1/(2*pi);
   end
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
   set(txt3,'string',sprintf('Frequey: %0.3g %s',w(wloc)*FrequencyUnit, FrequencyUnitStr));
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
end
