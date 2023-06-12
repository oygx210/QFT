function qelmtpopup(UIOperation)
% QELMTPOPUP Set up element popup. (Utility Function)
%          QADDELMT sets up the popup menu for the various elements allowed
%          within the IDEs.

% Author: Craig Borghesani
% 8/31/93
% Copyright (c) 2003, Terasoft, Inc.
%       $Revision: 1.5 $

grey = get(0,'defaultuicontrolbackground');
ltgrey = [0.5,0.5,0.5]*1.5;
dkgrey = [0.5,0.5,0.5]*0.5;

f = gcf;

bthan=get(f,'userdata');
infmat=get(bthan(16),'userdata');
hint_bar = get(bthan(36),'userdata');
T=get(bthan(13),'userdata');
ElementEnvironment = infmat(9,1);
PopupIndex = get(infmat(16,4),'value');
PopupString = get(infmat(16,4),'string');
QFTToolData = getappdata(f,'QFTToolData');
ElementsListbox = QFTToolData.Elements(17);
ElementsListboxValue = get(ElementsListbox,'value');
ElementsListboxInfo = get(ElementsListbox,'userdata');
ElementLocation = ElementsListboxInfo(ElementsListboxValue);

cont=get(bthan(19),'userdata');
lomat=get(bthan(20),'userdata');
if isempty(cont),
   cont = get(bthan(3),'userdata');
   lomat=get(bthan(1),'userdata');
end

if infmat(9,1)==3,
   set(infmat(6,1:2),'vis','off','userdata',[]);
else
   set(infmat(6,1),'vis','off','userdata',[]);
end

set(bthan(31),'userdata',{[], []});

% reseting button-state functions
set(f,'windowbuttonupfcn','');
if length(get(f,'windowbuttonmotionfcn'))>1,
   set(f,'windowbuttonmotionfcn','qmouse(''Floating'')');
end

% reset complex and frequency values
infmat(4,3)=0;
infmat(3,3)=0;
set(bthan(16),'userdata',infmat);

set(hint_bar,'string','Enter values into the appropriate edit boxes');
set(infmat(16,1:3),'string','','enable','off', 'back', ltgrey, 'callback','');
set(infmat(17,1:3),'enable','off','callback','');
set(infmat(19,1:3),'string','');
set(infmat(13,1:4),'enable','off');

switch PopupString{PopupIndex},
   case('Gain'),
      set(ElementsListbox,'value',1);
      if any(diff(lomat(2,:))) | ElementEnvironment > 1,
         set(infmat(13,4),'callback','qscrelmt(1)','enable','on');
      end
      qelmtlistbox('Edit');

   case('Real Pole'),
      set(infmat(16,1),'enable','on','back','w');
      set(infmat(19,1),'string','pole');
      set(infmat(13,1),'enable','on',...
                       'callback','qelmts(1,0,''Add'')');
      if any(diff(lomat(2,:)) > 0.1) | ElementEnvironment > 1,
         set(infmat(13,4),'callback','qscrelmt(2)','enable','on');
      end
      if strcmp(UIOperation,'Done'),
         set(infmat(16,1),'string',num2str(cont(ElementLocation,1)));
      end

   case('Real Zero'),
      set(infmat(16,1),'enable','on','back','w');
      set(infmat(19,1),'string','zero');
      set(infmat(13,1),'enable','on',...
                       'callback','qelmts(2,0,''Add'')');
      if any(diff(lomat(2,:)) > 0.1) | ElementEnvironment > 1,
         set(infmat(13,4),'callback','qscrelmt(2)','enable','on');
      end
      if strcmp(UIOperation,'Done'),
         set(infmat(16,1),'string',num2str(cont(ElementLocation,1)));
      end

   case('Complex Pole'),
      set(infmat(16,1:2),'enable','on','back','w');
      set(infmat(19,1),'string','zeta');
      set(infmat(19,2),'string','wn');
      set(infmat(13,1),'enable','on','callback','qelmts(3,0,''Add'')');
      if any(diff(lomat(2,:)) > 0.1) | ElementEnvironment > 1,
         set(infmat(13,4),'callback','qscrelmt(3)','enable','on');
      end
      if strcmp(UIOperation,'Done'),
         set(infmat(16,1),'string',num2str(cont(ElementLocation,1)));
         set(infmat(16,2),'string',num2str(cont(ElementLocation,2)));
      end

   case('Complex Zero'),
      set(infmat(16,1:2),'enable','on','back','w');
      set(infmat(19,1),'string','zeta');
      set(infmat(19,2),'string','wn');
      set(infmat(13,1),'enable','on','callback','qelmts(4,0,''Add'')');
      if any(diff(lomat(2,:)) > 0.1) | ElementEnvironment > 1,
         set(infmat(13,4),'callback','qscrelmt(3)','enable','on');
      end
      if strcmp(UIOperation,'Done'),
         set(infmat(16,1),'string',num2str(cont(ElementLocation,1)));
         set(infmat(16,2),'string',num2str(cont(ElementLocation,2)));
      end

   case('Super 2nd'),
      if any(diff(lomat(2,:)) > 0.1) & ElementEnvironment == 1,
         set(infmat(13,4),'callback','qscrelmt(6)','enable','on');
      end

   case('Integrator/Differentiator'),
      if T == 0,
         set(infmat(16,1:2),'enable','on','back','w','string','0');
         set(infmat(19,1),'string','Integ');
         set(infmat(19,2),'string','Diff');
         if all(cont(2,1:2)==0)
            set(infmat(13,1),'enable','on','callback','qelmts(0.7,0,''Add'')');
         elseif ~strcmp(UIOperation,'Done'),
            set(ElementsListbox,'value',2);
            qelmtlistbox('Edit');
         end
         if strcmp(UIOperation,'Done'),
            set(infmat(16,1),'string',num2str(cont(ElementLocation,1)));
            set(infmat(16,2),'string',num2str(cont(ElementLocation,2)));
         end
      else
         set(infmat(16,1:2),'enable','on','back','w','string','0');
         set(infmat(19,1),'string','Integ');
         set(infmat(19,2),'string','Diff');
         if all(cont(2,1:2)==0)
            set(infmat(13,1),'enable','on','callback','qelmts(0.5,0,''Add'')');
         elseif ~strcmp(UIOperation,'Done'),
            set(ElementsListbox,'value',2);
            qelmtlistbox('Edit');
         end
         if strcmp(UIOperation,'Done'),
            set(infmat(16,1),'string',num2str(cont(ElementLocation,1)));
            set(infmat(16,2),'string',num2str(cont(ElementLocation,2)));
         end
      end

   case('Delay/Predictor'),
      set(infmat(16,1:2),'enable','on','back','w','string','0');
      set(infmat(19,1),'string','Delay');
      set(infmat(19,2),'string','Pred');
      if all(cont(3,1:2)==0)
         set(infmat(13,1),'enable','on','callback','qelmts(0.6,0,''Add'')');
      elseif ~strcmp(UIOperation,'Done'),
         set(ElementsListbox,'value',3);
         qelmtlistbox('Edit');
      end
      if strcmp(UIOperation,'Done'),
         set(infmat(16,1),'string',num2str(cont(ElementLocation,1)));
         set(infmat(16,2),'string',num2str(cont(ElementLocation,2)));
      end

   case('Lead or Lag'),
      set(infmat(16,1:2),'enable','on','back','w');
      set(infmat(19,1),'string','phase');
      set(infmat(19,2),'string','w');
      set(infmat(13,1),'enable','on','callback','qelmts(5,0,''Add'')');
      if any(diff(lomat(2,:)) > 0.1) & ElementEnvironment == 1,
         set(infmat(13,4),'callback','qscrelmt(4)','enable','on');
      end
      if strcmp(UIOperation,'Done'),
         set(infmat(16,1),'string',num2str(cont(ElementLocation,1)));
         set(infmat(16,2),'string',num2str(cont(ElementLocation,2)));
      end

   case('Notch'),
      set(infmat(16,1:3),'enable','on','back','w');
      set(infmat(19,1),'string','zeta1(zero)');
      set(infmat(19,2),'string','zeta2(pole)');
      set(infmat(19,3),'string','wn');
      set(infmat(13,1),'enable','on','callback','qelmts(6,0,''Add'')');
      if any(diff(lomat(2,:)) > 0.1) | ElementEnvironment > 1,
         set(infmat(13,4),'callback','qscrelmt(5)','enable','on');
      end
      if strcmp(UIOperation,'Done'),
         set(infmat(16,1),'string',num2str(cont(ElementLocation,1)));
         set(infmat(16,2),'string',num2str(cont(ElementLocation,2)));
         set(infmat(16,3),'string',num2str(cont(ElementLocation,3)));
      end

   case('Complex Lead'),
      set(infmat(16,1:3),'enable','on','back','w');
      set(infmat(19,1),'string','phase');
      set(infmat(19,2),'string','zeta');
      set(infmat(19,3),'string','w');
      set(infmat(13,1),'enable','on','callback','qelmts(7,0,''Add'')');
      if any(diff(lomat(2,:)) > 0.1) & ElementEnvironment == 1,
         set(infmat(13,4),'callback','qscrelmt(7)','enable','on');
      end
      if strcmp(UIOperation,'Done'),
         set(infmat(16,1),'string',num2str(cont(ElementLocation,1)));
         set(infmat(16,2),'string',num2str(cont(ElementLocation,2)));
         set(infmat(16,3),'string',num2str(cont(ElementLocation,3)));
      end

end %switch
