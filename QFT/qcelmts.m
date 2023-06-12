function qcelmts
% QCELMTS Compute and convert lead/lag terms. (Utility Function)
%        QCELMTS computes the various frequency responses of the elements
%        that is selected from within the elements listbox.

% Author: Craig Borghesani
% Date: 1/26/02 2:02PM
% Copyright (c) 2003, Terasoft, Inc.
%       $Revision: 1.5 $

f=gcf;
bthan=get(f,'userdata');
infmat=get(bthan(16),'userdata');
QFTToolData = getappdata(f,'QFTToolData');

ElementEnvironment=infmat(9,1);
cont=get(bthan(19),'userdata');
if isempty(cont),
	cont = get(bthan(3),'userdata');
   set(bthan(19),'userdata',cont);
end
T=get(bthan(13),'userdata');
ElementsListbox = QFTToolData.Elements(17);

ListboxString = get(ElementsListbox,'string');
ListboxValue = get(ElementsListbox,'value');
ListboxInfo = get(ElementsListbox,'userdata');

SelectedElements = ListboxInfo(ListboxValue);

contzp=cont(SelectedElements,:);
cont(SelectedElements,:)=[];
[contldlg,msg]=zp2ldlg(contzp,T);
if length(contldlg),
   contnew=[cont;contldlg];
   set(bthan(19),'userdata',contnew);

elseif msg==3,
   errordlg('Unstable ZERO selected','Message','on');
   return;

elseif msg==4,
   errordlg('Unstable POLE selected','Message','on');
   return;

end

ElementLocation = size(contnew,1);
[ControllerString,ListboxInfo,ListboxValue] = cntstr(f,contnew,ElementLocation);
set(ElementsListbox, 'string', ControllerString,...
                     'userdata', ListboxInfo,...
                     'value',ListboxValue);
qelmtlistbox('Convert');
