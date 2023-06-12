function qdelmts
% QDELMTS Compute and delete individual terms. (Utility Function)
%        QDELMTS computes the various frequency responses of the element
%        that is selected from within the Add or Edit window.

% Author: Craig Borghesani
% Date: 1/26/02 2:02PM
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
if isempty(get(bthan(19),'userdata')),
	cont = get(bthan(3),'userdata');
	lomat=get(bthan(1),'userdata');
end
T=get(bthan(13),'userdata');
selctn=get(bthan(30),'userdata');
ElementsListbox = QFTToolData.Elements(17);

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


ListboxString = get(ElementsListbox,'string');
ListboxValue = get(ElementsListbox,'value');
ListboxInfo = get(ElementsListbox,'userdata');

% setting up for whether FSHAPE/DFSHAPE is being used
if any(ElementEnvironment==[1 3]), % SHAPE/DSHAPE/BODPLOT/DBODPLOT
 q=1; s=0;
elseif ElementEnvironment==2, % FSHAPE/DFSHAPE
 q=[1;1]; s=1;
end

del_cp = ones(size(lomat(1,:)));

SelectedElements = ListboxInfo(ListboxValue);
RemoveElements = [];
for SelectedElement = SelectedElements(:)',
   del_ele = cont(SelectedElement, :);
   del_cp = qcntbode(del_ele,lomat(1,:),T) .* del_cp;

   if cont(SelectedElement, 4) == 0.5,
      cont(SelectedElement,1:2) = 0;

   elseif any(cont(SelectedElement, 4) == [0.6, 0.7]),
      cont(SelectedElement,1:2) = 0;

   else
      RemoveElements = [RemoveElements, SelectedElement];

   end % if cont(SelectedElement, 4) == 0.5
end % for SelectedElement = SelectedElements(:)'
cont(RemoveElements, :) = [];
lomat(2:2+s,:)=lomat(2:2+s,:)./del_cp(q,:);
ElementLocation = size(cont,1);

if isempty(get(bthan(19),'userdata')),
   set(bthan(3),'userdata',cont);
   set(bthan(1),'userdata',lomat);

else
   set(bthan(19),'userdata',cont);
   set(bthan(20),'userdata',lomat);
end

if ElementEnvironment==1, qnicplt(f);
elseif ElementEnvironment==2, qmagplt(f);
elseif ElementEnvironment==3, mgphplot(f);
end
[ControllerString,ListboxInfo,ListboxValue] = cntstr(f,cont,ElementLocation);
set(ElementsListbox, 'string', ControllerString,...
                     'userdata', ListboxInfo,...
                     'value',ListboxValue);

qelmtlistbox('Delete');
