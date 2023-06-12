function cntrecl
% CNTRECL Recall controller. (Utility Function)
%         CNTRECL either clears the present design and returns the user
%         to what the IDE function was originally called with or returns
%         the user to what was saved with Save (pull down menu).

% Author: Craig Borghesani
% 9/2/93
% Copyright (c) 2003, Terasoft, Inc.

f=gcf;
bthan=get(f,'userdata');
infmat=get(bthan(16),'userdata');
QFTToolData = getappdata(f,'QFTToolData');
ElementsListbox = QFTToolData.Elements(17);

hint_bar = get(bthan(36),'userdata');
old_cont=get(bthan(3),'userdata');
old_lomat=get(bthan(1),'userdata');
set(bthan(27),'userdata',old_lomat);
set(bthan(28),'userdata',old_cont);
set(infmat(8,1),'enable','on');

% close any existing windows
qclswin(0);

cont=get(bthan(5),'userdata');

if length(cont),  % return to saved controller
 initwl2=get(bthan(7),'userdata');
 initlo2=get(bthan(9),'userdata');
 lomat=ones(3+(infmat(9,1)==2),length(initwl2));
 lomat(1,:)=initwl2;
 lomat(2:2+(infmat(9,1)==2),:)=initlo2;
 set(hint_bar,'string','Recalling stored design');

else   % return to what CAD environment was called with
 cont=get(bthan(4),'userdata');
 initwl=get(bthan(6),'userdata');
 initlo=get(bthan(8),'userdata');
 lomat=ones(3+(infmat(9,1)==2),length(initwl));
 lomat(1,:)=initwl;
 lomat(2:2+(infmat(9,1)==2),:)=initlo;
 set(hint_bar,'string','Recalling initial design');

end

ElementLocation = 1;
[ControllerString,ListboxInfo,ListboxValue] = cntstr(f,cont,ElementLocation);
set(ElementsListbox, 'string', ControllerString,...
                     'userdata', ListboxInfo,...
                     'value',ListboxValue);
qelmtlistbox('Edit');

set(bthan(1),'userdata',lomat);
set(bthan(3),'userdata',cont);
if infmat(9,1)==1, qnicplt(f);
elseif infmat(9,1)==2, qmagplt(f);
elseif infmat(9,1)==3, mgphplot(f);
end
