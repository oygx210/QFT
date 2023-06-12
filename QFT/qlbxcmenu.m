function qlbxcmenu(UIOperation)
%QLBXCMENU Listbox context menu. (Utility Function)

% Author: Craig Borghesani <cborg@terasoft.com>
% Date: 2/17/02 8:36PM
% Copyright (c) 2002, Terasoft, Inc.

ElementListbox = gco;
ElementListboxString = get(ElementListbox,'string');
ElementListboxUIContextMenu = get(ElementListbox,'uicontextmenu');
SelectAllOption = findobj(ElementListboxUIContextMenu,'label','Select All');

switch UIOperation,
   case('Check'),
      if length(ElementListboxString) > 1,
         set(SelectAllOption,'enable','on');
      else
         set(SelectAllOption,'enable','off');
      end

   case('Select All'),
      set(ElementListbox,'value',2:length(ElementListboxString));
      qelmtlistbox('Edit');

end %switch

