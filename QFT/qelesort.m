function qelesort(sortMethod)

f = findobj(allchild(0),'tag','CAD Window');
sortByType = findobj(f,'tag','sortByType');
sortAll = findobj(f,'tag','sortAll');

switch sortMethod,
   case('ByType'),
      if strcmp(get(sortByType,'checked'),'on')
         set(sortByType,'checked','off');
         set(sortAll,'checked','on');

      else
         set(sortByType,'checked','on');
         set(sortAll,'checked','off');

      end

   case('All'),
      if strcmp(get(sortAll,'checked'),'on')
         set(sortByType,'checked','on');
         set(sortAll,'checked','off');

      else
         set(sortByType,'checked','off');
         set(sortAll,'checked','on');

      end

end %switch

% update listing
cont=get(bthan(19),'userdata');
if isempty(cont),
   cont = get(bthan(3),'userdata');
end
ElementsListbox = QFTToolData.Elements(17);
ListboxValue = get(ElementsListbox,'value');
ListboxInfo = get(ElementsListbox,'userdata');
ElementLocation = ListboxInfo(ListboxValue);
[ControllerString,ListboxInfo,ListboxValue] = cntstr(f,cont,ElementLocation);
set(ElementsListbox, 'string', ControllerString,...
                     'userdata', ListboxInfo,...
                     'value',ListboxValue);

