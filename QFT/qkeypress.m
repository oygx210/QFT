function qkeypress

f = gcf;
QFTToolData = getappdata(f,'QFTToolData');
DeleteButton = QFTToolData.Elements(18);
currentObj = get(f,'currentobject');

keyPressed = get(f,'currentcharacter');

if abs(keyPressed) == 127 & strcmp(get(DeleteButton,'enable'),'on'),
   qdelmts;

elseif strcmp(get(currentObj,'type'),'uicontrol') & abs(keyPressed) >= 48 & abs(keyPressed) <= 57,
   set(currentObj,'string',keyPressed);

end
