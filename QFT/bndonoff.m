function bndonoff(UIOperation)
% BNDONOFF Bound buttons. (Utility Function)
%          BNDONOFF 'turns on' and 'turns off' the bound which the user
%          has selected using the push buttons on the  side of the IDE
%          environment in LPSHAPE and DLPSHAPE.

% Author: Craig Borghesani
% 10/10/93
% Copyright (c) 2003, Terasoft, Inc.

f=gcf;
QFTToolData = getappdata(f,'QFTToolData');
BoundsListbox = QFTToolData.Bounds(3);
BoundsListboxInfo = get(BoundsListbox,'userdata');
BoundsListboxValue = get(BoundsListbox,'value');
BoundsListboxString = get(BoundsListbox,'string');

selectionType = get(f,'selectiontype');
currentObj = get(f,'currentobject');

switch UIOperation,
   case('On/Off'),
      if strcmp(selectionType,'open') | strcmp(get(currentObj,'string'),'On/Off'),
         for BoundIndex = BoundsListboxValue(:)',
            bnddata = BoundsListboxInfo{BoundIndex};
            bndstr = BoundsListboxString{BoundIndex};
            if strcmp(get(bnddata(1),'vis'),'on'),
               set(bnddata(bnddata~=0),'vis','off');
               bndstr = strrep(bndstr,'on','off');
            else
               set(bnddata(bnddata~=0),'vis','on');
               bndstr = strrep(bndstr,'off','on');
            end
            BoundsListboxString{BoundIndex} = bndstr;
         end
      end

   case('All On'),
      for BoundIndex = 1:length(BoundsListboxString),
         bnddata = BoundsListboxInfo{BoundIndex};
         set(bnddata(bnddata~=0),'vis','on');
         bndstr = BoundsListboxString{BoundIndex};
         if length(findstr(bndstr, 'off')),
            bndstr = strrep(bndstr,'off','on');
            BoundsListboxString{BoundIndex} = bndstr;
         end

      end

   case('All Off'),
      for BoundIndex = 1:length(BoundsListboxString),
         bnddata = BoundsListboxInfo{BoundIndex};
         set(bnddata(bnddata~=0),'vis','off');
         bndstr = BoundsListboxString{BoundIndex};
         if length(findstr(bndstr, 'on')),
            bndstr = strrep(bndstr,'on','off');
            BoundsListboxString{BoundIndex} = bndstr;
         end
      end

end % switch UIOperation
set(BoundsListbox,'string',BoundsListboxString);
