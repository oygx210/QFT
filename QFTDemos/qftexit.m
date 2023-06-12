function qftexit(next_stage,info_btn)
% QFTEXIT Exit from shaping environment from within Demo facility.
%         QFTEXIT setups the EXIT within the CAD environment to call the
%         file specified in STRING

% Author: Craig Borghesani
% 10/20/94
% Copyright (c) 2002 by Terasoft, Inc.
%       $Revision: 1.6 $

% setup Exit from within shaping environment to move to next stage
f=findobj(allchild(0),'tag','CAD Window'); xitstr = [];
UIMenuHandle = findobj(f,'type','uimenu','label','QFile','parent',f);
UIMenuExit = findobj(UIMenuHandle,'label','Exit');
f_str = int2str(f);
xitstr=['xitcade(',f_str,',3);',next_stage];
set(UIMenuExit,'callback',['close(',f_str,')']);
set(f,'closerequestfcn',xitstr);

%ui_han=get(f,'children');
%for k=1:length(ui_han),
% if strcmp(get(ui_han(k),'type'),'uimenu'),
%  if strcmp(get(ui_han(k),'label'),'QFile'),
%   ui_han2=get(ui_han(k),'children');
%   for k2=1:length(ui_han2),
%    if strcmp(get(ui_han2(k2),'label'),'Exit'),
%     f_str = int2str(f);
%     xitstr=['xitcade(',f_str,',3);',next_stage];
%     set(ui_han2(k2),'callback',['close(',f_str,',''force'')']);
%     set(f,'closerequestfcn',xitstr);
%    end
%   end
%  end
% end
%end

% setup continue button to exit shaping environment
close_str = ['close(',f_str,')'];
set(info_btn(3),'enable','on',...
    'callback',['set(gco,''enable'',''off'');set(findobj(allchild(0),''tag'',''QuitButton''),''enable'',''off'');eval(''',close_str,''',''qfterror(2)'');set(findobj(allchild(0),''tag'',''QuitButton''),''enable'',''on'')']);
set(info_btn(2),'enable','off');
