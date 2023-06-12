function qzoomplt
% QZOOMPLT Zoom tool for the QFT Toolbox. (Utility Function)
%          QZOOMPLT handles the zoom option offered for PLOTTMPL, PLOTBNDS,
%          CHKSISO, and CHKGEN.

% Author: Craig Borghesani
% Date: 10/10/93
% Revised: 2/22/96 1:45 PM V1.1 updates
% Copyright (c) 2003, Terasoft, Inc.


f = gcf;
a = gca;
set(f,'units','pixels');

old_axis = get(f,'userdata');
xlim = get(a,'xlim');
ylim = get(a,'ylim');
pt0 = get(a,'currentpoint');
mouse_pos = get(f,'currentpoint');

if any(pt0(1,1:2) < [xlim(1),ylim(1)]) | ...
   any(pt0(1,1:2) > [xlim(2),ylim(2)]),
 set(a,'xlim',old_axis(1:2),'ylim',old_axis(3:4));
else
 rbbox([mouse_pos 0 0],[mouse_pos]);
 drawnow;
 pt1 = get(a,'currentpoint');
 minx = min(pt1(1,1),pt0(1,1));
 miny = min(pt1(1,2),pt0(1,2));
 maxx = max(pt1(1,1),pt0(1,1));
 maxy = max(pt1(1,2),pt0(1,2));
 if minx < maxx & miny < maxy,
  set(a,'xlim',[minx,maxx],'ylim',[miny,maxy]);
 end
end
set(f,'units','norm');
