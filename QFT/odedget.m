function outCont = odedget
% ODEDGET Pull out controller matrix from design environment.

% Author: Craig Borghesani <cborg@terasoft.com>
% Date: 1/19/03 11:00AM

loopfig = findobj(allchild(0),'tag','CAD Window');
bthan=get(loopfig,'userdata');

cont = get(bthan(3),'userdata');


