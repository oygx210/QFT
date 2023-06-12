function qbtnkill
% QBTNKILL Remove button down functions. (Utility Function)
%          QBTNKILL resets the button down function of all line objects to
%          the null setting ''.

% Author: Craig Borghesani
% 5/12/93
% Copyright (c) 2003, Terasoft, Inc.



f2=gcf;
f=get(f2,'userdata');
bthan=get(f,'userdata');
a=gca;
han = get(bthan(26),'userdata');

sel_mat = get(han(6),'userdata');
sz=size(sel_mat);

set(sel_mat(2:3:sz(1),2),'buttondownfcn','');
