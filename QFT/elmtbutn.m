function elmtbutn(flag,f,f2)
% ELMTBUTN Element buttons. (Utility Function)
%          ELMTBUTN creates the uicontrols for the Add window.  The
%          selection of buttons varys depending upon which IDE and whether
%          it is discrete or continuous

% Author: Craig Borghesani
% 9/5/93
% Copyright (c) 2003, Terasoft, Inc.
%       $Revision: 1.5 $

bthan=get(f,'userdata');
T=get(bthan(13),'userdata');

% determine which buttons to display
if any(flag==[1 3]) & T > 0,
 vec=[1 2 3 4 5 6 7 8];
 ele=[1,2,3,4,0.6,0.5,5,6];
elseif any(flag==[1 3]),
 vec=[1 2 3 4 6 7 8];
 ele=[1,2,3,4,0,0.7,5,6];
else
 vec=[1 2 3 4 6 8];
 ele=[1,2,3,4,0,0.5*(T > 0)+0.7*(T == 0),0,6];
end

top=20*(length(vec)-1)+65; dely=20;
str=str2mat('Real Pole','Real Zero','Complex Pole','Complex Zero',...
            'Delay/Pred','Integ/Diff','Lead/Lag','Notch');
len=[9,9,12,12,10,10,8,5];

for k=vec,
 opr=['qaddelmt(',num2str(ele(k)),')'];
 uicontrol(f2,'style','push','pos',[35 top 160 dely],...
           'horizontalalignment','center','string',str(k,1:len(k)),...
           'callback',opr);
 top=top-dely;
end
