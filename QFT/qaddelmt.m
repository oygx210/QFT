function qaddelmt
% QADDELMT Set up for Add interface. (Utility Function)
%          QADDELMT sets up the ADD window for the various elements allowed
%          within the IDEs.

% Author: Craig Borghesani
% 8/31/93
% Copyright (c) 2003, Terasoft, Inc.
%       $Revision: 1.5 $

f = gcf;

bthan=get(f,'userdata');
infmat=get(bthan(16),'userdata');
hint_bar = get(bthan(36),'userdata');
PopupIndex = get(infmat(16,4),'value');
PopupString = get(infmat(16,4),'string');

lomat=get(bthan(1),'userdata');
cont=get(bthan(3),'userdata');

set(bthan(19),'userdata',cont);
set(bthan(20),'userdata',lomat);

set(hint_bar,'string','Enter values into the appropriate edit boxes');
set(infmat(16,1:3),'string','','enable','off');
set(infmat(19,1:3),'string','');
set(infmat(13,1),'enable','on','callback',['qelmts(',num2str(PopupIndex),',0)']);

switch PopupString{PopupIndex},
  case('Gain'),
  case('Real Pole'),
	 str=['pole';'zero'];
 set(infmat(16,1),'enable','on');
 set(infmat(19,1),'string','pole');

	case('Complex Pole', 'Complex Zero'}),
	case('Super 2nd'),
	case('Integrator/Differentiator'),
	case('Lead or Lag'),
	case('Notch'),
end; %switch


if any(PopupIndex==[1 2]),   % add pole or zero

elseif any(PopupIndex==[3 4]),  % add second order
 set(infmat(16,1:2),'enable','on');
 if PopupIndex==3,
  set(infmat(19,2),'string','zeta(pole)');
  set(infmat(19,1),'string','wn(pole)');
 else
  set(infmat(19,2),'string','zeta(zero)');
  set(infmat(19,1),'string','wn(zero)');
 end

elseif any(PopupIndex==[0.5,0.6,0.7]),  % add integrator/differentiator/delay/pred
 set(infmat(16,1),'enable','on');
 set(infmat(19,1),'string','n');

elseif PopupIndex==5,   % add lead/lag
 set(infmat(16,1:2),'enable','on');
 set(infmat(19,2),'string','phase');
 set(infmat(19,1),'string','w');

elseif PopupIndex==6,   % add notch
 set(infmat(16,1:3),'enable','on');
 set(infmat(19,3),'string','zeta1(zero)');
 set(infmat(19,2),'string','zeta2(pole)');
 set(infmat(19,1),'string','wn');

end
