function qfterror(flag,info)
% QFTERROR Caught error function. (Utility Function)
%          QFTERROR is executed if there is an error with SISOBNDS and GENBNDS.
%          It is also used by the demo facility to catch errors that may
%          occur.

% Author: Craig Borghesani
% 7/18/94
% Copyright (c) 2003, Terasoft, Inc.

if flag == 1, % An error occurred in a bound manager function

 close(info(1));

elseif flag == 2, % An error occurred during a demo session

 h = errordlg('An error has occurred in this demo','message','on');
 butn = findobj(h,'style','pushbutton');
 callbk = get(butn,'callback');
 set(butn,'callback',[callbk,';close([info_win]);']);

end
