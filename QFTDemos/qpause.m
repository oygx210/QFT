function qpause
% QPAUSE Used instead of PAUSE to allow graph window to remain active
%        inside an M file.  QPAUSE also makes sure to prompt user to
%        exit shaping environment before proceeding

% Author: Craig Borghesani
% 8/23/94
% Copyright (c) 2002 by Terasoft, Inc.
%       $Revision: 1.4 $

x = 1;
datal = 37;
while datal == 37,

% special pause that allows use of handle graphics while in a script file
 while ~isempty(x),
  x = input(' ');
 end

 close(findobj(allchild(0),'tag','CAD Window'),'force');
 datal = 2;

% software 'catch' to make sure that user has exited from shaping environment
% chil = get(0,'children');
% if length(chil),
%  for k = 1:length(chil),
%   data = get(chil(k),'userdata');
%   datal = length(data);
%   if datal == 37,
%    disp('Please exit from the shaping environment first');
%    x = 1;
%    break;
%   end
%  end
% else
%  datal = 0;
% end

end
