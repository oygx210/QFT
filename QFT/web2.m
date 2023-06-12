function web2(html_file)
% WEB2 Wrapper function for WEB.
%
%      WEB2 correctly modifies the file locator
%      string so that Netscape can find the
%      requested HTML file on the MAC.

% Author: Craig Borghesani
% Date: 11/13/96 10:46AM
% Copyright (c) 2003, Terasoft, Inc.
% $Revision: 1.3 $ $Date: 1997/05/15 18:02:51 $

if ismac,
   html_file = ['file:/',html_file];
end

web(html_file);
