function qstabeta(flag)
%
% Utility function: QSTABETA
%
% The purpose of this function is to compute and display the
% stabilizing beta.

% Author: Craig Borghesani <cborg@terasoft.com>
% Date: 10/20/98 10:19PM
% Copyright (c) 2003, Terasoft, Inc.

if flag==0,
% qclswin(0);

 f=gcf;
 bthan=get(f,'userdata');
 infmat=get(bthan(16),'userdata');
 T=get(bthan(13),'userdata');
 proc_str=[];
 if infmat(25,2)>1, proc_str=['(',int2str(infmat(25,2)),')']; end

 fig_color=[128/255,128/255,128/255];
 proc_num = int2str(infmat(25,2));
 win_tag = findobj('tag',['qft10',proc_num]);
 if ~length(win_tag),
  infmat(25,4)=figure('name',['Stabilizing Beta ',proc_str],'numbertitle','off',...
                      'position',[100,100,220,60],...
                      'menubar','none','color',fig_color,...
                      'resize','off',...
                      'vis','off','userdata',f,'tag',['qft10',proc_num],...
                      'closerequestfcn','set(gcf,''vis'',''off'')',...
                      'handlevisibility','callback');
  set(gca,'vis','off');
  txt(1)=uicontrol('style','text','pos',[10,30,140,17],...
      'string','Stabilizing Beta:','horiz','left');

  edt(1)=uicontrol('style','edit','pos',[155,30,60,20],...
                          'background',[1,1,1],'string','1');
  edt(2)=uicontrol('style','edit','pos',[10,5,80,20],...
                          'background',[1,1,1]);
  edt(3)=uicontrol('style','edit','pos',[95,5,80,20],...
                          'background',[1,1,1]);

  set(bthan(16),'userdata',infmat);
  set(bthan(39),'userdata',edt);
  drawnow;
  set(infmat(25,4),'vis','on');
 else
  set(infmat(25,4),'vis','on');
 end
end
