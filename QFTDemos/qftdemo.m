function qftdemo(flag)
% QFTDEMO Demonstration facility.
%         QFTDEMO provides a GUI for the demonstration files included with
%         the toolbox.
%
%         Each demonstration has a matching example file QFTEX# where #
%         ranges from 1-15.

% Author: Craig Borghesani
% 10/6/93
% Copyright (c) 2002 by Terasoft, Inc.
%       $Revision: 1.7 $

% updated to run the matlab demo facility
demo

%if nargin==0,
% title = 'QFT Control Design Demos';
% figcolor = [128/255,128/255,128/255];

% f = colordef('new','none');
% set(f,'name',title,'numbertitle','off','pos',[10,100,270,382],...
%            'resize','off','color',figcolor,'vis','off');

% DemoList = { ...
%   'Main Example',                     'qft_val=1;qftdm1';
%   'QFT Tracking',                     'qft_val=2;qftdm2';
%   'Non-Parametric Uncertainty',       'qft_val=3;qftdm3';
%   'Classical Design for Fixed Plant', 'qft_val=4;qftdm4';
%   'ACC Benchmark',                    'qft_val=5;qftdm5';
%   'Missile',                          'qft_val=6;qftdm6';
%   'Cascaded: Inner-Outer',            'qft_val=7;qftdm7';
%   'Cascaded: Outer-Inner',            'qft_val=8;qftdm8';
%   'Uncertain Flexible Mechanism',     'qft_val=9;qftdm9';
%   'Inverted Pendulum',                'qft_val=10;qftdm10';
%   'Active Vibration Isolation',       'qft_val=11;qftdm11';
%   'Main Example (Discrete-time)'      'qft_val=12;qftdm12';
%   'QFT Tracking (Discrete-time)',     'qft_val=13;qftdm13';
%   'CD Mechanism (Sampled-data)',      'qft_val=14;qftdm14';
%   '2x2 MIMO',                         'qft_val=15;qftdm15'};

% uicontrol(f,'style','push','string','Close','pos',[10,5,250,20],...
%           'horizontalalignment','center',...
%           'callback','close(gcf)');
%%'clear;close(gcf);if exist(''uservals.mat''),load uservals;delete uservals.mat;end;');

% for k=15:-1:1,
%  str=DemoList{k,1};
%  h=uicontrol(f,'style','push','pos',[10,28+(15-k)*23,250,20],...
%              'callback',['eval(''',DemoList{k,2},''',''qfterror(2)'')'],...
%              'userdata',k);
%  set(h,'string',str,'horizontalalignment','center');
% end

% set(f,'visible','on');
%end
