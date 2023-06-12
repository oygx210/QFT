function tbxStruct = demos
% DEMOS Return demo information to the Matlap Expo.
% Copyright (c) 2002 by Terasoft, Inc.
% $Revision: 1.8 $ $Date: 1998/04/10 15:36:23 $

if nargout == 0, demo toolbox; return; end

tbxStruct.Name = 'QFT Feedback Control Design';
tbxStruct.Type = 'toolbox';
tbxStruct.Help = { ...
   ' The QFT Control Design Toolbox '
   ' is a collection of MATLAB functions for '
   ' designing robust feedback systems '
   ' using Quantitative Feedback Theory.  '
   ' '
   ' QFT is an engineering method devoted to '
   ' the practical design of feedback systems. '
   ' It uses frequecy domain concepts to satisfy '
   ' performance specifications and to handle '
   ' plant uncertainty.  '
   ' '
   ' The QFT Control Toolbox includes a convenient '
   ' GUI that facilitates loop shaping for the '
   ' formulation of simple, low-order controllers to '
   ' meet design requirements.'
   ' '
   ' Make sure to bring the command window forward '
   ' when you start the examples.'};
tbxStruct.DemoList = { ...
   'Main Example',                     'qftex1';
   'QFT Tracking',                     'qftex2';
   'Non-Parametric Uncertainty',       'qftex3';
   'Classical Design for Fixed Plant', 'qftex4';
   'ACC Benchmark',                    'qftex5';
   'Missile',                          'qftex6';
   'Cascaded: Inner-Outer',            'qftex7';
   'Cascaded: Outer-Inner',            'qftex8';
   'Uncertain Flexible Mechanism',     'qftex9';
   'Inverted Pendulum',                'qftex10';
   'Active Vibration Isolation',       'qftex11';
   'Main Example - Discrete-time'      'qftex12';
   'QFT Tracking - Discrete-time',     'qftex13';
   'CD Mechanism - Sampled-data',      'qftex14';
   '2x2 MIMO',                         'qftex15'};
