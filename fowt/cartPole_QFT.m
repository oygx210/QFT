% Linearized inverted pendulum QFT control
%
%   AUTHOR  : Mohammad Odeh
%   DATE    : Jun.  5th, 2023
%
% CHANGELOG :
%   Jun.  5th, 2023
%       - Initial script
%
%   Jun. 19th, 2023
%       - Working controller is MGS_linearizedInvertedPendulum_Cart_V2.shp
%

%% Setup environment
clear all; close all; clc;
set( groot, 'defaultLineLineWidth', 1.5 );	% Set default line width of plots
set( 0, 'DefaultAxesFontSize', 12 );        % Set a default axes font size
% Change default interpreter (affects text, labels, etc...)
set( 0, 'DefaultTextInterpreter', 'latex' );
set( 0, 'DefaultLegendInterpreter', 'latex' );
set( 0, 'DefaultAxesTickLabelInterpreter', 'latex' );
format compact;
fontsize = 12;

% Flags
CNTR = 1;                                   % Figure handle counter
PLOT = true;                                %#ok<NASGU> If true, plot figures!
% PLOT = false;                               % COMMENT OUT TO PRINT FIGURES
PRNT = true;                                %#ok<NASGU>
% PRNT = false;                               % COMMENT OUT TO PRINT FIGURES

%% Add folders/files to path
% Get current path to working directory and split
pathParts = strsplit( pwd, filesep );
% Go up one level and generate new path
src = fullfile( pathParts{1:end-1} );

% If on a UNIX machine (i.e. macOS, Ubuntu, etc...), fix path since
% strsplit() removes the leading '/' from the path.
if( isunix )
    src = [ filesep src ];
end

% Add QFT2 to path
addpath( genpath(src) );

%% Generate SS model using analytical approach

M_0 = 2.0   ;                   % Mass of cart                  [  kg  ]
m_0 = 0*0.075 ;                   % Mass of rod                   [  kg  ]
g   = 9.8   ;                   % Gravitational acceleration    [m.s^-2]
h   = 0.5   ;                   % Rod length                    [  m   ]

% ------------------------------
A_theory = [ 0      1       0                   0   ;
             0      0   -m_0*g/M_0              0   ;
             0      0       0                   1   ;
             0      0   (M_0+m_0)*g/(M_0*h)     0 ] ;
% ------------------------------
B_theory = [ 0              ;
             1/M_0          ;
             0              ;
             -1/(M_0*h) ]   ;
% ------------------------------
C_theory = [ 1  0   0   0   ;
             0  0   1   0 ] ;
% ------------------------------
D_theory = [ 0  ;
             0 ];

stateNames = {'x' 'x_dot' 'phi' 'phi_dot'};
inputNames = {'u'};
outputNames = {'x'; 'phi'};

sys_theory = ss( A_theory, B_theory, C_theory, D_theory , ...
                 'StateName'    ,   stateNames          , ...
                 'InputName'    ,   inputNames          , ...
                 'OutputName'   ,   outputNames );

% --- Generate TF from SS model
TF_theory = tf( sys_theory );

%% Cart-Pole tranfer functions
p_11 = tf( [-h, 0, g], [-M_0*h, 0, (M_0+m_0)*g, 0] );   % Cart TF
p_21 = tf( 1, [-M_0*h, 0, (M_0+m_0)*g] );               % Pole TF

%% Get feedback controller, G(s), and pre-filter, F(s), generated using QFT

% [INFO] ...
fprintf( 'Retrieving stored G(s) and F(s)...' );

% --- Directory where QFT generated controllers are stored
src = './controllerDesigns/';

% --- Pole controller, G_theta(s)
% G_file  = [ src 'linearInvPend_Pole_Simplified_V2.shp' ];
G_file  = [ src 'linearInvPend_Pole_V2.shp' ];
G_theta = getqft( G_file );

% --- Cart controller, G_x(s)
G_file  = [ src 'linearInvPend_Cart_Simplified_V2.shp' ];
G_x     = getqft( G_file );

% --- Cart controller, G_x(s)
F_file  = [ src 'linearInvPend_Cart_Simplified_V2.fsh' ];
F_x     = getqft( F_file );

fprintf( 'DONE!\n' );

%% Misc. Operations for verification purposes

syms s;
I = eye(size(A_theory));
Ps = C_theory*(s*I - A_theory)^-1*B_theory;
clear s;


