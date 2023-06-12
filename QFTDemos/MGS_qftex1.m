% QFTEX1 Main Example.

% Author: Y. Chait
% 10/10/94
% 4.18.03 revised for V2 (YC)
% Copyright (c) 2002 by Terasoft, Inc.

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
PRNT = true;                                %#ok<NASGU>
% PRNT = false;                               % COMMENT OUT TO PRINT FIGURES

% Plot line color/style
c_line = [ 'r', 'g', 'b', 'c', 'm', 'k' ];
C_line = [ c_line, c_line, c_line, c_line ];
C_line = [ C_line, C_line, C_line, C_line ];

% Get current path to working directory and split
pathParts = strsplit( pwd, filesep );
% Go up one level and generate new path
src = fullfile( pathParts{1:end-1} );
% Add QFT2 to path
addpath( genpath(src) );

%% [INFO] Strings
ACK = 'COMPLETED\n\n';

%%
% Example 2.1 from Mario Garcia-Sanz's:
%   "Robust Control Engineering: Practical QFT Solutions"
%
%
%	Consider a continuous-time, siso, negative unity feedback system
%
%                                        |V(s)                |D(s)
%          ----  +           ----        v        ----        v
%   ----> |F(s)|--> O ----> |G(s)| ----> O ----> |P(s)| ----> O ------>
%   R(s)   ----     ^ -      ----                 ----        |  Y(s)
%                   |                  ----                   |
%                   ----------------- |H(s)| -----------------
%                                      ----
%
%
%       The plant P(s) has parametric a uncertainty model:
%	                                 KT                    
%	   P(s)  =  ------------------------------------------- 
%	            (J.La)s^2 + (J.Ra + b.La)s + (b.Ra + KT.Ke)     
%
%   Where
%       J  in [ J_min,  J_max], b  in [ b_min,  b_max]
%       La in [La_min, La_max], Ra in [Ra_min, Ra_max]
%       KT in [KT_min, KT_max],         Ke = KT
%
%
%	The performance specifications are: design a controller
%       G(s) such that it achieves
%
%  1) Type 1: Stability specification
%       ==> Corresponds to sisobnds( 1, ... )
%
%                   |  P(jw)G(jw)  |
%       |T_1(jw)| = |--------------| <= del_1(w) == W_s == 1.46
%                   |1 + P(jw)G(jw)|
%
%       w in [1, 10000] rad/s
%
%  2) Type 3: Sensitivity OR Disturbance at plant output specification 
%       ==> Corresponds to sisobnds( 2, ... )
%
%                   |       1      |                  (s/a_d)
%       |T_3(jw)| = |--------------| <= del_3(w) == -----------
%                   |1 + P(jw)G(jw)|                (s/a_d) + 1
%
%       w in [1, 1000] rad/s ; a_d = 1,000
%
%  3) Type 6: Reference tracking specification
%       ==> Corresponds to sisobnds( 7, ... )
%
%                                 |        P(jw)G(jw)  |
%       del_6_lo(w) < |T_3(jw)| = |F(jw) --------------| <= del_6_hi(w)
%                                 |      1 + P(jw)G(jw)|
%
%       w in [1, 5000] rad/s
%
%   Where
%
%                     |        1        |
%       del_6_lo(w) = |-----------------| ; a_L = 1000
%                     | ((s/a_L) + 1)^2 |
%
%                     |       ((s/a_U) + 1)       |
%       del_6_hi(w) = |---------------------------| ; a_U = 1000
%                     | (s/w_n)^2 + (2zs/w_n) + 1 |
%
%       z = 0.8 , w_n = 1.25a_u/z
%
%

%% Step 1: Plant Modeling & Uncertainty

% --- Plant parameters
%   min_    : Minimum value
%   max_    : Maximum value
%   grid_   : Gridding
%
min_J   = 4e-5; max_J   = 8e-5; grid_J  = 5;
min_b   = 6e-5; max_b   = 9e-5; grid_b  = 3;
min_La  = 1e-3; max_La  = 3e-3; grid_La = 5;
min_Ra  = 0.5 ; max_Ra  = 0.7 ; grid_Ra = 4;
min_KT  = 0.04; max_KT  = 0.06; grid_KT = 3;


% --- Gridding
%   ***NOTE: Can grid using logspace() or linspace()
%   _g  : Gridded variable
%
J_g = logspace( log10(min_J)    , log10(max_J)  , grid_J  );
b_g = logspace( log10(min_b)    , log10(max_b)  , grid_b  );
La_g= logspace( log10(min_La)   , log10(max_La) , grid_La );
Ra_g= logspace( log10(min_Ra)   , log10(max_Ra) , grid_Ra );
KT_g= logspace( log10(min_KT)   , log10(max_KT) , grid_KT );


% --- Plant generation
%   *** Note on transfer function generation:
%       The first two indices represent the number of outputs and
%       inputs for the models, while the third index is the number
%       of models in the array.
%
%       i.e. => P( 1, 1, 900 ) == SISO with 900 TFs
%
n_Plants = grid_J*grid_b*grid_La*grid_Ra*grid_KT;   % Number of plants
P = tf( zeros(1,1,n_Plants) );                      % Pre-allocate memory

% [INFO] ...
fprintf( 'Step 1:' );
fprintf( '\tComputing QFT templates using %3i plants...', n_Plants );

NDX = 1;                                            % Plant counter
for var1 = 1:grid_J                                 % Loop over J
    J = J_g( var1 );                                % ....

    for var2 = 1:grid_b                             % Loop over b
        b = b_g( var2 );                            % ....

        for var3 = 1:grid_La                        % Loop over La
            La = La_g( var3 );                      % ...

            for var4 = 1:grid_Ra                    % Loop over Ra
                Ra = Ra_g( var4 );                  % ...

                for var5 = 1:grid_KT                % Loop over KT
                    KT = KT_g( var5 );              % ...

                    % --- Here we create the plant TF
                    num = KT;
                    den = [ J*La (J*Ra + b*La) (b*Ra + KT*KT) ];
                    P(1, 1, NDX) = tf( num, den );  % Transfer Function
                    NDX = NDX + 1;                  % Incerement counter
                end
            end
        end
    end
end

% [INFO] ...
fprintf( ACK );

%% Step 2: The Nominal Plant
J_0     = 6e-5;
b_0     = 6e-5;
La_0    = 1e-3;
Ra_0    = 0.7 ;
KT_0    = 0.04;

% [INFO] ...
fprintf( 'Step 2:' );
fprintf( '\tComputing nominal plant...' );

% --- Generate nominal plant TF
%   Any one of the models above can be used as the nominal plant.
%   We just happened to chose this one.
%
num = KT_0;
den = [ J_0*La_0 (J_0*Ra_0 + b_0*La_0) (b_0*Ra_0 + KT_0*KT_0) ];
P_0(1, 1, 1) = tf( num, den );          % Nominal Transfer Function

% --- Append to the end of the gridded plants
P( 1, 1, end+1 ) = P_0;
% --- Define nominal plant case
nompt = length( P );

% [INFO] ...
fprintf( ACK );

% --- Plot bode diagram
w = [ 1 5 10 50 100 200 300 500 700 1000 2000 5000 10000 ];
figure( CNTR ); CNTR = CNTR + 1;
bode( P_0, w ); grid on;
[p0, theta0] = bode( P_0, w );

make_nice_plot();

%% Step 3: QFT Template

% [INFO] ...
fprintf( 'Step 3:' );
fprintf( '\tPlotting QFT templates...' );

% --- Working frequencies
w = [ 1 5 10 50 100 200 300 500 700 1000 2000 5000 10000 ]; 

% --- Plot QFT templates
plottmpl( w, P, nompt );

% --- Change legend position
hLegend = findobj(gcf, 'Type', 'Legend');   % Get legend property
set( hLegend, 'location', 'southeast' );    % Access and change location
xlim( [-200 10] );
title('Plant Templates')

% --- Beautify plot
make_nice_plot();

% [INFO] ...
fprintf( ACK );

%% Step 4: Define Stability Specifications
% --- The stability margins, gain and phase, are specified here
%

% L(s) = G(s)P(s) where G(s) == control, P(s) == plant
%
%                L(s)        G(s)P(s)
% Then H(s) = --------- = --------------
%              1 + L(s)    1 + G(s)P(s)  
%

% [INFO] ...
fprintf( 'Step 4:' );
fprintf( '\tDefining stability specifications...' );

% --- Type 1
% Frequencies of interest
omega_1 = [ 1 5 10 50 100 200 300 500 700 1000 2000 5000 10000 ];
% Restriction
W_s         = 1.46;
del_1       = W_s;
PM          = 180 -2*(180/pi)*acos(0.5/W_s);         % In deg
GM          = 20*log10( 1+1/W_s );                   % In dB

% [INFO] ...
fprintf( ACK );

%% Step 5: Define Performance Specifications

% [INFO] ...
fprintf( 'Step 5:' );
fprintf( '\tDefining performance specifications...' );

% --- Type 3
% Frequencies of interest
omega_3 = [ 1 5 10 50 100 200 300 500 700 1000 ];

% Restriction
a_d     = 1e3;
num     = [ 1/a_d 0 ];
den     = [ 1/a_d 1 ];
del_3   = tf( num, den );


% --- Type 6
% Frequencies of interest
omega_6 =[ 1 5 10 50 100 200 300 500 700 1000 2000 5000 ];

% Restriction
% Upper bound
a_U = 1000; zeta = 0.8; wn = 1.25*a_U/zeta;
num = [ 1/a_U 1 ];
den = [ (1/wn)*(1/wn) (2*zeta/wn) 1 ];
del_6_hi = tf( num, den );
% Lower bound
a_L = 1000;
num = 1;
den = [ 1e-6 0.002 1 ];
del_6_lo = tf( num, den );
% Tracking weight
del_6 = [ del_6_hi  ;
          del_6_lo ];

% [INFO] ...
fprintf( ACK );

%% Step 6: Calculate Staibility QFT Bounds

% --- Example 2.1 continued (Pg. 36)
%   - Type 1: Stability specification
%       > Corresponds to sisobnds( 1, ... )
%   - Type 3    //
%       > Corresponds to sisobnds( 2, ... )
%   - Type 6    //
%       > Corresponds to sisobnds( 7, ... )
%

% --------------------------------------------------
% ----      Type 1: Stability specification     ----
% --------------------------------------------------
spec = 1;

% [INFO] ...
fprintf( 'Step 6:' );
fprintf( '\tCalculating stability QFT bounds\n' );
fprintf( '\tComputing bounds: ' );
fprintf( 'bdb%i = sisobnds( %i, ... )\n', spec, spec );
fprintf( '\t\t > ' );

% --- Compute bounds
bdb1 = sisobnds( spec, omega_1, del_1, P, [], nompt );
% R = 0; bdb1 = sisobnds( spec, omega_1, del_1, P, R, nompt );

% [INFO] ...
fprintf( 'Plotting bounds...' );

% --- Plot bounds
plotbnds( bdb1 );
title( 'Robust Stability Bounds' );
xlim( [-360 0] ); ylim( [-15 25] );
make_nice_plot();

% [INFO] ...
fprintf( ACK );

%% Step 7: Calculate Performance QFT Bounds

% -------------------------------------------
% ---- Type 3: Sensitivity specification ----
% -------------------------------------------
spec = 2;

% [INFO] ...
fprintf( '\tComputing bounds: ' );
fprintf( 'bdb%i = sisobnds( %i, ... )\n', spec, spec );
fprintf( '\t\t > ' );

% --- Compute bounds
bdb2 = sisobnds( spec, omega_3, del_3, P, [], nompt );

% [INFO] ...
fprintf( 'Plotting bounds...' );

% --- Plot bounds
plotbnds(bdb2);
title('Sensitivity Reduction Bounds');
make_nice_plot();

% [INFO] ...
fprintf( ACK );

% --------------------------------------------------
% ---- Type 6: Reference tracking specification ----
% --------------------------------------------------
spec = 7;

% [INFO] ...
fprintf( '\tComputing bounds: ' );
fprintf( 'bdb%i = sisobnds( %i, ... )\n', spec, spec );
fprintf( '\t\t > ' );

% --- Compute bounds
bdb7 = sisobnds( spec, omega_6, del_6, P );

% [INFO] ...
fprintf( 'Plotting bounds...' );

% --- Plot bounds
plotbnds(bdb7);
title('Robust Tracking Bounds');
make_nice_plot();

% [INFO] ...
fprintf( ACK );


%% Step 8: Intersection of QFT Bounds and Compatibility

% [INFO] ...
fprintf( 'Step 8:' );
fprintf( '\tGrouping bounds...' );

% --- Grouping bounds
bdb = grpbnds( bdb1, bdb2, bdb7 );
plotbnds(bdb);
title('All Bounds');

% [INFO] ...
fprintf( ACK );
fprintf( '\tIntersection of bounds...' );

% --- Find bound intersections
ubdb=sectbnds(bdb);
plotbnds(ubdb);
title('Intersection of Bounds');

% [INFO] ...
fprintf( ACK );

%% Step 9: Synthesize Feedback Controller G(s)

% [INFO] ...
fprintf( 'Step 9:' );
fprintf( '\tSynthesize G(s)...' );

% --- Design TF of G(s)
syms s;
num = sym2poly( 1000*(s/55 + 1)*(s/1700 + 1)  );  % Get coefficients
den = sym2poly( s*(s/25000 + 1)               );  % ...
clear s;

% Construct controller TF
G = tf( num, den );

% Define a frequency array for loop shaping
wl = linspace( 1, 500000 );
L0 = P( 1, 1, nompt );
L0.ioDelay = 0; % no delay
lpshape( wl, ubdb, L0, G );

% [INFO] ...
fprintf( ACK );


%% Step 10: Synthesize Prefitler F(s)

% [INFO] ...
fprintf( 'Step 10:' );
fprintf( '\tSynthesize F(s)...' );

syms s;
num = 1;
den = sym2poly( s/550 + 1 );
clear s;

F = tf( num, den );

pfshape( 7, wl, del_6, P, [], G, [], F );

% [INFO] ...
fprintf( ACK );

%% Step 11-13: ANALYSIS

disp(' ')
disp('chksiso(1,wl,del_1,P,R,G); %margins spec')
chksiso( 1, wl, del_1, P, [], G );
ylim( [0 3.5] );

disp(' ')
disp('chksiso(2,wl,del_3,P,R,G); %Sensitivity reduction spec')
ind = find(wl <= 1000);
chksiso( 2, wl(ind), del_3, P, [], G );
ylim( [-90 10] );

disp(' ')
disp('chksiso(7,wl,W3,P,R,G); %input disturbance rejection spec')
drawnow
% chksiso(7,wl(ind),W7,P,[],G);
chksiso( 7, wl, del_6, P, [], G, [], F );
% ylim( [-0.1 1.3] );
