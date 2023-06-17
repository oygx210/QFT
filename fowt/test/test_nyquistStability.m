% Testing Nyquist stability code
%
%   AUTHOR  : Mohammad Odeh
%   DATE    : Jun. 16th, 2023
%
% CHANGELOG :
%   Jun. 16th, 2023
%       - Initial script
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
PRNT = true;                                %#ok<NASGU>
% PRNT = false;                               % COMMENT OUT TO PRINT FIGURES

% Get current path to working directory and split
pathParts = strsplit( pwd, filesep );
% Go up one level and generate new path
src = fullfile( pathParts{1:end-1} );

% If on a UNIX machine (i.e. macOS, Ubuntu, etc...), fix path since
% strsplit() removes the leading '/' from the path.
src = [ filesep src ];

% Add QFT2 to path
addpath( genpath(src) );

%% Nyquist plot testing of cases

% Testing L1
syms s;
num = 140.* sym2poly( (-0.5*s+1)*(-0.5714*s+1)*(-0.6*s+1)*(-3.5*s+1) );
den = sym2poly( s^5*(-5*s+1) );
clear s;

L1 = tf( num, den );
% figure(); nyquist( L1 );
% figure(); nichols( L1 );

[zc, N, num_p_RHP, Na, Nb, Nc, Nd, zpCancel, k, sigm, alpha, gamma] = nyquistStability(L1)

% Testing L8
syms s;
num = (0.0071) .* sym2poly( (-5*s+1) );
den = sym2poly( s^3*(-0.5*s+1)*(-0.5714*s+1)*(-0.6*s+1)*(-3.5*s+1) );
clear s;

L8 = tf( num, den );
% figure(); nyquist( L8 );
% figure(); nichols( L8 );

[zc, N, num_p_RHP, Na, Nb, Nc, Nd, zpCancel, k, sigm, alpha, gamma] = nyquistStability(L8)
