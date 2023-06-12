function defs = qftdefs
% QFTDEFS  User-defined QFT Toolbox default values.
%
%          The purpose of this file is to allow the user to change defaults
%          related to bound computations and the shaping environment.
%
%          Use your text editor to change the values within this file for
%          default values which you wish to customize.

% Author: Craig Borghesani
% 10/1/94
% Copyright (c) 2003, Terasoft, Inc.


% Phase vector (must be between -360 and 0 degrees)
% default = [0 : -5 : -360]

min_phase = -360;  % minimum phase value
max_phase =    0;  % maximum phase value
phase_res =   -5;  % phase resolution

% Nominal plant and controller indices
% default = [1, 1]

plant_index      = 1;
controller_index = 1;

% Controller type (specify which controller is UNKNOWN in the loop)
% 1 = G
% 2 = H
% default = 1

controller_type = 1;

% Frequency vector (used in shaping environments)
% values are in powers of 10
% default = logspace(-2,3,100)

min_w =   -2;
max_w =    3;
w_len =  100;

%
% END OF DEFAULT SECTION
%

defs = [max_phase, phase_res, min_phase;
        plant_index, controller_index, 0;
        controller_type, 0, 0;
        min_w, max_w, w_len];
