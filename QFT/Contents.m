% QFT Feedback Control Design Toolbox.
% Version 2.6 (R14) 24-Jun-2004
%
% Specialized X-Y Graphs.
%   plotbnds  - Nichols plot of bounds.
%   plottmpl  - Nichols plot of templates.
%
% Bound Computation Managers.
%   sisobnds  - Single-Input/Single-Output bound computation.
%   genbnds   - General bound computation.
%
% Integrated Development Environments.
%   lpshape   - Controller design.
%   pfshape   - Pre-filter design.
%
% Analysis.
%   chksiso   - Compare design to SISO specifications.
%   chkgen    - Compare design to general specifications.
%
% Arithmetic.
%   addtmpl   - Addition of LTI sets.
%   cltmpl    - Closed-loop of LTI sets.
%   multmpl   - Multiplication of LTI sets.
%
% General Utility.
%   freqcp    - Continuous-time frequency response of num/den matrix.
%   dfreqcp   - Discrete-time frequency response of num/den matrix.
%   qftdefs   - User-defined defaults.
%   putqft    - Interactively set IDE MAT-file.
%   getqft    - Interactively retrieve IDE MAT-file.

% Copyright (c) 2003, Terasoft, Inc.

%All other non-user functions:
%
%MOUSE
%=====
%modisp   - motion function displaying current pointer location
%mogain   - motion function to edit the response gain
%mofirst  - motion function to add a real pole or zero term
%mosecond - motion function to add a complex pole or zero term
%moldlg   - motion function to add a lead/lag term
%montch   - motion function to add a notch term
%mo2ovr2  - motion function to add a second-order over second-order term
%qbtnkill - reset mouse functions while in Freq Weight IDE
%qbtnpres - mouse manager function while in Freq Weight IDE
%qbtnup   - button up manager in Freq Weight IDE
%qaddobj  - down function used for adding frequency weights in the FW IDE
%qmoveobj - motion function used for moving frequency weights in the FW IDE
%qzoomplt - zoom tool for PLOTBDS, PLOTTMP, CHEKGEN, and CHEKSISO
%qscrelmt - manager for mo* functions
%
%FREQUENCY WEIGHT BALANCING (Pepijn's)
%=====================================
%qde2sn    - transforms the (D,E) form to the (S,n) form
%qdf2sn    - transforms the complex (D,F) form to the complex (S,n) form
%qf2e      - transform from LBC to XBC form
%qlyaps    - solve the general form of the Lyapunov matrix equation
%qsimlbal  - minimal PQ-balanced realization of G
%qsnsys    - transforms a system to the (S,n) mu-format
%qrpk2gf   - SISO system to complex modal system realization
%qsym2def  - factorizes a symmetric matrix
%qunpack   - unpacks a system matrix in mu-format into equivalent forms
%xxxpep    - calculates Z = diax(diax(X)*diax(Y)) (diax does not exist
%            because it was merged into xxxpep)
%zpk2schr - zeros-poles-gain to real Schur system realization form
%
%INTERFACE SETUP
%===============
%qbtnsetup - interface for all IDE functions
%cntdisp   - interface for delete, edit, iterate, convert options
%            also, display present controller terms
%qmodlred  - interface for model order reduction
%qwatecad  - interface for frequency weight IDE
%qaxischng - interface for changing axis limits for IDE
%qfreqchng - interface for changing frequency vector for IDE
%envtowks  - interface for placing values into the workspace
%cntopen   - uigetfile used for retrieval of controller matrix
%cntsave   - uiputfile used to save a controller matrix
%
%SLIDER CALLBACKS
%================
%qitrelop  - iterate element values in frequency response
%
%BUTTON/MENU CALLBACKS
%=====================
%qdelelmt - delete elements from frequency response
%qedtelmt - edit elements in frequency response
%qitrelmt - iterate elements in frequency response
%cntcnvt  - convert pole/zero elements to lead/lag elements
%qundoit  - undo last frequency response operation
%qcadevuw - manager for QNICVUW and QMAGVUW
%qnicvuw  - axis limits manipulator for LPSHAPE
%qmagvuw  - axis limits manipulator for PFSHAPE
%bndonoff - toggle bound visibility property
%cntstor  - store present controller matrix in memory
%cntrecl  - recall stored controller matrix
%qelmts   - compute and store selected element from QADDELMT
%
%DEFAULT SETTING
%===============
%bndsdef   - defaults for SISOBNDS and GENBNDS
%chkdef    - defaults for CHKSISO and CHKGEN
%clcpdef   - defaults for CLCP
%clnddef   - defaults for CLND
%lpshpdef  - defaults for LPSHAPE
%pfshpdef  - defaults for PFSHAPE
%qplotdef  - defaults for PLOTBNDS
%
%FREQUENCY RESPONSE COMPUTATION
%==============================
%qcntbode  - response of controller matrix in complex
%rlroot    - first order pole or zero term in complex
%cproot    - second order pole or zero term in complex
%cintegtr  - continuous-time integrator term in complex
%dintegtr  - discrete-time integrator term in complex
%ldlgcplx  - lead/lag term in complex
%ntchcplx  - notch term in complex
%qcpqft    - response of num/den in complex
%
%FREQUENCY RESPONSE PLOTTING
%===========================
%qnicplt   - nichols response
%qmagplt   - magnitude response
%
%GENERAL PLOTTING
%================
%qplotbd   - plot bounds for PLOTBNDS, LPSHAPE, DLPSHAPE
%
%GENERAL COMPUTATION
%===================
%chkzp     - check for any unstable, repeated, or jw-axis poles or zeros
%chkstab   - stability of frequency response using state-space
%cnt2zpk   - controller matrix format to zero/pole/gain format
%cntcvrt   - convert notch and lead/lag terms to zero and pole terms
%cntdcgn   - D.C. gain of a controller matrix
%cntextr   - controller matrix format to num/den format
%cntpars   - num/den format to controller matrix format
%csecond   - continuous zeta and wn from phase/magnitude change
%dsecond   - discrete zeta and wn from phase/magnitude change
%qfindfrq  - frequency of response closest to mouse pointer
%qfindinf  - retrieve information stored at bottom of bound vector
%qfixfase  - unwrap frequency response depending upon current axis limits
%qfrqenh   - augment frequency vector depending second order terms
%qatan4    - 4-quadrant arctangent. returns value between -360 and 0
%qclassfy  - classify computed bounds to speed intersection
%qrobust   - robust stability bound for non-parametric uncertainty
%quadrtic  - specialized quadratic roots computation
%seconsec  - second-order over second-order from phase/magnitude change
%zp2cnt    - zero/pole/gain format to controller matrix format
%zp2ldlg   - zero and pole to lead/lag term (phase and frequency)
%qsubset   - determine intersection of two vectors
%
%HELP FILES
%==========
%qfile.hlp   - options under QFile menu
%qtools.hlp  - options under QTools menu
%qview.hlp   - options under QView menu
%bndshelp.txt - definition of a bound (used in demonstration facility)
%tmplhelp.txt - definition of a template (used in demonstration facility)
%
%OTHER
%=====
%elmtbutn - various element button choices in QADDELMT depending upon IDE
%qputobj  - uiputfile/uigetfile for frequency weights
%qdelobj  - function used for deleting frequency weights in the FW IDE
%qedtobj  - function used for editing frequency weights in the FW IDE
%repltest - test for what direction MATLAB replicates using ':'
%copybnds - copy bounds to other phase if axis limits changed
%mesgbnds - different messages displayed concerning status of bound
%           computation
%qclswin  - reset interface within LPSHAPE/PFSHAPE
%xitcade  - close all windows related to current IDE
%qpause   - special pause that allows for use of IDE within example files
%nxtstage - script file that sets up call to next stage of demonstration
%presexit - setup for exiting from IDE within demonstration
%qfterror - executed if an error occurs in a bound computation function
