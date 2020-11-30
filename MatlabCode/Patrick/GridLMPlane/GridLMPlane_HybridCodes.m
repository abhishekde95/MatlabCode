function [C, trialParamsHeader, exptParamsHeader, additionalRasterFields, badLogic, otherInstructions] = GridLMPlane_HybridCodes()

% This is a collection of constant definitions for use with analysis
% programs.  These codes are specific to the GridLMPlane paradigm

% 10/11     Created             JPW
% 2/9/12    Updated (added DM)  JPW
% 10/12/12  Updated (added STMPARTPCD) JPW
% 10/22/12  Updated (Correcting L/M switch that occurred in datasets
%       prior to 10/12/12)      JPW
% 10/23/12  Updated (Trying a simpler fix to previous problem. L and M
%       codes now switched (permanently)    JPW
% 10/23/12  New File Created for GridLMPlane_Hybrid    JPW


% The offset that is added to the value ecodes that follow the
% identifier ecodes (7000-9000 range).  Note, this value must be
% identical to the definition of VALODFFSET in rig1.h on the REX machine.
C.VALOFFSET = 4000;

% Ecodes that identify the header parameters
% (Parameters that do not change on a trial to trial basis
% and so are dropped once, at the begining of the run)
C.FIXXCD =          8000;
C.FIXYCD =          8001;
C.FPSIZECD =        8002;
C.FPRCD =           8003;
C.FPGCD =           8004;
C.FPBCD =           8005;
C.EYEWINXCD =       8006;
C.EYEWINYCD =       8007;
C.FRAMERATECD =     8008;
C.BKGNDRGBCD =      8009;
C.PIXPERDEGCD =     8010;
C.GAMMATABLECD =    8011;
C.FUNDAMENTALSCDC = 8012;
C.MONSPDCD =        8013;
C.RFXCD =           8014;
C.RFYCD =           8015;
C.THETACD =         8016;
C.SIGMAXCD =        8017;
C.SIGMAYCD =        8018;
C.DRIFTRATECD =     8019;
C.NSTDCD =          8020;
C.NEXPANSECD =      8021;

C.HDRCOMPLETECD =   8998;

% Ecodes that identify the trial parameters
% (Parameters that change on a trial to trial basis
% and are therefore dropped at the end of each trial)

C.LCCCD =       7000;
C.MCCCD =       7001;
C.SCCCD =       7002;
C.THETATPCD =   7003;
C.SIGMAXTPCD =  7004;
C.SIGMAYTPCD =  7005;
C.DRTPCD =      7006;
C.NSTDTPCD =    7007;
C.NEXPTPCD =    7008;
C.ADAPTIVECD =  7009;

% Ecodes that are used to mark the times of various trial events
C.FPONCD =      1004;
C.FPOFFCD =     1006;
C.TARGONCD =	1008;
C.TARGOFFCD =	1009;
C.ABORTCD =     1013;
C.REWCD =       1015;
C.EOTCD =       1016;
C.STIMONCD =    1030;
C.STIMOFFCD =	1031;
C.FPACQCD =     1032;

%******************%
%   for nex2stro   %
%******************%

% includes the trial column identifiers, the ecode id's, the data type, and
% bookending info.
trialParamsHeader = {
    'fpon_t',           'time', 0, C.FPONCD
    'fpacq_t',          'time', 0, C.FPACQCD
    'fpoff_t',          'time', 0, C.FPOFFCD
    'stimon_t',         'time', 0, C.STIMONCD
	'Lcc',      		'float', 0, C.LCCCD
    'Mcc',              'float', 0, C.MCCCD
    'Scc',              'float', 0, C.SCCCD
    'Theta',            'float', 0, C.THETATPCD
    'Sigma X',          'float', 0, C.SIGMAXTPCD
    'Sigma Y',          'float', 0, C.SIGMAYTPCD
    'Drift Rate',       'float', 0, C.DRTPCD
    'nstd',             'float', 0, C.NSTDTPCD
    'nexp',             'float', 0, C.NEXPTPCD
    'Adaptive Paradigm','float', 0, C.ADAPTIVECD
    }';

%for the stro.sum.expParams array (includes things that get dropped only
%once):
exptParamsHeader = {
    'fp_x',         'int', 0, C.FIXXCD
    'fp_y',         'int', 0, C.FIXYCD
    'fp_size',      'int', 0, C.FPSIZECD
    'fp_r',         'int', 0, C.FPRCD
    'fp_g',         'int', 0, C.FPGCD
    'fp_b',         'int', 0, C.FPBCD
    'eyewin_x',     'int', 0, C.EYEWINXCD
    'eyewin_y',     'int', 0, C.EYEWINYCD
    'rf_x',         'float', 0, C.RFXCD
    'rf_y',         'float', 0, C.RFYCD
    'theta',        'float', 0, C.THETACD
    'sigma_X',      'float', 0, C.SIGMAXCD
    'sigma_Y',      'float', 0, C.SIGMAYCD
    'driftrate',    'float', 0, C.DRIFTRATECD
    'nstd',         'float', 0, C.NSTDCD
    'nexpanse',     'float', 0, C.NEXPANSECD
    }';

additionalRasterFields = {};

%define how good and bad trials are to be distinguished
badLogic = 'find(trialCodes(:,2) == C.ABORTCD, 1)';

%define other routines as an array of function handles. the first row of
%the array should be a name of the handle, and the second row should be the
%actual function handle. The functions should be included as non-nested
%subfunctions in this file. DTStimCodes.m serves as a template
otherInstructions = {}; %none for now

end
