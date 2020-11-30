function [C, trialParamsHeader, exptParamsHeader, additionalRasterFields, badLogic, otherInstructions] = MemsaccCodes()

% This is a collection of constant definitions for use with analysis
% programs.  These codes are specific to the WhiteNoise paradigm

% The offset that is added to the value ecodes that follow the
% identifier ecodes (7000-9000 range).  Note, this value must be
% identical to the definition of VALODFFSET in rig1.h on the REX machine.
C.VALOFFSET = 4000;

% Ecodes that identify the header parameters
% (Parameters that do not change on a trial to trial basis
% and so are dropped once, at the begining of the run)
C.STIMSIZECD = 8000;
C.FPSIZECD = 8001;
C.FPRCD = 8002;
C.FPGCD = 8003;
C.FPBCD = 8004;
C.EYEWINXCD = 8005;
C.EYEWINYCD = 8006;
C.RFXCD = 8007;
C.RFYCD = 8008;
C.NTRLSPERCONDCD = 8009;
C.STIMRCD = 8010;
C.STIMGCD = 8011;
C.STIMBCD = 8012;

% Ecodes that identify the trial parameters
% (Parameters that change on a trial to trial basis
% and so are dropped at the end of each trial)
C.STIMXCD = 7000;
C.STIMYCD = 7001;

% Ecodes that are used to mark the times of various trial events
C.FPONCD = 1004;
C.FPOFFCD = 1006;
C.TARGONCD = 1008;
C.TARGOFFCD = 1009;
C.ABORTCD = 1013;
C.REWCD = 1015;
C.EOTCD = 1016;
C.SACMADCD = 1026;
C.STIMONCD = 1030;
C.STIMOFFCD = 1031;
C.FPACQCD = 1032;

%******************%
%   for nex2stro   %
%******************%

% includes the trial column identifiers, the ecode id's, the data type, and
% bookending info.
trialParamsHeader = {
    'fpon_t',     'time', 0, C.FPONCD
    'fpacq_t',    'time', 0, C.FPACQCD
    'fpoff_t',    'time', 0, C.FPOFFCD
    'stimon_t',   'time', 0, C.STIMONCD
    'stimoff_t',  'time', 0, C.STIMOFFCD
    'saccmade_t', 'time', 0, C.SACMADCD
    'stim_x',     'int',  0, C.STIMXCD
    'stim_y',     'int',  0, C.STIMYCD
    }';

%for the stro.sum.expParams array (includes things that get dropped only
%once):
exptParamsHeader = {
    'fp_size',     'int', 0, C.FPSIZECD
    'fp_r',        'int', 0, C.FPRCD
    'fp_g',        'int', 0, C.FPGCD
    'fp_b',        'int', 0, C.FPBCD
    'eyewin_x',    'int', 0, C.EYEWINXCD
    'eyewin_y',    'int', 0, C.EYEWINYCD
    'stim_size',   'int', 0, C.STIMSIZECD
    'stim_r',      'int', 0, C.STIMRCD
    'stim_g',      'int', 0, C.STIMGCD
    'stim_b',      'int', 0, C.STIMBCD
    'rf_x',      'int', 0, C.RFXCD
    'rf_y',      'int', 0, C.RFYCD
    'trlspercond', 'int', 0, C.NTRLSPERCONDCD
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
