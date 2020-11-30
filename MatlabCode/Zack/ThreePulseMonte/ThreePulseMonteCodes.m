function [C, trialParamsHeader, exptParamsHeader, additionalRasterFields, badLogic, otherInstructions] = ThreePulseMonteCodes()

% The offset that is added to the value ecodes that follow the
% identifier ecodes (7000-9000 range).  Note, this value must be
% identical to the definition of VALOFFSET in rigconsts.h on the REX machine.
C.VALOFFSET = 4000;

% Ecodes that identify the header parameters
% (Parameters that do not change on a trial to trial basis
% and so are dropped once, at the begining of the run)
C.FIXXCD =     8000;
C.FIXYCD =     8001;
C.FPSIZECD =   8002;
C.EYEWINWCD =  8003;
C.EYEWINHCD =  8004;
C.TARGXCD =    8005;
C.TARGYCD =    8006;
C.TARGSIZECD = 8007;
C.TARGWINWCD = 8008;
C.TARGWINHCD = 8009;
C.STIMDURCD =  8010;
C.NPULSESCD =  8011;
C.OPTPULSECD = 8012;
C.NBLOCKSCD =  8013;
C.MMATRIXCD =  8014;
C.BKGNDRGBCD = 8015;

% Ecodes that identify the trial parameters
% (Parameters that change on a trial to trial basis
% and so are dropped at the end of each trial)
C.PULSE2OFFSETCD =  7000;
C.CORRECTTRIALCD =  7001;
C.THRESHMULTCD   =  7002;
C.LCCCD =           7003;
C.MCCCD =           7004;
C.CORRTARGXCD =     7005;
C.CORRTARGYCD =     7006;
C.RIGHTWARDSACCCD = 7007;
C.TFCD =            7008;

% Ecodes that are used to mark the times of various trial events
C.FPONCD =    1004;
C.FPOFFCD =   1006;
C.TARGONCD =  1008;
C.ABORTCD =   1013;
C.ERRCD =     1014;
C.REWCD =     1015;
C.EOTCD =     1016;
C.CORRECTCD = 1018;
C.SACMADCD =  1026;
C.FPACQCD =   1032;
C.TPULSE1CD = 1040;
C.TPULSE2CD = 1041;
C.TPULSE3CD = 1042;

%******************%
%   for nex2stro   %
%******************%

% includes the trial column identifiers, the ecode id's, the data type, and
% bookending info.
trialParamsHeader = {
    'fpon_t',        'time', 0, C.FPONCD
    'fpacq_t',       'time', 0, C.FPACQCD
    'fpoff_t',       'time', 0, C.FPOFFCD
    'pulse1_t',      'time', 0, C.TPULSE1CD
    'pulse2_t',      'time', 0, C.TPULSE2CD
    'pulse3_t',      'time', 0, C.TPULSE3CD
    'targon_t',      'time', 0, C.TARGONCD
    'saccmade_t',    'time', 0, C.SACMADCD
    'rew_t',         'time', 0, C.REWCD
    'pulse_offset',   'int', 0, C.PULSE2OFFSETCD
    'correct',        'int', 0, C.CORRECTTRIALCD
    'thresh_mult',  'float', 0, C.THRESHMULTCD
    'lcc',          'float', 0, C.LCCCD
    'mcc',          'float', 0, C.MCCCD
    'tf',           'float', 0, C.TFCD
    'corrtarg_x',     'int', 0, C.CORRTARGXCD
    'corrtarg_y',     'int', 0, C.CORRTARGYCD
    'rightward_sacc', 'int', 0, C.RIGHTWARDSACCCD
    }';

%for the stro.sum.exptParams array (includes things that get dropped only
%once):
exptParamsHeader = {
    'fp_x',            'int', 0, C.FIXXCD
    'fp_y',            'int', 0, C.FIXYCD
    'fp_size',         'int', 0, C.FPSIZECD
    'eyewin_w',        'int', 0, C.EYEWINWCD
    'eyewin_h',        'int', 0, C.EYEWINHCD
    'targ_size',       'int', 0, C.TARGSIZECD
    'targwin_w',       'int', 0, C.TARGWINWCD
    'targwin_h',       'int', 0, C.TARGWINHCD
    'targ_x',          'int', 0, C.TARGXCD
    'targ_y',          'int', 0, C.TARGYCD
    'nblocks',         'int', 0, C.NBLOCKSCD
    'npulses',         'int', 0, C.NPULSESCD
    'optional_pulse',  'int', 0, C.OPTPULSECD
    'stim_dur',        'int', 0, C.STIMDURCD
    'M',            'double', 1, C.MMATRIXCD
    'bkgndrgb',     'double', 1, C.BKGNDRGBCD
    }';

additionalRasterFields = {};

%define how good and bad trials are to be distinguished
badLogic = 'find(trialCodes(:,2) == C.ABORTCD, 1)';

%define other routines as an array of function handles. the first row of
%the array should be a name of the handle, and the second row should be the
%actual function handle. The functions should be included as non-nested
%subfunctions in this file. DTStimCodes.m serves as a template
otherInstructions = {}; %none for now
