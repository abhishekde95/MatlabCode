function [C, trialParamsHeader, exptParamsHeader, additionalRasterFields, badLogic, otherInstructions] = SCStimcueCodes()

% This is a collection of constant definitions for use with analysis
% programs.  These codes are specific to the WhiteNoise paradigm

% The offset that is added to the value ecodes that follow the
% identifier ecodes (7000-9000 range).  Note, this value must be
% identical to the definition of VALODFFSET in rig1.h on the REX machine.
C.VALOFFSET = 4000;

% Ecodes that identify the header parameters
% (Parameters that do not change on a trial to trial basis
% and so are dropped once, at the begining of the run)
C.FIXXCD         = 8000;
C.FIXYCD         = 8001;
C.FPSIZECD       = 8002;
C.FPRCD          = 8003;
C.FPGCD          = 8004;
C.FPBCD          = 8005;
C.EYEWINXCD      = 8006;
C.EYEWINYCD      = 8007;
C.TARGSIZECD     = 8008;
C.TARGRCD        = 8009;
C.TARGGCD        = 8010;
C.TARGBCD        = 8011;
C.TARGWINXCD     = 8012;
C.TARGWINYCD     = 8013;
C.RFXCD          = 8014;
C.RFYCD          = 8015;
C.NTRLSPERCONDCD = 8016;
C.MAXCTACD       = 8017;
C.STIMTYPECD     = 8018;
C.DOTSPEEDCD     = 8019;
C.DOTSIZECD      = 8020;
C.NDOTSCD        = 8021;
C.APSIZECD       = 8022;

% Ecodes that identify the trial parameters
% (Parameters that change on a trial to trial basis
% and so are dropped at the end of each trial)
C.OPTSTIMFREQCD = 7000;
C.TARGXCD       = 7001;
C.TARGYCD       = 7002;
C.CTACD         = 7003;
C.OPTSTIMCD     = 7004;
C.VISSTIMCD     = 7005;
C.DOTDIRCD      = 7006;
C.CORTRIALCD    = 7007;

% Ecodes that are used to mark the times of various trial events
C.FPONCD    = 1004;
C.FPOFFCD   = 1006;
C.TARGONCD  = 1008;
C.TARGOFFCD = 1009;
C.ABORTCD   = 1013;
C.REWCD     = 1015;
C.EOTCD     = 1016;
C.SACMADCD  = 1026;
C.BROKEFIX  = 1029;
C.STIMONCD  = 1030;
C.STIMOFFCD = 1031;
C.FPACQCD   = 1032;

%******************%
%   for nex2stro   %
%******************%

% includes the trial column identifiers, the ecode id's, the data type, and
% bookending info.
trialParamsHeader = {
    'fpon_t',       'time', 0, C.FPONCD
    'fpacq_t',      'time', 0, C.FPACQCD
    'fpoff_t',      'time', 0, C.FPOFFCD
    'stimon_t',     'time', 0, C.STIMONCD
    'stimoff_t',    'time', 0, C.STIMOFFCD
    'targon_t',     'time', 0, C.TARGONCD
    'targoff_t',    'time',	0, C.TARGOFFCD
    'saccinit_t',   'time', 0, C.BROKEFIX
    'saccmade_t',   'time', 0, C.SACMADCD
    'optstim',       'int', 0, C.OPTSTIMCD
    'optstim_freq',  'int',	0, C.OPTSTIMFREQCD
    'targ_x',        'int',	0, C.TARGXCD
    'targ_y',        'int',	0, C.TARGYCD
    'cta',           'int', 0, C.CTACD
    'visstim',       'int', 0, C.VISSTIMCD
    'direction',     'int', 0, C.DOTDIRCD
    'correct',       'int', 0, C.CORTRIALCD
    }';


%for the stro.sum.expParams array (includes things that get dropped only
%once):
exptParamsHeader = {
    'fp_x',        'int', 0, C.FIXXCD
    'fp_y',        'int', 0, C.FIXYCD
    'fp_size',     'int', 0, C.FPSIZECD
    'fp_r',        'int', 0, C.FPRCD
    'fp_g',        'int', 0, C.FPGCD
    'fp_b',        'int', 0, C.FPBCD
    'eyewin_x',    'int', 0, C.EYEWINXCD
    'eyewin_y',    'int', 0, C.EYEWINYCD
    'targ_size',   'int', 0, C.TARGSIZECD
    'targ_r',      'int', 0, C.TARGRCD
    'targ_g',      'int', 0, C.TARGGCD
    'targ_b',      'int', 0, C.TARGBCD
    'targwin_x',   'int', 0, C.TARGWINXCD
    'targwin_y',   'int', 0, C.TARGWINYCD
    'rf_x',        'int', 0, C.RFXCD
    'rf_y',        'int', 0, C.RFYCD
    'trlspercond', 'int', 0, C.NTRLSPERCONDCD
    'maxcta',      'int', 0, C.MAXCTACD
    'stimtype',    'int', 0, C.STIMTYPECD
    'dotspeed',    'int', 0, C.DOTSPEEDCD
    'dotsize',     'int', 0, C.DOTSIZECD
    'ndots',       'int', 0, C.DOTSIZECD
    'ap_size',     'int', 0, C.APSIZECD  
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
