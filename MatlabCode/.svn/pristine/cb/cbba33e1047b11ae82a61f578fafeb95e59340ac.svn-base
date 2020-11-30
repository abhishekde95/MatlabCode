function [C, trialParamsHeader, exptParamsHeader, additionalRasterFields, badLogic, otherInstructions] = SMurrayCodes()

% This is a collection of constant definitions for use with analysis
% programs.  These codes are specific to the WhiteNoise paradigm

% The offset that is added to the value ecodes that follow the
% identifier ecodes (7000-9000 range).  Note, this value must be
% identical to the definition of VALODFFSET in rig1.h on the REX machine.
C.VALOFFSET = 4000;

% Ecodes that identify the header parameters
% (Parameters that do not change on a trial to trial basis
% and so are dropped once, at the begining of the run)
C.FIXXCD = 8000;
C.FIXYCD = 8001;
C.FPSIZECD = 8002;
C.FPRCD = 8003;
C.FPGCD = 8004;
C.FPBCD = 8005;
C.EYEWINXCD = 8006;
C.EYEWINYCD = 8007;
C.RFXCD = 8008;
C.RFYCD = 8009;
C.BKGNDRGBCD = 8010; %double,1
C.FRAMERATECD = 8011; %double, 0
C.TRLSPERCONDCD = 8012;
C.NFRAMESPERDISPLAYCD = 8013;
C.PIXPERDEGCD = 8014; % double, 1
C.STIMRCD = 8015; 
C.STIMGCD = 8016; 
C.STIMBCD = 8017; 
C.PREFORIENTCD = 8018;
C.PREFSFCD = 8019;
C.PREFDIAMCD = 8020;
C.NSIGMASPERDIAMCD = 8021;
C.PHASECD = 8022;


% Ecodes that identify the trial parameters
% (Parameters that change on a trial to trial basis
% and so are dropped at the end of each trial)
C.TARGSHOWNCD = 7000;
C.FLANKSHOWNCD = 7001;
C.FLANKINLINECD = 7002;
C.CONFIGCD = 7003;
C.CONTRASTCD = 7004;
C.FLANKERDISTCD = 7005;
C.FLANKFLANKSHOWNCD = 7006;
C.FLANKFLANKINLINECD = 7007;

% Ecodes that are used to mark the times of various trial events
C.FPONCD = 1004;
C.FPACQCD = 1032;
C.MASKONCD = 1010;
C.STIMONCD = 1030;
C.ALLOFFCD = 1012;
C.FPOFFCD = 1006;
C.REWCD = 1015;
C.EOTCD = 1016;
C.ABORTCD = 1013;

%******************%
%   for nex2stro   %
%******************%

% includes the trial column identifiers, the ecode id's, the data type, and
% bookending info.
trialParamsHeader = {'fp_on', 'fp_acq', 'flankers_on', 'stim_on', 'stim_off', 'fp_off', 'targ', 'flank', 'flankinline', 'flankflank', 'flankflankinline', 'contrast', 'config', 'flankerdist';...
                     'time',  'time',   'time',        'time',    'time',     'time',   'int',  'int',   'int',         'int',        'int',              'float',    'int',    'int';...
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;  ...
                       C.FPONCD, C.FPACQCD, C.MASKONCD, C.STIMONCD, C.ALLOFFCD, C.FPOFFCD, C.TARGSHOWNCD, C.FLANKSHOWNCD, C.FLANKINLINECD, C.FLANKFLANKSHOWNCD, C.FLANKFLANKINLINECD, C.CONTRASTCD, C.CONFIGCD, C.FLANKERDISTCD};

%for the stro.sum.expParams array (includes things that get dropped only
%once):
exptParamsHeader = {'fp_x', 'fp_y', 'fp_size', 'fp_r', 'fp_g', 'fp_b', 'eyewin_x', 'eyewin_y', 'rf_x', 'rf_y', 'bkgndrgb', 'framerate', 'trialspercond','nframes', 'pixperdeg', 'orient', 'sf',   'diam', 'nsigmasperdiam', 'phase';...
                    'int',  'int',  'int',     'int',  'int',  'int',  'int',      'int',      'int',  'int',  'double',   'double',    'int',          'int',     'double',    'float',  'float','float','float',          'float';...
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0 ,0 0;...
                    C.FIXXCD, C.FIXYCD, C.FPSIZECD, C.FPRCD, C.FPGCD, C.FPBCD, C.EYEWINXCD, C.EYEWINYCD, C.RFXCD, C.RFYCD, C.BKGNDRGBCD, C.FRAMERATECD, C.TRLSPERCONDCD, C.NFRAMESPERDISPLAYCD, C.PIXPERDEGCD, C.PREFORIENTCD, C.PREFSFCD, C.PREFDIAMCD, C.NSIGMASPERDIAMCD, C.PHASECD};

additionalRasterFields = {};


%define how good and bad trials are to be distinguished
badLogic = 'find(trialCodes(:,2) == C.ABORTCD, 1)';

%define other routines as an array of function handles. the first row of
%the array should be a name of the handle, and the second row should be the
%actual function handle. The functions should be included as non-nested
%subfunctions in this file. DTStimCodes.m serves as a template
otherInstructions = {}; %none for now

end

