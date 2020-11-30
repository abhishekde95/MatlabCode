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
C.MONDISTCD = 8021;
C.FPMOVECD = 8023;
C.STIMXNEARCD = 8024;
C.STIMYNEARCD = 8025;
C.STIMXFARCD = 8026;
C.STIMYFARCD = 8027;
C.ALLOWEQRADIICD = 8028;
C.TARGWINXCD = 8029;
C.TARGWINYCD = 8030;
C.NDIAMSCD = 8031;
C.WRONGTARGCONTCD = 8032;
C.STIMEQUIDISTCD = 8033;
C.STIMMOVECD = 8034;
C.EYEPOSCONTCD = 8035;

% Ecodes that identify the trial parameters
% (Parameters that change on a trial to trial basis
% and so are dropped at the end of each trial)
C.FPXCD = 7000;
C.FPYCD = 7001;
C.STIMIRADCD = 7002;
C.STIMORADCD = 7003;
C.STIMCONFIGCD = 7004;
C.STIMTYPECD = 7005;
C.STIMIDNEARCD = 7006;
C.STIMODNEARCD = 7007;
C.STIMIDFARCD = 7008;
C.STIMODFARCD = 7009;
C.STIMXACT_NEARCD = 7010;
C.STIMYACT_NEARCD = 7011; 
C.STIMXACT_FARCD = 7012;
C.STIMYACT_FARCD = 7013;

% Ecodes that are used to mark the times of various trial events
C.FPONCD = 1004;
C.TARGONCD = 1008;
C.TARGOFFCD = 1009;
C.FPACQCD = 1032;
C.STIMONCD = 1030;
C.STIMOFFCD = 1031;
C.FPOFFCD = 1006;
C.REWCD = 1015;
C.EOTCD = 1016;
C.CORRECTCD = 1018;
C.ABORTCD = 1013;
C.SACMADCD = 1026;
C.SACMADENEARCD = 1050;
C.SACMADEFARCD = 1051;

%******************%
%   for nex2stro   %
%******************%

% includes the trial column identifiers, the ecode id's, the data type, and
% bookending info.
trialParamsHeader = {
    'fp_on',       'time', 0, C.FPONCD
    'fp_acq',      'time', 0, C.FPACQCD
    'stim_on',      'time', 0, C.STIMONCD
    'stim_off',     'time', 0, C.STIMOFFCD
    'targ_on', 'time', 0, C.TARGONCD
    'fp_off',    'time', 0, C.FPOFFCD
    'sacc_time', 'time', 0, C.SACMADCD
    'sacc_time_near', 'time', 0, C.SACMADENEARCD
    'sacc_time_far', 'time', 0, C.SACMADEFARCD
    'fp_x',     'int', 0, C.FPXCD
    'fp_y',    'int',	0, C.FPYCD
    'stim_ir',    'float',	0, C.STIMIRADCD
    'stim_or',    'float',	0, C.STIMORADCD
    'stimx_near', 'int', 0, C.STIMXACT_NEARCD
    'stimy_near', 'int', 0, C.STIMYACT_NEARCD
    'stimx_far', 'int', 0, C.STIMXACT_FARCD
    'stimy_far', 'int', 0, C.STIMYACT_FARCD
    'bkgnd', 'int',	0, C.STIMCONFIGCD
    'stim_ir_near', 'float', 0, C.STIMIDNEARCD
    'stim_or_near', 'float', 0, C.STIMODNEARCD
    'stim_ir_far', 'float', 0, C.STIMIDFARCD
    'stim_or_far', 'float', 0, C.STIMODFARCD
    'correct_resp', 'time', 0, C.CORRECTCD
    }';

%for the stro.sum.expParams array (includes things that get dropped only
%once):
exptParamsHeader = {
    'fp_size',        'int', 0, C.FPSIZECD
    'fp_r',        'int', 0, C.FPRCD
    'fp_g',     'int', 0, C.FPGCD
    'fp_b',        'int', 0, C.FPBCD
    'eyewin_x',        'int', 0, C.EYEWINXCD
    'eyewin_y',        'int', 0, C.EYEWINYCD
    'targwin_x', 'int', 0, C.TARGWINXCD
    'targwin_y', 'int', 0, C.TARGWINYCD
    'rf_x',    'int', 0, C.RFXCD
    'rf_y',    'int', 0, C.RFYCD
    'bkgndrgb',   'double', 1, C.BKGNDRGBCD
    'framerate',      'double', 0, C.FRAMERATECD
    'trialspercond',      'int', 0, C.TRLSPERCONDCD
    'nframes',      'int', 0, C.NFRAMESPERDISPLAYCD
    'pixperdeg',   'double', 0, C.PIXPERDEGCD
    'ndiams', 'int', 0, C.NDIAMSCD
    'stimr',   'int', 0, C.STIMRCD
    'stimg',        'int', 0, C.STIMGCD
    'stimb',        'int', 0, C.STIMBCD
    'mondist', 'int', 0, C.MONDISTCD
    'fpmove', 'int', 0, C.FPMOVECD
    'stimx_near', 'int', 0, C.STIMXNEARCD
    'stimy_near', 'int', 0, C.STIMYNEARCD
    'stimx_far', 'int', 0, C.STIMXFARCD
    'stimy_far', 'int', 0, C.STIMYFARCD
    'allow_eq_radii', 'int', 0, C.ALLOWEQRADIICD
    'wrong_targ_cont', 'int', 0, C.WRONGTARGCONTCD
    'stim_equidist', 'int', 0, C.STIMEQUIDISTCD
    'stim_move', 'int', 0, C.STIMMOVECD
    'eye_pos_control', 'int', 0, C.EYEPOSCONTCD
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