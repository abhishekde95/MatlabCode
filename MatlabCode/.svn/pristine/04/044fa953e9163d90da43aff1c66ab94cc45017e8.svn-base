function [C, trialParamsHeader, exptParamsHeader, additionalRasterFields, badLogic, otherInstructions] = DTonelocCodes()

% This is a collection of constant definitions for use with analysis
% programs.  These codes are specific to the WhiteNoise paradigm

% The offset that is added to the value ecodes that follow the
% identifier ecodes (7000-9000 range).  Note, this value must be
% identical to the definition of VALODFFSET in rig1.h on the REX machine.
C.VALOFFSET = 4000;

% Ecodes that identify the header parameters
% (Parameters that do not change on a trial to trial basis
% and so are dropped once, at the begining of the run)
C.FIXXCD =           8000;
C.FIXYCD =           8001;
C.FPSIZECD =         8002;
C.EYEWINWCD =        8003;
C.EYEWINHCD =        8004;
C.TARGSIZECD =       8005;
C.TARGWINWCD =       8006;
C.TARGWINHCD =       8007;
C.RFXCD =            8008;
C.RFYCD =            8009;
C.NSTIMCD =          8010;
C.TRIALSPERBLOCKCD = 8011;
C.NBLOCKSCD =        8012;
C.FRAMERATECD =      8013; % double, 0
C.PIXPERDEGCD =      8014; % double, 0
C.BKGNDRGBCD =       8015; % double, 1
C.GAMMATABLECD =     8016; % double, 1
C.FUNDSCD =          8017; % double, 1
C.MONSPDCD =         8018; % double, 1
C.LMSSTAIRCD =       8019; % float, 1
C.WRONGTARGCONTCD =  8020; % float, 0
C.THETACD =          8023;
C.SFCD =             8024;
C.PHICD =            8025;
C.GAMMACD =          8026;
C.SIGMACD =          8027;
C.NSIGMASCD =        8028;
C.MODULATORCD =      8029;
C.RAMPTIMECD =       8030;
C.CRREWMULTCD =      8031;
C.HEADERREADYCD =    8998;

% Ecodes that identify the trial parameters
% (Parameters that change on a trial to trial basis
% and so are dropped at the end of each trial)
C.STIMPRESENTCD =  7000;
C.LCCCD =          7001;
C.MCCCD =          7002;
C.SCCCD =          7003;
C.TFCD =           7004;
C.OOGCD =          7005;
C.STIMIDXCD =      7006;
C.LASERCD =        7007;
C.CORRECTTRIALCD = 7008;
C.STIMXCD =        7009;
C.STIMYCD =        7010;

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
C.EYINWDCD =  1027;
C.STIMONCD =  1030;
C.STIMOFFCD = 1031;
C.FPACQCD =   1032;
C.MSTIMONCD = 1033;
C.MSTIMOFFCD= 1034;
C.EXTRAREWCD= 1035;

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
    'saccstart_t',  'time', 0, C.SACMADCD
    'saccend_t',    'time', 0, C.EYINWDCD
    'rew_t',        'time', 0, C.REWCD
    'laseron_t',     'time', 0, C.MSTIMONCD
    'laseroff_t',    'time', 0, C.MSTIMOFFCD
    'stimpresent',   'int', 0, C.STIMPRESENTCD
    'lcc',        'double', 0, C.LCCCD
    'mcc',        'double', 0, C.MCCCD
    'scc',        'double', 0, C.SCCCD
    'tf',         'double', 0, C.TFCD
    'oog',           'int', 0, C.OOGCD
    'stim_idx',      'int', 0, C.STIMIDXCD
    'optstim',       'int', 0, C.LASERCD
    'correct',       'int', 0, C.CORRECTTRIALCD
    'stim_x',        'int', 0, C.STIMXCD
    'stim_y',        'int', 0, C.STIMYCD
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
    'rf_x',            'int', 0, C.RFXCD
    'rf_y',            'int', 0, C.RFYCD
    'n_stim',          'int', 0, C.NSTIMCD
    'trials_block',    'int', 0, C.TRIALSPERBLOCKCD
    'nblocks',         'int', 0, C.NBLOCKSCD
    'framerate',    'double', 0, C.FRAMERATECD
    'pixperdeg',    'double', 0, C.PIXPERDEGCD
    'bkgndrgb',     'double', 1, C.BKGNDRGBCD
    'gamma_table',  'double', 1, C.GAMMATABLECD
    'fundamentals', 'double', 1, C.FUNDSCD
    'mon_spd',      'double', 1, C.MONSPDCD
    'theta',        'float',  0, C.THETACD
    'sf',           'float',  0, C.SFCD
    'phi',          'float',  0, C.PHICD
    'gamma',        'float',  0, C.GAMMACD
    'sigma',        'float',  0, C.SIGMACD
    'nsigmas',      'float',  0, C.NSIGMASCD
    'wrongtargscale','float', 0, C.WRONGTARGCONTCD
    'use_modulator',   'int', 0, C.MODULATORCD
    'cr_rew_mult',   'float', 0, C.CRREWMULTCD
    }';

additionalRasterFields = {};

%define how good and bad trials are to be distinguished
badLogic = 'find(trialCodes(:,2) == C.ABORTCD, 1)';

%define other routines as an array of function handles. the first row of
%the array should be a name of the handle, and the second row should be the
%actual function handle. The functions should be included as non-nested
%subfunctions in this file. DTStimCodes.m serves as a template
otherInstructions = {}; %none for now
