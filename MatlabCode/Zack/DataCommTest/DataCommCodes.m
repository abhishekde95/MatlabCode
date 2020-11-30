function [C, trialParamsHeader, exptParamsHeader, additionalRasterFields, badLogic, otherInstructions] = DataCommCodes()
% This is a collection of constant definitions for use with analysis
% programs.  These codes are specific to the WhiteNoise paradigm

% The offset that is added to the value ecodes that follow the
% identifier ecodes (7000-9000 range).  Note, this value must be
% identical to the definition of VALODFFSET in rig1.h on the REX machine.
C.VALOFFSET = 4000;

% Ecodes that identify the header parameters
% (Parameters that do not change on a trial to trial basis
% and so are dropped once, at the begining of the run)

C.HD.ICD =    8000;
C.HD.DCD =    8001;
C.HD.FCD =    8002;
C.HD.LCD =    8003;
C.HD.CCD =    8004;
C.HD.SCD =    8005;
C.HD.USCD =   8006;
C.HD.ULCD =   8007;
C.HD.UICD =   8008;
C.HD.LLCD =   8009;
C.HD.ULLCD =  8010;
C.HD.AICD =   8011;
C.HD.ADCD =   8012;
C.HD.AFCD =   8013;
C.HD.ALCD =   8014;
C.HD.ACCD =   8015;
C.HD.ASCD =   8016;
C.HD.AUSCD =  8017;
C.HD.AULCD =  8018;
C.HD.AUICD =  8019;
C.HD.ALLCD =  8020;
C.HD.AULLCD = 8021;

% Ecodes that identify the trial parameters
% (Parameters that change on a trial to trial basis
% and so are dropped at the end of each trial)
C.TP.ICD =    7000;
C.TP.DCD =    7001;
C.TP.FCD =    7002;
C.TP.LCD =    7003;
C.TP.CCD =    7004;
C.TP.SCD =    7005;
C.TP.USCD =   7006;
C.TP.ULCD =   7007;
C.TP.UICD =   7008;
C.TP.LLCD =   7009;
C.TP.ULLCD =  7010;
C.TP.AICD =   7011;
C.TP.ADCD =   7012;
C.TP.AFCD =   7013;
C.TP.ALCD =   7014;
C.TP.ACCD =   7015;
C.TP.ASCD =   7016;
C.TP.AUSCD =  7017;
C.TP.AULCD =  7018;
C.TP.AUICD =  7019;
C.TP.ALLCD =  7020;
C.TP.AULLCD = 7021;

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
    'fpon_t',    'time', 0, C.FPONCD
    'fpacq_t',   'time', 0, C.FPACQCD
    'fpoff_t',   'time', 0, C.FPOFFCD
    'stimon_t',  'time', 0, C.STIMONCD
    'stimoff_t', 'time', 0, C.STIMOFFCD
    'targon_t',  'time', 0, C.TARGONCD
    'targoff_t', 'time', 0, C.TARGOFFCD
    'i',		  'int', 0, C.TP.ICD
    'd',       'double', 0, C.TP.DCD
    'f',        'float', 0, C.TP.FCD
    'l',         'long', 0, C.TP.LCD
    'c',         'char', 0, C.TP.CCD
    's',        'short', 0, C.TP.SCD
    'us',      'ushort', 0, C.TP.USCD
    'ul',       'ulong', 0, C.TP.ULCD
    'ui',        'uint', 0, C.TP.UICD
    % a 64-bit integer trialparam will be cast to double in the stro. you may
    % lose precision if your value is > 2^53; use additionalRasterFields
    }';

%for the stro.sum.expParams array (includes things that get dropped only
%once):
exptParamsHeader = {
    'i',	   'int', 0, C.HD.ICD
    'd',    'double', 0, C.HD.DCD
    'f',     'float', 0, C.HD.FCD
    'l',      'long', 0, C.HD.LCD
    'c',      'char', 0, C.HD.CCD
    's',     'short', 0, C.HD.SCD
    'us',   'ushort', 0, C.HD.USCD
    'ul',    'ulong', 0, C.HD.ULCD
    'ui',     'uint', 0, C.HD.UICD
    'll',    'llong', 0, C.HD.LLCD
    'ull',  'ullong', 0, C.HD.ULLCD
    'is',      'int', 1, C.HD.AICD
    'ad',   'double', 1, C.HD.ADCD
    'af',    'float', 1, C.HD.AFCD
    'al',     'long', 1, C.HD.ALCD
    'ac',     'char', 1, C.HD.ACCD
    'as',    'short', 1, C.HD.ASCD
    'aus',  'ushort', 1, C.HD.AUSCD
    'aul',   'ulong', 1, C.HD.AULCD
    'aui',    'uint', 1, C.HD.AUICD
    'all',   'llong', 1, C.HD.ALLCD
    'aull', 'ullong', 1, C.HD.AULLCD
    }';

additionalRasterFields = {
    'll',    'llong', 0, C.TP.LLCD
    'ull',  'ullong', 0, C.TP.ULLCD
    'is',      'int', 1, C.TP.AICD
    'ad',   'double', 1, C.TP.ADCD
    'af',    'float', 1, C.TP.AFCD
    'al',     'long', 1, C.TP.ALCD
    'ac',     'char', 1, C.TP.ACCD
    'as',    'short', 1, C.TP.ASCD
    'aus',  'ushort', 1, C.TP.AUSCD
    'aul',   'ulong', 1, C.TP.AULCD
    'aui',    'uint', 1, C.TP.AUICD
    'all',   'llong', 1, C.TP.ALLCD
    'aull', 'ullong', 1, C.TP.AULLCD
    }';

%define how good and bad trials are to be distinguished
badLogic = 'find(trialCodes(:,2) == C.ABORTCD, 1)';

%define other routines as an array of function handles. the first row of
%the array should be a name of the handle, and the second row should be the
%actual function handle. The functions should be included as non-nested
%subfunctions in this file. DTStimCodes.m serves as a template
otherInstructions = {}; %none for now
