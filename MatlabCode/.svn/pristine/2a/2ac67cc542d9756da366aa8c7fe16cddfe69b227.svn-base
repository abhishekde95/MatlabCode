function units = InitializeUnitConversionMatrices(expoDataSet)
% function [units] = InitializeUnitConversionMatrices(expoDataSet)
%
% An Expo helper function used for unit conversions
%
% Used in GetEvents, GetAnalog, GetSpikeTimes, GetWaveforms etc.
%
%   Author:      Julian Brown
%   Last updated:  2004-12-28
%   E-mail:      julianb@stanford.edu

    tickDur = expoDataSet.environment.Conversion.TickDuration * 1e-6;
    tickFreq = 1/tickDur;
    pixTocm = expoDataSet.environment.Conversion.PixelSize;
    
    analogRangeVolts = expoDataSet.environment.Conversion.AnalogRangeVolt;
    VToNormV = 1/analogRangeVolts;
    
    analogRangeDeg = expoDataSet.environment.Conversion.AnalogRangeDeg;
    
    cmToPix = 1/pixTocm;
    cmToDeg = 180 * atan(1/expoDataSet.environment.Conversion.ViewDistance)/pi;
    cmToNormVolt = cmToDeg/analogRangeDeg;
    cmToVolt = cmToNormVolt * analogRangeVolts;
    
    pixToDeg = pixTocm*cmToDeg;
    pixToNormVolt = pixToDeg/analogRangeDeg;
    pixToVolt = analogRangeVolts * pixToNormVolt;
    
    degToNormVolt = 1/analogRangeDeg;
    degToVolt = degToNormVolt * analogRangeVolts;
    degToRad = pi/180;
    degToPix = 1/pixToDeg;
    degTocm = 1/cmToDeg;
    
    sizeOfConversionMatrix = 29;
    
    % unit enumeration from Expo unit.h
    U_ANY = 0; U_NONE = 1; U_TICK = 2; U_SEC = 3; U_PIX = 4; U_CM = 5; U_DEGD = 6; U_DEGOR = 7; U_RAD = 8; U_PER_TICK = 9;
    U_PER_SEC = 10; U_C_SEC = 11; U_PER_PIX = 12; U_PER_CM = 13; U_PER_DEG = 14; U_C_DEG = 15; U_C_CM = 16; U_VOLT = 17; U_NORMVOLT = 18;
    U_IMP_SEC = 19; U_IMP_BIN = 20; U_IMP = 21; U_PER_DEG_L = 22; U_C_DEG_L = 23; U_DEG_SEC = 24; U_DEG_TICK =25; U_PIX_SEC = 26; U_PIX_TICK = 27;
    U_BASETIME = 28; U_MSEC=29;
    
    units.U_NORMVOLT = U_NORMVOLT;
    units.U_BASETIME = U_BASETIME;
    
    % set an array of units that can be used to identify each group of units
    units.GroupRepUnits = [U_TICK U_PIX U_DEGOR U_PER_TICK U_PER_PIX U_PER_DEG_L U_DEG_SEC];
    
    % create sparse matrices representing conversions factors s and whether the conversions require a recipricol r
    i = [U_TICK  U_PIX   U_PIX    U_PIX     U_PIX         U_CM    U_CM     U_CM         U_DEGD    U_DEGD        U_VOLT     U_DEGOR  U_PER_TICK U_PER_TICK U_PER_SEC]; 
    j = [U_SEC   U_CM    U_DEGD   U_VOLT    U_NORMVOLT    U_DEGD  U_VOLT   U_NORMVOLT   U_VOLT    U_NORMVOLT    U_NORMVOLT U_RAD    U_PER_SEC  U_C_SEC    U_C_SEC]; 
    s = [tickDur pixTocm pixToDeg pixToVolt pixToNormVolt cmToDeg cmToVolt cmToNormVolt degToVolt degToNormVolt VToNormV   degToRad tickDur    tickDur   1 ]; 
    r = [0       0       0        0         0             0       0        0            0         0             0          0        0          1          1 ];

    i = [i, [U_PER_PIX U_PER_PIX U_PER_PIX U_PER_PIX U_PER_CM  U_PER_CM U_PER_CM U_PER_DEG U_PER_DEG U_C_DEG ]];
    j = [j, [U_PER_CM  U_PER_DEG U_C_DEG   U_C_CM    U_PER_DEG U_C_DEG  U_C_CM   U_C_DEG   U_C_CM    U_C_CM  ]];
    s = [s, [pixTocm   pixToDeg  pixToDeg  cmToPix   cmToDeg   cmToDeg  1        1         cmToDeg   cmToDeg ]];
    r = [r, [0         0         1         1         0         1        1        1         1         0       ]];
    
    i = [i, [U_PER_DEG_L U_DEG_SEC  U_DEG_SEC U_DEG_SEC        U_DEG_TICK        U_DEG_TICK U_PIX_SEC  U_TICK        U_TICK       U_SEC      U_SEC  U_MSEC]];
    j = [j, [U_C_DEG_L   U_DEG_TICK U_PIX_SEC U_PIX_TICK       U_PIX_SEC         U_PIX_TICK U_PIX_TICK U_BASETIME    U_MSEC       U_BASETIME U_MSEC U_BASETIME]];
    s = [s, [1           tickDur    degToPix  degToPix*tickDur degToPix*tickFreq degToPix   tickDur    tickDur*10000 tickDur*1000 10000      1000   10]];
    r = [r, [1           0          0         0                0                 0          0          0             0            0          0      0]];
    
    halfConversionMatrix = sparse(i, j, s, sizeOfConversionMatrix, sizeOfConversionMatrix);
    
    M = full(halfConversionMatrix);
    
    % to complete the matrix we need the opposite conversions eg U_SEC -> U_TICK 
    % get transpose 
    MT = M';
    
    % we need to get recipricol of all non-zero values
    % to avoid division by zero subtract 1 from all 0 elements 
    MT = MT - (MT==0);
    
    %get the recipricol of each element
    MT = 1./MT;
    
    %now add 1 to all of the elements that are -1 to revert them back to 0
    MT = MT + (MT==-1);
    
    % add the processed transpose and add the identity matrix (to give identity converions such as U_TICK -> U_TICK) to the original matrix
    units.ConversionMatrix = sparse(M + MT + eye(sizeOfConversionMatrix, sizeOfConversionMatrix));
    
    % create the matrix that specifies which conversions require a recipricol value e.g. U_PER_SEC -> U_C_SEC 
    halfRecipricolMatrix = sparse(i, j, r, sizeOfConversionMatrix, sizeOfConversionMatrix);
    
    R = full(halfRecipricolMatrix);
    
    units.ReciprocalMatrix = sparse(R + R');
    
    units.Names = {'ticks', 'sec', 'pix', 'cm', 'deg', 'deg', 'rad', 'period ticks', 'period sec', ...
                'cyc/sec', 'period pix', 'period cm', 'period deg', 'cyc/deg', 'cyc/cm', 'volts', ...
                'norm volts', 'impulses/sec',  'impulses/bin', 'impulses', 'period deg', 'cyc/deg', ...
                'deg/sec', 'deg/tick', 'pix/sec', 'pix/tick', '1/10msec', 'msec' };

return
