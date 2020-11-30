function conversionFactor = GetConversionFactor(expoDataSet, fromUnit, toUnit)
%function conversionFactor = GetConversionFactor(expoDataSet, fromUnit, toUnit)
%
% An Expo function used for calculating the conversion fact between
% different units.
%
% Available unit types and groups are:
%               ticks                 sec            1/10msec                msec
%                 pix                  cm                 deg               volts          norm volts
%                 deg                 rad
%        period ticks          period sec             cyc/sec
%          period pix           period cm          period deg             cyc/deg              cyc/cm
%          period deg             cyc/deg
%             deg/sec            deg/tick             pix/sec            pix/tick
%
%
% See also ReadExpoXML, GetSlots, GetPasses, GetEvents, GetSpikeTimes, GetAnalog, 
% GetPSTH, PlotPSTH, GetWaveforms, GetStartTimes, GetEndTimes, GetDuration,
%
%   Author:      Julian Brown
%   Last updated:  2004-12-28
%   E-mail:      julianb@stanford.edu

units = expoDataSet.environment.Conversion.units;

fromUnitNum = GetUnitNumFromName(units, fromUnit);
toUnitNum = GetUnitNumFromName(units, toUnit);

conversionFactor = ConvertNum(units, 1, fromUnitNum, toUnitNum);

return
