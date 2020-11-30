function CheckExpoVersion(expoDataSet, matlabImportVersion);
%function CheckExpoVersion(expoDataSet, matlabImportVersion);
%
% An Expo helper function used in all of the main Expo Matlab utility
% functions to check version compatibilty.  The dataset object contains two
% version numbers - the Expo version that generated the XML and the version 
% of the Matlab code - MatlabImportVersion - that converted the XML to a Matlab data file.  
% The version numbers are compared with the locally stored ExpoVersion # and the
% MatlabImportVersion # provided as a parameter.
%
% Used in ReadExpoXML, GetEvents, GetPasses, GetEvents, GetSpikeTimes,
% GetWaveformes, GetAnalog, GetDots
%
%   Author:      Julian Brown
%	Version:     1.1
%   Last updated:  2005-03-28
%   E-mail:      julianb@stanford.edu
%
% See also ReadExpoXML, GetSlots, GetPasses, GetEvents, GetSpikeTimes, GetWaveformes, GetAnalog, GetDots

    expoVersion = '1.1';

    if ~isfield(expoDataSet, 'ExpoVersion')
        error('The expoDataSet object you provided has no ExpoVersion field and is therefore invalid.')
    end

    if ~strcmp(expoDataSet.ExpoVersion, expoVersion)
        if sscanf(expoDataSet.ExpoVersion, '%f') > sscanf(expoVersion,'%f')
            error('Incorrect version of expoDataSet object. It was produced using Expo v%s but this routine was built for v%s.',expoDataSet.ExpoVersion, expoVersion)
        else
            warning('The expoDataSet object was generated from version %s of Expo whereas this Matlab code was written for version %s.',expoDataSet.ExpoVersion, expoVersion)
        end
    end

    if ~strcmp(expoDataSet.MatlabImportVersion, matlabImportVersion)
        if sscanf(expoDataSet.ExpoVersion, '%i') > sscanf(matlabImportVersion,'%i')
            error('Incorrect MatlabImportVersion of expoDataSet object (v%s). You must regenerate using ReadExpoXML v%s.',expoDataSet.MatlabImportVersion, matlabImportVersion)
        else
            warning('The expoDataSet object was generated using version %s of ExpoMatlab''s ReadExpo whereas this is ExpoMatlab code is version %s.',expoDataSet.MatlabImportVersion, matlabImportVersion)
        end
    end
    
return
