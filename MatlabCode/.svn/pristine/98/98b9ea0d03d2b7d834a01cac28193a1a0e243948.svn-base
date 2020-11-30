
function [stro,filename,no_subunit] = library()
% This function contains a repository of all the files collected and allows the user to select any 1 file of interest. 
filename = 'NULL';
no_subunit = 0;
% ********************************************************************************
%                            Declaring the library here 
%                               (List of good files)
% ********************************************************************************


%**************************
% Files from  NUT
%**************************
% stro = nex2stro(findfile('N021915003.nex')); % This file has just 1 subunit
% stro = nex2stro(findfile('N022715002.nex')); % simple luminance cell 
% stro = nex2stro(findfile('N040415002.nex')); % This has 2 subunits and the background as the 3rd subinit for modes 4,5,6
% filename = 'N050415001.nex'; stro = nex2stro(findfile(filename)); % excellent isolation, has 1 subunit, luminance cell
% stro = nex2stro(findfile('N041115001.nex')); % good isolation, works fine, low firing neuron, raw ensemble hist not circularly symmetric
% stro = nex2stro(findfile('N050215001-header-merge.nex')); % excellent isolation, works fine, raw ensemble hist not circularly symmetric
% stro = nex2stro(findfile('N050215002.nex')); % excellent isolation, works fine, Greg thinks this could be a single opponent color cell, raw ensemble hist not circularly symmetric
% stro = nex2stro(findfile('N050415002.nex')); % good isolation,works fine, subunit has a center-surround configuration, raw ensemble hist not circularly symmetric
% stro = nex2stro(findfile('N060515001.nex')); % channel 'a' : poor isolation, has structure in PC1, channel 'b': good isolation, also has nice structure 
% stro = nex2stro(findfile('N060515002.nex')); % excellent isolation, has structure in PC1
% stro = nex2stro(findfile('N060615002.nex')); % poor/good isolation, blue-yellow STA. blue-yellow DO cell, excellent data
% stro = nex2stro(findfile('N060615005.nex')); % good isolation, small blue spot in the STA, luminance structure in PC
% stro = nex2stro(findfile('N060815001.nex')); % poor/good isolation, blue spot in STA, 'AND' comp
% stro = nex2stro(findfile('N060815002.nex')); % poor isolation, same as the previous cell
% stro = nex2stro(findfile('N060815003.nex')); % poor isolation, luminance cell, good data
% stro = nex2stro(findfile('N060815006.nex')); % good isolation, could be a Red-Green DO cell, is also very similar to a luminance cell. 
% filename = 'N061515001.nex'; stro = nex2stro(findfile(filename)); % poor isolation, complex luminance cell. good data. STA not very clear.
% stro = nex2stro(findfile('N061615001.nex')); % good isolation, red-green DO STA 
% stro = nex2stro(findfile('N061615002.nex')); % poor isolation. Has structure in PC1
% stro = nex2stro(findfile('N061615003.nex')); % excellent isolation, has PC1 
% stro = nex2stro(findfile('N061715002.nex')); % excellent isolation, STA has a blue spot. Interesting
% stro = nex2stro(findfile('N061715003.nex')); % poor/good isolation, STA, has a luminance signature (accidently lost the plx file)
% stro = nex2stro(findfile('N061715004.nex')); % good isolation, luminance cell, linear cell , (accidently lost the plx file)
% stro = nex2stro(findfile('N062115001.nex')); % poor isolation, complex luminance cell; beautiful data
% stro = nex2stro(findfile('N062315002.nex')); % good isolation, yellow-orange STA, good isolation
% stro = nex2stro(findfile('N062915002.nex')); % excellent isolation, Beautiful complex luminance cell
% stro = nex2stro(findfile('N063015001.nex')); % poor isolation, luminance signature in PC1. raw STA unclear
% stro = nex2stro(findfile('N063015002.nex')); % good isolation, Blue spot in STA, signature in luminance in PC1 and PC2, good data
% stro = nex2stro(findfile('N063015003.nex')); % poor/good isolation, complex luminance cell
% stro = nex2stro(findfile('N070115001.nex')); % good isolation, Blue spot in STA, complex luminance cell, good data
% stro = nex2stro(findfile('N070315002.nex')); % poor isolation, L-M cell, SO
% stro = nex2stro(findfile('N070715001.nex')); % poor isolation, Blue Spot in STA, STA very clear, good data, could be a red-green DO cell
% stro = nex2stro(findfile('N070715005.nex')); % excellent isolation. Could be a direction selective cell.
% stro = nex2stro(findfile('N070815001.nex')); % poor isolation, blue-purple center-surround mask, linear cell, interesting data
% stro = nex2stro(findfile('N070815002.nex')); % poor isolation, luminance cell,interesting STA, contrast changes from yellow to green, interesting data
% stro = nex2stro(findfile('N070815004.nex')); % poor isolation, blue-purple center-surround cell, interesting data
% stro = nex2stro(findfile('N071015002.nex')); % poor/good isolation, beautiful Complex luminance cell with a blue signature in STA. good data
% stro = nex2stro(findfile('N071015003.nex')); % good/ excellent isolation, interesting luminance in STA, PC1, PC2. Good data 
% stro = nex2stro(findfile('N071015004.nex')); % good isolation, complex luminance cell
% stro = nex2stro(findfile('N071215001.nex')); % poor isolation, green spot in the STA. simple luminance cell
% stro = nex2stro(findfile('N071615002.nex')); % good isolation, simple luminance cell which changes contrast polarity,  excellent data, low spatial frequency
% stro = nex2stro(findfile('N071715001.nex')); % excellent isolation, complex luminance cell, PC1 and PC2 are interesting
% stro = nex2stro(findfile('N071815001.nex')); % excellent isolation. less data. only one subunit is playing a role.
% stro = nex2stro(findfile('N071815003.nex')); % poor isolation, center surround mask. blue spot in STA. good data. red-green DO
% stro = nex2stro(findfile('N071815004.nex')); % poor isolation, complex luminance cell.
% stro = nex2stro(findfile('N072215003.nex')); % excellent isolation, more of a complex luminance cell, OKish data
% stro = nex2stro(findfile('N072315002.nex')); % poor isolation, simple luminance cell. OKish data. 
% stro = nex2stro(findfile('N072615002.nex')); % excellent isolation. simple luminance cell. center-surround organisation. good data.
% stro = nex2stro(findfile('N072815001.nex')); % good isolation, orange-yellow smudge in STA. good data
% stro = nex2stro(findfile('N073015002.nex')); % poor isolation. blue smudge in STA. structure in PC1 also. good data
stro = nex2stro(findfile('N073015003.nex')); % excellent isolation. structure in STA. red-green DO cell. good data.
% stro = nex2stro(findfile('N073015004.nex')); % poor isolation. low spatial-frequency simple luminance cell. good data.
% stro = nex2stro(findfile('N080315001.nex')); % poor isolation. red-green DO cell. good data.
% stro = nex2stro(findfile('N080515002.nex')); % poor isolation, good data. looks like red-green DO cell
% stro = nex2stro(findfile('N080515003.nex')); % poor isolation, good data. looks like a red-green DO cell
% stro = nex2stro(findfile('N080715002.nex')); % excellent isolation, good data, STA has a green-magenta DO shade. Greg also has a cell. could be a blue-yellow DO.
% stro = nex2stro(findfile('N080715003.nex')); % good isolation, good data, complex luminance cell oriented at 45 degrees
% stro = nex2stro(findfile('N081015001.nex')); % excellent isolation, good data, simple luminance cell
% stro = nex2stro(findfile('N081115001.nex')); % poor isolation, good data, center-surround luminance cell.
% stro = nex2stro(findfile('N081115002.nex')); % excellent isolation, good data, Luminance cell. Center-surround have same polarity
% stro = nex2stro(findfile('N081215002.nex')); % excellent isolation, good data, complex luminance cell.
% stro = nex2stro(findfile('N081415001.nex')); % excellent isolation, bad data, great simple luminance STA. could not hold the cell long enough for subunit paradigm
% stro = nex2stro(findfile('N081515001.nex')); % good isolation, good data, less data. subunit should have been center surround and not adjacent dark-light bars
% stro = nex2stro(findfile('N081515003.nex')); % excellent isolation. good data Great file. simple luminance cell oriented at 135 degress. adjacent light-dark bars.
% stro = nex2stro(findfile('N081715001.nex')); % good isolation, good data, center-surround luminance cell.
% stro = nex2stro(findfile('N081715002.nex')); % excellent isolation, good data, center-surround luminance cell. structure in STV and PC1.
% stro = nex2stro(findfile('N081815002.nex')); % poor isolation, good data, center-surround color cell. Yet to characterize which kind of color cell is it.
% stro = nex2stro(findfile('N082015002.nex')); % poor/good isolation, good data, direction selective simple luminance cell. very clear STA

%**************************
% Files from  PANGU
%**************************

% Monkey data - neurophysiology (good files: "Neurothresh")
% stro = nex2stro(findfile('P050116001.nex')); % good isolation, simple cell, good WhiteNoise Data, maybe good neurothresh data, .
% stro = nex2stro(findfile('P052016002.nex')); % poor isolation, LGN afferent, center-surround organization, 
% stro = nex2stro(findfile('P052316005.nex')); % excellent isolation, oriented low sf rf, good Neurothresh Data, linear cell
% stro = nex2stro(findfile('P052616001.nex')); % excellent isolation, green subunit, good Neurothresh Data, non-linear cell
% stro = nex2stro(findfile('P052816003.nex')); % excellent isolation,simple cell, Neurothresh Data
% stro = nex2stro(findfile('P052916001.nex')); % excellent isolation, red color subunit cell, Neurothresh Data
% stro = nex2stro(findfile('P072716002.nex')); % excellent isolation, oriented luminance cell, all out of gamut
% stro = nex2stro(findfile('P072816001.nex')); % excellent isolation, has STA,PC1 and PC2. 2 TFR, Super interesting data
% stro = nex2stro(findfile('P072816004.nex')); % excellent isolation, oriented luminance cell. AND operation, good data
% stro = nex2stro(findfile('P080116002.nex')); % excellent isolation, oriented simple luminance cell, 2 TFR
% stro = nex2stro(findfile('P080216003.nex')); % excellent isolation, oriented luminance cell, 2 TFR 
% stro = nex2stro(findfile('P080316001.nex')); % excellent isolation, red STA, 2TFR
% stro = nex2stro(findfile('P080516003.nex')); % good isolation, STA and PC1, 2 TFR, all out of gamut
% stro = nex2stro(findfile('P080716001.nex')); % poor/good isolation, oriented luminance cell, Good WhiteNoise, Neurothresh all out of gamut, 2TFR, high baseline, LGN afferent  
% stro = nex2stro(findfile('P080816001.nex')); % excellent isolation, oriented luminance cell, Good WhiteNoise, Neurothresh all out of gamut 
% stro = nex2stro(findfile('P080816003.nex')); % good isolation, center-surround LGN afferent, Good WhiteNoise, Neurothresh all out of gamut
% stro = nex2stro(findfile('P081216004.nex')); % excellent isolation, STA and PC1, Good WhiteNoise, STA vs PC search, All of of gamut
% stro = nex2stro(findfile('P081616002.nex')); % excellent isolation, blue-yellow STA, 2 TFR, good data


% To be tested
%---------------------------------------------


%**************************************************************************
%                                  3 subunits
%**************************************************************************
% stro = nex2stro(findfile('N050115002.nex')); % working fine, just need to extend the analysis for 3 subunits

%**************************************************************************
%                                Bad/Non-relevant files 
%**************************************************************************
% stro = nex2stro(findfile('N050115001.nex')); % header absent
% stro = nex2stro(findfile('N050215001.nex')); % header absent
% stro = nex2stro(findfile('N041115002.nex'));% has no mask/subunit
% stro = nex2stro(findfile('N041215003.nex')); % really bad data file.
% stro = nex2stro(findfile('N050115001-header-merge.nex')); % has no mask
% stro = nex2stro(findfile('N041215001.nex')); % works fine, junk
% stro = nex2stro(findfile('N041215002.nex')); % works fine, junk
% stro = nex2stro(findfile('N042815001.nex')); % works fine, junk
% stro = nex2stro(findfile('N060615003.nex')); % no mask selected, lost the cell
% stro = nex2stro(findfile('N060615006.nex')); % no mask selected, lost the cell
% stro = nex2stro(findfile('N061715001.nex')); % STA not very clear. not interesting
% stro = nex2stro(findfile('N061815001.nex')); % linear along S2 unit. Nothing special. Would consider a bad data
% stro = nex2stro(findfile('N061815002.nex')); % poor isolation, luminance signature in PC1. bad data
% stro = nex2stro(findfile('N061815005.nex')); % luminance structure in frame 4 STA. Not very clear though.
% stro = nex2stro(findfile('N062315004.nex')); % noisy data. looks like a luminance cell
% stro = nex2stro(findfile('N062915003.nex')); % No prominent structure in PCs and the STA
% stro = nex2stro(findfile('N070715004.nex')); % excellent isolation, looks like a luminance cell. OKish data
% stro = nex2stro(findfile('N071715003.nex')); % poor isolation, blue spot in the STA, not so good data
% stro = nex2stro(findfile('N071915003.nex')); % PC1 not clear. could be a complex luminance cell. not good data. good isolation.
% stro = nex2stro(findfile('N071515001.nex')); % interesting cell, luminance structure in STA and PC1. could be an artifact because of Nut's erratic behavior. bad data
% filename = 'N072215002.nex'; stro = nex2stro(findfile('N072215002.nex')); % excellent isolation blue stripe in the STA, bad data

%% REX testing files
% stro = nex2stro(findfile('N102215001.nex')); 
% stro = nex2stro(findfile('N102315001.nex'));
% stro = nex2stro(findfile('N102915001.nex'));
% stro = nex2stro(findfile('N110115002.nex'));
% stro = nex2stro(findfile('N110115004.nex'));
% stro = nex2stro(findfile('N110115005.nex'));
% stro = nex2stro(findfile('N012616002.nex'));
% stro = nex2stro(findfile('N012616004.nex'));

%% 
% Use the script 'WN_nosubunit.m' to run these scripts.
%**************************************************************************
% Example cells collected by Gregory Horwitz which can be used for references
%**************************************************************************
% stro = nex2stro(findfile('K110408003.nex')); no_subunit = 1;
% stro = nex2stro(findfile('K080808006.nex')); no_subunit = 1;
% stro = nex2stro(findfile('K032708001.nex')); no_subunit = 1;
% stro = nex2stro(findfile('K110508001.nex')); no_subunit = 1;
% stro = nex2stro(findfile('K042809001.nex')); no_subunit = 1;
% stro = nex2stro(findfile('K110608001.nex')); no_subunit = 1; % blue-yellow DO
% stro = nex2stro(findfile('A052215001.nex')); no_subunit = 1; % This cell has green to yellow contrast transition similar to N070815002.nex
% stro = nex2stro(findfile('K043008002.nex')); no_subunit = 1; % red-green DO
% stro = nex2stro(findfile('K020209005.nex')); no_subunit = 1; % red-green DO
% stro = nex2stro(findfile('K021309004.nex')); no_subunit = 1; % red-green DO, STA quite faint, not a great example

%**************************************************************************
%             Example LGN cells collected by Emily Gelfand
%**************************************************************************
% stro = nex2stro(findfile('A052215003.nex')); no_subunit = 1;
% stro = nex2stro(findfile('A061215001.nex')); no_subunit = 1;
% stro = nex2stro(findfile('A062315003.nex')); no_subunit = 1;
% stro = nex2stro(findfile('A062315006.nex')); no_subunit = 1;

%**************************************************************************
%             Example LGN cells collected by Greg
%**************************************************************************
% stro = nex2stro(findfile('A062217004.nex')); no_subunit = 0;
% stro = nex2stro(findfile('A062317001.nex')); no_subunit = 1;
% stro = nex2stro(findfile('A062317003.nex')); no_subunit = 1;
% stro = nex2stro(findfile('A062317005.nex')); no_subunit = 1;
% stro = nex2stro(findfile('A062317006.nex')); no_subunit = 0;
% stro = nex2stro(findfile('A062817006.nex')); no_subunit = 0;
end

