
% Author - Abhishek De
% WNthresh_ST function for analysing files obtained by running the WhiteNoiseOnlinethresh_AD.m paradigm
% 05/2016 - Made some changes. Added the flexibility to choose the mex files to compute STC. use_STCOVmex_st
% is variable that decides which file to use.

clearvars;
close all;
echo off;
plot_counter = 1; % defining a counter to keep a track on the index of the figure
% ********************************************************************************
% REX testing
% ********************************************************************************
% stro = nex2stro(findfile('N020516001.nex'));
% stro = nex2stro(findfile('N020716001.nex'));
% stro = nex2stro(findfile('N020816001.nex'));
% stro = nex2stro(findfile('N020916001.nex'));
% stro = nex2stro(findfile('N020916002.nex'));
% stro = nex2stro(findfile('N020916003.nex'));
% stro = nex2stro(findfile('N020916004.nex'));
% stro = nex2stro(findfile('N020916006.nex'));
% stro = nex2stro(findfile('N020916007.nex'));
% stro = nex2stro(findfile('N021116001.nex'));
% stro = nex2stro(findfile('N021716002.nex'));
% stro = nex2stro(findfile('N022216001.nex')); % test file collected for ~ 27 minutes, 705 trials, 64 directions
% stro = nex2stro(findfile('N022316001.nex')); % test file collected for ~ 39 minutes, 1048 trials, 96 directions
% stro = nex2stro(findfile('N022416001.nex')); % test file ~ 35 minutes, "subunits" are basis vectors, 979 trials, 120 directions
% stro = nex2stro(findfile('N022516001.nex')); % test file ~ 9 minutes, "subunits", 2 targetspikerates, 229 trials
% stro = nex2stro(findfile('N022516002.nex')); % test file ~ 16 minutes, "subunits", 2 targetspikerates, 522 trials
% stro = nex2stro(findfile('N022516003.nex')); % test file ~ 23 minutes , "STAvsPC", 3 targetspikerates, 653 trials
% stro = nex2stro(findfile('N030716002.nex')); % test file to analyze the bug the basis vec was not getting passed into
% the file while the monkey makes a saccade out of the fixation window, resulting in an ABORTCD
% stro = nex2stro(findfile('N030716003.nex')); % same as above
% stro = nex2stro(findfile('N030716006.nex')); % same as above
% stro = nex2stro(findfile('N030716007.nex')); % same as above
% stro = nex2stro(findfile('N040116001.nex'));
% stro = nex2stro(findfile('N040116002.nex'));
% stro = nex2stro(findfile('N040116003.nex'));
% stro = nex2stro(findfile('N040116004.nex'));
% stro = nex2stro(findfile('N040116005.nex'));
% stro = nex2stro(findfile('N040116006.nex'));
% stro = nex2stro(findfile('N040116007.nex'));
% stro = nex2stro(findfile('N040216002.nex'));
% stro = nex2stro(findfile('N040216003.nex'));
% stro = nex2stro(findfile('N040216004.nex'));
% stro = nex2stro(findfile('N040216005.nex'));
% stro = nex2stro(findfile('N040216006.nex'));
% stro = nex2stro(findfile('N040216007.nex'));
% stro = nex2stro(findfile('N040216008.nex'));
% stro = nex2stro(findfile('N040316001.nex'));
% stro = nex2stro(findfile('N040316002.nex'));
% stro = nex2stro(findfile('N040316003.nex'));
% stro = nex2stro(findfile('N040316004.nex'));
% stro = nex2stro(findfile('N040316005.nex'));
% stro = nex2stro(findfile('N040316006.nex'));
% stro = nex2stro(findfile('N040316007.nex'));
% stro = nex2stro(findfile('N040316008.nex'));
% stro = nex2stro(findfile('N040416001.nex')); % subunits
% stro = nex2stro(findfile('N040416002.nex')); % STA vs PC1
% stro = nex2stro(findfile('N040416003.nex')); % subunits
% stro = nex2stro(findfile('N040616001.nex')); % subunits
% stro = nex2stro(findfile('N041116001.nex')); % subunits
% stro = nex2stro(findfile('N041116002.nex')); % subunits
% stro = nex2stro(findfile('N041116003.nex')); % subunits
% stro = nex2stro(findfile('N041216001.nex')); % subunits
% stro = nex2stro(findfile('N041316001.nex')); % subunits
% stro = nex2stro(findfile('N041316002.nex')); % subunits
% stro = nex2stro(findfile('N041316004.nex')); % subunits
% stro = nex2stro(findfile('N041316005.nex')); % subunits
% stro = nex2stro(findfile('N041316006.nex')); % subunits, strange error noticed, switched to a non-neurothresh trial but droppped a correct code
% stro = nex2stro(findfile('N041416001.nex'));
% stro = nex2stro(findfile('N041416002.nex')); % subunits, with REVERSALCD implemented
% stro = nex2stro(findfile('N041516001.nex')); % subunits, with REVERSALCD implemented
% stro = nex2stro(findfile('P041516002.nex'));
% stro = nex2stro(findfile('P041516003.nex'));
% stro = nex2stro(findfile('P041516004.nex'));
% stro = nex2stro(findfile('P041616001.nex'));
% stro = nex2stro(findfile('P041616002.nex'));
% stro = nex2stro(findfile('P041616003.nex'));
% stro = nex2stro(findfile('N041716001.nex'));
% stro = nex2stro(findfile('N041716002.nex'));
% stro = nex2stro(findfile('N041916001.nex')); % subunits, introduced a PARENTVERTICESCD, forgot to bookend TP_PARENTVERTICES in REX code
% stro = nex2stro(findfile('N041916002.nex'));
% stro = nex2stro(findfile('P042216001.nex')); % subunits, reversalflagrequest in PLEXON will drop the actual number of reversals instead of just 0 or 1.
% stro = nex2stro(findfile('P042216002.nex')); % Online and offline latency don't match up
% stro = nex2stro(findfile('P042216003.nex'));
% stro = nex2stro(findfile('P072916001.nex')); % Implementing search in first quadrant
% stro = nex2stro(findfile('M121016001.nex')); % Implementing angle bisection instead of line bisection
% stro = nex2stro(findfile('M121016002.nex'));
% stro = nex2stro(findfile('M121316001.nex')); % Angle bisection along with search along 2 dirs in one round
% stro = nex2stro(findfile('M121316002.nex'));
% stro = nex2stro(findfile('M121316003.nex'));
% stro = nex2stro(findfile('M121316005.nex')); % Something's wrong in this file
% stro = nex2stro(findfile('M121316006.nex')); % Still something wrong, hopefully have figured out next
% stro = nex2stro(findfile('M121316007.nex')); % Error fixed, STA vs PC1, 1 TFR
% stro = nex2stro(findfile('M121316008.nex')); % subunit, 1 TFR
% stro = nex2stro(findfile('M121416001.nex')); % subunit, 2 TFR
% stro = nex2stro(findfile('M022117001.nex')); % 3 subunits, Neurothresh, 1TFR, bug in the code
% stro = nex2stro(findfile('M022117002.nex')); % 3 subunits, Neurothresh, 1TFR, code fixed
% stro = nex2stro(findfile('M022117003.nex')); % Checking if code works with 2 subunits, Neurothresh, 1TFR,
% stro = nex2stro(findfile('M101517001.nex')); % file with new calibration structure
% stro = nex2stro(findfile('Junk082618001.nex')); % file with some more additions - RF structure and cone/gun weights

% ********************************************************************************
% Monkey Data - behavior
% ********************************************************************************
% stro = nex2stro(findfile('M030316001.nex')); % test file ~ 75 minutes, "STAvsPC", 1 targetspikerate, 1079 trials, 128 weights
% stro = nex2stro(findfile('P030716006.nex')); % test file, "STAvsPC", 1 targetspikerate, 434 trials, 32 weights
% stro = nex2stro(findfile('P030816001.nex')); % test file, "subunits", 1 targetspikerate, 510 trials,
% stro = nex2stro(findfile('P030916001.nex')); % test file, "STA vs PC", 1 targetspikerate, 1228 trials
% stro = nex2stro(findfile('P031116001.nex')); % test file, "STA vs PC", 1 targetspikerate, 1238 trials
% stro = nex2stro(findfile('P031616001.nex')); % test file, "STA vs PC", 1 targetspikerate, 1361 trials

% ********************************************************************************
% Monkey data - neurophysiology (bad files)
% ********************************************************************************
% stro = nex2stro(findfile('P031816002.nex')); % poor isolation, green center in STA, luminance cell, probed in 4 directions
% stro = nex2stro(findfile('P032216001.nex')); % poor isolation, center- surround green STA, luminance cell, probed in 8 directions
% stro = nex2stro(findfile('P032216002.nex')); % poor isolation, center-surround blue-green STA, Neurothresh broke down
% stro = nex2stro(findfile('P033116002.nex')); % good isolation, simple cell, Neurothresh broke down
% stro = nex2stro(findfile('P040516002_2.nex')); % good isolation,
% stro = nex2stro(findfile('P040516003.nex')); % good isolation, lost the cell
% stro = nex2stro(findfile('P040516004.nex')); % good isolation,
% stro = nex2stro(findfile('P041016001.nex')); % multi-unit
% stro = nex2stro(findfile('P042016001.nex')); % noisy WhiteNoise STA, great subunit STA, red-green opponent signal, probed completely in first 8 direction, next 8 directions weren't probed completely
% stro = nex2stro(findfile('P042016002.nex')); % noisy WhiteNoise and subunit STA, descent isolation
% stro = nex2stro(findfile('P042016003.nex')); % noisy WhiteNoise STA, great subunit STA, red-green opponent signal
% stro = nex2stro(findfile('P042716001.nex')); % simple cell, poor isolation,
% stro = nex2stro(findfile('P050216001.nex')); % excellent isolation, no clear structure in STA, might have been a complex cell (look STV)
% stro = nex2stro(findfile('P050616001.nex')); % good isolation, no clear structure in STA
% stro = nex2stro(findfile('P052316001.nex')); % poor isolation, could be an LGN afferent, green oriented structure in the STA, good WhiteNoise Data
% stro = nex2stro(findfile('P052316003.nex')); % excellent isolation, green spot in the STA, no subunit data
% stro = nex2stro(findfile('P070416002.nex')); % excellent isolation, nothing in STA, interesting Neurothresh Data
% stro = nex2stro(findfile('P080416004.nex')); % poor isolation, luminance complex cell, Just the WhiteNoise
% stro = nex2stro(findfile('P072616004.nex')); % good isolation, seems like a luminance cell (not sure), something is weird about this data
% stro = nex2stro(findfile('M122016002.nex')); % excellent isolation, bad file
% stro = nex2stro(findfile('M050517002.nex')); % excellent isolation, bad data set, Neurothresh completed for 1 TFR, STAvsPC1
% stro = nex2stro(findfile('M060917001.nex')); % excellent isolation, no STA, Bad WhiteNoise Data
% stro = nex2stro(findfile('M092617003.nex')); % excellent isolation, bad neurothresh file
% stro = nex2stro(findfile('P092717002.nex')); % excellent isolation, center surround luminance cell, Not Good WhiteNoise subunit data

% ********************************************************************************
% Monkey data - neurophysiology (good files: "WhiteNoise Checkerboard")
% ********************************************************************************
% stro = nex2stro(findfile('M113016005.nex')); % excellent isolation, complex cell, Good WhiteNoise Checkerboard Data
% stro = nex2stro(findfile('M120716002.nex')); % good isolation, good WhiteNoise checkerboard data, oriented luminance cell
% stro = nex2stro(findfile('M121816003.nex')); % poor/good isolation, even symmetric red-cyan cell, good WhiteNoise checkerboard data
% stro = nex2stro(findfile('M121916002.nex')); % excellent isolation, luminance cell, Good WhiteNoise checkerboard data
% stro = nex2stro(findfile('M122416001.nex')); % poor/good isolation, Good WhiteNoise checkerboard Data
% stro = nex2stro(findfile('M122616002.nex')); % excellent isolation, Good WhiteNoise checkerboard Data
% stro = nex2stro(findfile('M122516003.nex')); % excellent isolation, Good WhiteNoise Checkboard Data.
% stro = nex2stro(findfile('M122516001.nex')); % good isolation, complex cell, Good WhiteNoise checkerboard data
% stro = nex2stro(findfile('M010117001.nex')); % poor/good isolation, complex cell, Good WhiteNoise checkerboard Data
% stro = nex2stro(findfile('M010317001.nex')); % good isolation, yellowish-orange tinge in STA, Good WhiteNoise Checkerboard Data
% stro = nex2stro(findfile('M010317002.nex')); % good isolation, yellowish-orange tinge in STA, Good WhiteNoise Checkerboard Data
% stro = nex2stro(findfile('M011417001.nex')); % good isolation, simple luminance cell, Good WhiteNoise checkerboard Data
% stro = nex2stro(findfile('M011517001.nex')); % good isolation, blue-yellow cell, Good WhiteNoise Checkerboard Data
% stro = nex2stro(findfile('M013117002.nex')); % good isolation, luminance cell, Good WhiteNoise checkerboard data
% stro = nex2stro(findfile('M021617002.nex')); % good isolation, oriented luminance cell, Good WhiteNoise checkerboard Data
% stro = nex2stro(findfile('M021717002.nex')); % excellent isolation, luminance cell, Good WhiteNoise Checkerboard Data
% stro = nex2stro(findfile('M022017001.nex')); % excellent isolation, vertically oriented luminance cell, Good WhiteNoise check data
% stro = nex2stro(findfile('P080216002.nex')); % good isolation, Good Whitenoise file without subunit, blue STA
% stro = nex2stro(findfile('P082216002.nex')); % excellent isolation, red-green DO cell, Only WhiteNoise Data (not enough spikes from subunit)
% stro = nex2stro(findfile('M070217004.nex')); % good isolation, good whitenoise check data
% stro = nex2stro(findfile('P090617004.nex')); % good isolation, BY color cell, Good WhiteNoise checkerboard data


% ********************************************************************************
% Monkey data - neurophysiology (good files: "WhiteNoise Subunit")
% ********************************************************************************
% stro = nex2stro(findfile('P032916001.nex')); % excellent isolation, simple cell, Neurothresh broke down (didn't send the 3rd and the 4th weights at all, 1 targetspikerate, "Subunits", weights never made it into the file
% stro = nex2stro(findfile('P033116001.nex')); % good isolation, simple cell, probed in 8 directions
% stro = nex2stro(findfile('P042816002.nex')); % simple cell, good WhiteNoise Data,
% stro = nex2stro(findfile('P051916001.nex')); % good isolation, color cell (orange)
% stro = nex2stro(findfile('P051616005.nex')); % excellent isolation, complex cell
% stro = nex2stro(findfile('P052316004.nex')); % excellent isolation, green spot in the STA, relevant WhiteNoise Data, Pangu shook while the Neurothresh button was pressed
% stro = nex2stro(findfile('P052816001.nex')); % good isolation, single opponent red-green cell with a single subunit. Could potentially be a good WhiteNoise Data
% stro = nex2stro(findfile('P072816002.nex')); % Poor/good isolation, oriented simple luminance cell, good WhiteNoise Data
% stro = nex2stro(findfile('P080216001.nex')); % good isolation, Good WhiteNoise file, blue STA
% stro = nex2stro(findfile('P080616001.nex')); % good isolation, Good WhiteNoise Data, oriented luminance cell
% stro = nex2stro(findfile('M113016006.nex')); % excellent isolation, high sf luminance cell, Good WhiteNoise Data, No Neurothresh
% stro = nex2stro(findfile('P080816003.nex')); % good isolation, center-surround LGN afferent, Good WhiteNoise, Neurothresh all out of gamut
% stro = nex2stro(findfile('P081016004.nex')); % good isolation, luminance cell, Good WhiteNoise, no neurothresh
% stro = nex2stro(findfile('P081116001.nex')); % excellent isolation, luminance cell, Good WhiteNoise, no neurothresh
% stro = nex2stro(findfile('P081716001.nex')); % poor/good isolation, oriented luminance cell, STA and PC1, incomplete search Neurothresh
% stro = nex2stro(findfile('M120916002.nex')); % good isolation, blue-yellow center surround cell, Good WhiteNoise Data
% stro = nex2stro(findfile('M121716002.nex')); % good isolation, oriented complex cell, good WhiteNoise Data,
% stro = nex2stro(findfile('M121716003.nex')); % good isolation, oriented complex cell, good WhiteNoise Data
% stro = nex2stro(findfile('M121816002.nex')); % good isolation, oriented luminance cell, Good WhiteNoise Data
% stro = nex2stro(findfile('M121816006.nex')); % excellent isolation, complex cell, Good WhiteNoise Data
% stro = nex2stro(findfile('M121916003.nex')); % excellent isolation, oriented luminance cell, Good WhiteNoise Data
% stro = nex2stro(findfile('M122316001.nex')); % poor/good isolation, oriented luminance complex cell, Good WhiteNoise Data
% stro = nex2stro(findfile('M122416005.nex')); % excellent isolation, 2 units, 1st unit is the good one, Good WhiteNoise Data
% stro = nex2stro(findfile('M122516004.nex')); % excellent isolation, STA and PC1, Good WhiteNoise Data, only 4 points in Neurothresh, STA vs PC1
% stro = nex2stro(findfile('M010117003.nex')); % good isolation, luminance simple cell, low SF, Good WhiteNoise subunit Data
% stro = nex2stro(findfile('P081516002.nex')); % excellent isolation, excellent WhiteNoise data, center-surround luminance cell, Neurothresh: completed search in 1 dir only
% stro = nex2stro(findfile('M010617002.nex')); % good/excellent isolation, has STA, can be used as a Good WhiteNoise File, Neurothresh data few points, STA vs PC1
% stro = nex2stro(findfile('M012417001.nex')); % good isolation, 2 units, luminance cells, Good WhiteNoise Data
% stro = nex2stro(findfile('M020517003.nex')); % excellent isolation, luminance STA, Good WhiteNoise Data
% stro = nex2stro(findfile('M020617002.nex')); % excellent isolation, luminance cell, center-surround, Good WhiteNoise Data
% stro = nex2stro(findfile('M021017002.nex')); % excellent isolation, complex luminance cell, Good WhiteNoise data
% stro = nex2stro(findfile('M021117001.nex')); % excellent isolation, center-surround luminance cell, Good WhiteNoise Data
% stro = nex2stro(findfile('M021217001.nex')); % excellent isolation, STA and PC1, Good WhiteNoise Data
% stro = nex2stro(findfile('M021217005.nex')); % excellent isolation, STA and PC1, Good WhiteNoise Data
% stro = nex2stro(findfile('M022217002.nex')); % excellent isolation, luminance complex cell, Good WhiteNoise Data
% stro = nex2stro(findfile('M022317001.nex')); % excellent isolation, luminance complex cell, Good WhiteNoise Data, Over 15000 spikes, classic textbook example cell
% stro = nex2stro(findfile('M022417001.nex')); % excellent isolation, luminance complex cell, Good WhiteNoise Data
% stro = nex2stro(findfile('M022817002.nex')); % excellent isolation, has STA and PC1, Good WhiteNoise data
% stro = nex2stro(findfile('M030417001.nex')); % excellent isolation, color opponent cell, Good WhiteNoise Data
% stro = nex2stro(findfile('M040417001.nex')); % excellent isolation, luminance cell, Good WhiteNoise Data
% stro = nex2stro(findfile('M040917001.nex')); % excellent isolation, complex cell, Good WhiteNoise Data, Less Neurothresh data
% stro = nex2stro(findfile('M051917001.nex')); % excellent isolation, luminance STA and PC1 , Good WhiteNoise subunit Data
% stro = nex2stro(findfile('M060817001.nex')); % excellent isolation, luminance STA, Good WhiteNoise Data
% stro = nex2stro(findfile('M061617001.nex')); % excellent isolation, BY opponent STA, Good WhiteNoise Data
% stro = nex2stro(findfile('M062217004.nex')); % excellent isolation, center surround luminance simple cell, Good WhiteNoise Data
% stro = nex2stro(findfile('M062417001.nex')); % excellent isolation, blue -yellow cell, Good WhiteNoise Data
% stro = nex2stro(findfile('M062417002.nex')); % excellent isolation, blue center surround cell, Good WhiteNoise Data
% stro = nex2stro(findfile('P082216001.nex')); % poor isolation, orange single opponent cell, only WhiteNoise data
% stro = nex2stro(findfile('P052016002.nex')); % poor isolation, LGN afferent, center-surround organization, not good neurothresh data
% stro = nex2stro(findfile('P052816003.nex')); % excellent isolation,simple cell, Neurothresh Data, all out of gamut
% stro = nex2stro(findfile('P052916001.nex')); % excellent isolation, red color subunit cell, Neurothresh Data, all out of gamut
% stro = nex2stro(findfile('P072716002.nex')); % excellent isolation, oriented luminance cell, all out of gamut
% stro = nex2stro(findfile('P080116002.nex')); % excellent isolation, oriented simple luminance cell, 2 TFR, something wrong with this file, basis vec was never passed to the plx file
% stro = nex2stro(findfile('P080516003.nex')); % good isolation, STA and PC1, 2 TFR, all out of gamut
% stro = nex2stro(findfile('P080716001.nex')); % poor/good isolation, oriented luminance cell, Good WhiteNoise, Neurothresh all out of gamut, 2TFR, high baseline, LGN afferent
% stro = nex2stro(findfile('P080816001.nex')); % excellent isolation, oriented luminance cell, Good WhiteNoise, Neurothresh all out of gamut
% stro = nex2stro(findfile('P081216004.nex')); % excellent isolation, STA and PC1, Good WhiteNoise, STA vs PC search, All of of gamut
% stro = nex2stro(findfile('P081716002.nex')); % poor/good isolation, STA and PC1, STA vs PC Neurothresh, 1 TFR, not good Neurothresh Data
% stro = nex2stro(findfile('M062917003.nex')); % excellent isolation, oriented luminance simple cell, Good WhiteNoise subunit Data
% stro = nex2stro(findfile('M070117001.nex')); % Excellent isolation, oriented color cell, potential BY cell, Good WhiteNoise Data
% stro = nex2stro(findfile('M070117002.nex')); % Good isolation, center surround (color?) cell, Good WhiteNoise Data
% stro = nex2stro(findfile('P070417001.nex')); % good isolation, luminance simple cell, Good WhiteNoise Data
% stro = nex2stro(findfile('M080217002.nex')); % excellent isolation, center-surround luminance cell, Good WhiteNoise subunit data
% stro = nex2stro(findfile('P080817003.nex')); % poor/good isolation, center-surround luminance cell, Good WN data
% stro = nex2stro(findfile('M082217003.nex')); % excellent isoltion, luminance simple cell, Good WhiteNoise subunit Data
% stro = nex2stro(findfile('M082517003.nex')); % excellent isolation, center-surround color cell, Good WhiteNoise subunit data
% stro = nex2stro(findfile('P090217002.nex')); % good isolation, single-opponent color cell, Good WhiteNoise subunit data
% stro = nex2stro(findfile('M090617001.nex')); % good/excellent isolation, color cell, Good WhiteNoise subunit data
% stro = nex2stro(findfile('P090617005.nex')); % good isolation, DO color cell, Good WhiteNoise subunit data
% stro = nex2stro(findfile('M090717001.nex')); % good isolation, color cell, Good WhiteNoise subunit data
% stro = nex2stro(findfile('P091217003.nex')); % excellent isolation, luminance s11imple cell, Good WhiteNoise subunit data
% stro = nex2stro(findfile('M092217007.nex')); % excellent isolation, luminance simple cell, Good WhiteNoise subunit data
% stro = nex2stro(findfile('P092217003.nex')); % excellent isolation, color cell, Good WhiteNoise subunit Data
% stro = nex2stro(findfile('M092617004.nex')); % excellent isolation, luminance simple cell (could be direction selective), good WhiteNoise data
% stro = nex2stro(findfile('M092817001.nex')); % excellent isoltion, luminance complex cell, Good Whitenoise subunit data
% stro = nex2stro(findfile('M100717003.nex')); % excellent isolation, luminance simple cell, Good WhiteNoise subunit Data
% stro = nex2stro(findfile('M101717001.nex')); % good isolation, multi-unit, single opponent color cell, Good Whitenoise subunit data
% stro = nex2stro(findfile('P102517002.nex')); % excellent isolation, BY color cell, Good Whitenoise subunit data
% stro = nex2stro(findfile('P102517005.nex')); % excellent isolation, BY color cell, Good Whitenoise subunit data
% stro = nex2stro(findfile('P103017006.nex')); % good isolation, center-surround luminance cell, Good Whitenoise data
% stro = nex2stro(findfile('M110617001.nex')); % excellent isolation, color cell has STA and PC1, Good Whitenoise subunit data
% stro = nex2stro(findfile('M111317002.nex')); % excellent isolation, single opponent color cell, Good WhiteNoise Data
% stro = nex2stro(findfile('P111317003.nex')); % excellent isolation, single opponent color cell, has STA and PC1, Good WhiteNoise Data
% stro = nex2stro(findfile('M111517005.nex')); % excellent isolation, luminance simple cell, Good WhiteNoise data
% stro = nex2stro(findfile('P122017001.nex')); % good isolation, luminance simple cell, Good WhiteNoise Data,  2 subunits + surround
% stro = nex2stro(findfile('P122817002.nex')); % excellent isolation, luminance simple cell, Good WhiteNoise data, 2 subunits + surround
% stro = nex2stro(findfile('P103118002.nex')); % excellent isolation, lime-magenta DO cell, Good WhiteNoise data,

% ********************************************************************************
% Monkey data - neurophysiology (good files: "Neurothresh")
% ********************************************************************************
% stro = nex2stro(findfile('P050116001.nex')); % good isolation, simple cell, good WhiteNoise Data, maybe good neurothresh data, .
% stro = nex2stro(findfile('P052316005.nex')); % excellent isolation, oriented low sf rf, good Neurothresh Data, linear cell
% stro = nex2stro(findfile('P052616001.nex')); % excellent isolation, green subunit, good Neurothresh Data, non-linear cell
% stro = nex2stro(findfile('P072816001.nex')); % excellent isolation, has STA,PC1 and PC2. 2 TFR, Super interesting data
% stro = nex2stro(findfile('P072816004.nex')); % excellent isolation, oriented luminance cell. AND operation, good data
% stro = nex2stro(findfile('P080216003.nex')); % excellent isolation, oriented luminance cell, 2 TFR, good file, non-linear cell
% stro = nex2stro(findfile('P080316001.nex')); % excellent isolation, red STA, 2TFR, non-linear cell
% stro = nex2stro(findfile('P081616002.nex')); % excellent isolation, blue-yellow STA, 2 TFR, good data
% stro = nex2stro(findfile('P081716003.nex')); % excellent isolation, yellow-pink cell, 2TFR, subunit Neurothresh, good data
% stro = nex2stro(findfile('P082616002.nex')); % excellent isolation, excellent WhiteNoise data, Neurothresh: completed in 2 directions
% stro = nex2stro(findfile('M120116002.nex')); % excellent isolation, oriented luminance cell, Neurothresh: completed for 2 target firing rates, STA vs PC1
% stro = nex2stro(findfile('M120216001.nex')); % excellent isolation, has center surround luminance STA and PC1, Neurothresh: completed for 2 target firing rates, subunits
% stro = nex2stro(findfile('M120216002.nex')); % excellent isolation, has oriented luminance complex edge, Neurothresh: completed for 2 target firing rates, STA vs PC1
% stro = nex2stro(findfile('M120916003.nex')); % excellent isolation, blue-yellow center surround cell, Neurothresh: completed for 2 target firing rates, subunit
% stro = nex2stro(findfile('M121216001.nex')); % excellent isolation, oriented luminance cell, Neurothresh: completed for 2 target firing rates, subunit
% stro = nex2stro(findfile('M121616004.nex')); % excellent isolation, center-surround Red-Green cell, has STA and PC1,Neurothresh: completed for 2 target firing rates, subunit
% stro = nex2stro(findfile('M121916001.nex')); % excellent isolation, oriented luminance cell, STA, PC1 and PC2, Neurothresh completed for 1 TFR
% stro = nex2stro(findfile('M122016003.nex')); % excellent isolation, red-cyan edge detector, STA, Neurothresh completed for 1 TFR, linear cell
% stro = nex2stro(findfile('M122116001.nex')); % excellent isolation, red-cyan center surround cell, has STA and PC1, Neurothresh completed for 2 TFRs, subunit
% stro = nex2stro(findfile('M122216001.nex')); % excellent isolation, complex cell, Neurothresh completed for 2 TFRs, STA vs PC1
% stro = nex2stro(findfile('M122316002.nex')); % excellent isolation, oriented luminance cell, Neurothresh completed for 2 TFRs, subunit
% stro = nex2stro(findfile('M122416004.nex')); % excellent isolation, blue-yellow center surround cell, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('M122616001.nex')); % excellent isolation, luminance simple cell, Neurothresh completed for 2 TFR, subunit, linear cell
% stro = nex2stro(findfile('M122716008.nex')); % excellent isolation, oriented complex cell, Neurothresh completed for 1 TFR, STA vs PC1
% stro = nex2stro(findfile('M122816002.nex')); % good isolation, sig001b, complex cell, has STA and PC1, Neurothresh completed for 2TFRs, STA vs PC1
% stro = nex2stro(findfile('M122816003.nex')); % excellent isolation, complex cell, Neurothresh completed for 1 TFR, STA vs PC1
% stro = nex2stro(findfile('M122916001.nex')); % good isolation, simple cell, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('M122916002.nex')); % good isolation, complex cell, Neurothresh completed for 1 TFR, STA vs PC1
% stro = nex2stro(findfile('M123116004.nex')); % excellent isolation, color cell, linear, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('M010117004.nex')); % excellent isolation, simple cell, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('M010217001.nex')); % excellent isolation, red-green color cell, Neurothresh completed for 1 TFR, subunit
stro = nex2stro(findfile('M010317003.nex')); % excellent isolation, blue - yellow cell, linear cell, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('M010317004.nex')); % excellent isolation, STA, color PC1, Neurothresh completed for 1 TFR, less data points, subunit
% stro = nex2stro(findfile('M010417004.nex')); % excellent isolation, STA, PC1, complex cell, Neurothresh completed for 1 TFR, STA vs PC1
% stro = nex2stro(findfile('M010617003.nex')); % excellent isolation, has STA and lum PC, Neurothresh completed for 1 TFR, STA vs PC1
% stro = nex2stro(findfile('M010917001.nex')); % excellent isolation, luminance simple cell, linear cell, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('M011217001.nex')); % poor/good isolation, luminance STA and PC1, Neurothresh completed for 1 TFR, STA vs PC1
% stro = nex2stro(findfile('M011217002.nex')); % excellent isolation, luminance simple cell, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('M011417002.nex')); % good/excellent isolation, simple luminance cell, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('M011417003.nex')); % excellent isolation, simple luminance cell, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('M011517002.nex')); % excellent isolation, complex luminance cell, Neurothresh completed for 1 TFR, STAvsPC1
% stro = nex2stro(findfile('M011717001.nex')); % excellent isolation, luminance cell, STA and PC1, Neurothresh completed for 2 TFRs, STAvsPC1
% stro = nex2stro(findfile('M011817008.nex')); % excellent isolation, luminance simple cell, STA, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('M012217002.nex')); % excellent isolation, luminance cell, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('M012417007.nex')); % excellent isolation, sig100a multi unit, oriented complex cell, Neurothresh for 1 TFR, STA vs PC1
% stro = nex2stro(findfile('M012517001.nex')); % excellent isolation, complex cell, Neurothresh completed for 1 TFR, STA vs PC1
% stro = nex2stro(findfile('M013017001.nex')); % good isolation, red-cyan color cell, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('M013117001.nex')); % excellent isolation, kinda luminance/color, longer latency, STA and PC1, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('M020117001.nex')); % excellent isolation, STA and PC1, STA not very clear, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('M020317005.nex')); % excellent isolation, luminance STA, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('M020517001.nex')); % excellent isolation, has luminance STA and color PC1, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('M020617001.nex')); % excellent isolation, color STA & luminance PC1, Neurothresh completed for 1 TFR, STAvsPC1
% stro = nex2stro(findfile('M021317001.nex')); % excellent isolation, color cell, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('M021417004.nex')); % excellent isolation, simple luminance cell, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('M021617003.nex')); % excellent isolation, oriented luminance cell, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('M021717001.nex')); % excellent isolation, oriented luminance cell, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('M021817001.nex')); % excellent isolation, color cell, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('M021817003.nex')); % excellent isolation, has STA and PC1, Neurothresh completed for 1 TFR, STAvsPC1
% stro = nex2stro(findfile('M022217003.nex')); % excellent isolation, luminance complex cell, Neurothresh completed for 1 TFR, STAvsPC1
% stro = nex2stro(findfile('M022317002.nex')); % excellent isolation, luminance complex cell, Neurothresh completed for 1 TFR, STAvsPC1, really good data
% stro = nex2stro(findfile('M022717003.nex')); % excellent isolation, luminance cell, has STA and PC1, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('M022817001.nex')); % excellent isolation, luminance cell, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('M030117001.nex')); % excellent isolation, center-surround luminance cell, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('M030717003.nex')); % excellent isolation, multi-unit, luminance complex cell, Neurothresh completed for 1 TFR, STAvsPC1
% stro = nex2stro(findfile('M040217001.nex')); % excellent isolation, has STA and PC1, Neurothresh completed for 1 TFR, STAvsPC1
% stro = nex2stro(findfile('M040717001.nex')); % excellent isolation, complex cell, Neurothresh completed for 1 TFR, STAvsPC1
% stro = nex2stro(findfile('M041117001.nex')); % excellent isolation, looks like single opponent, has STA and PC1, Neurothresh completed for 1 TFR, STAvsPC1
% stro = nex2stro(findfile('M041917001.nex')); % excellent isolation, blue-yellow color cell, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('M050317001.nex')); % excellent isolation, one subunit, luminance cell, Neurothresh completed for 2 TFR, subunit
% stro = nex2stro(findfile('M050517001.nex')); % excellent isolation, luminance simple cell, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('M050617003.nex')); % excellent isolation, center-surround color cell, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('M050617004.nex')); % excellent isolation, low SF luminance cell, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('M050717001.nex')); % excellent isolation, center surround cell, Neurothresh, completed for 1 TFR, subunit
% stro = nex2stro(findfile('M050717002.nex')); % excellent isolation, center surround cell, Neurothresh completed for 1 TFR, subunit, noisy data, Good WhiteNoise Data
% stro = nex2stro(findfile('M051917003.nex')); % excellent isolation, center surround cell, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('M061217001.nex')); % excellent isolation, center surround luminance cell, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('M061317005.nex')); % good/excellent isolation, simple luminance cell, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('M061417001.nex')); % excellent isolation, oriented luminance simple cell, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('M061517001.nex')); % excellent isolation, red-green DO cell, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('M012417002.nex')); % excellent isolation, sig100a, luminance cell, Neurothresh completed for 1 TFR, STA vs PC1, descent Neurothresh Data
% stro = nex2stro(findfile('M070117003.nex')); % Excellent isolation, oriented luminance simple cell, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('M070217003.nex')); % excellent isolation,oriented luminance simple cell, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('M073117002.nex')); % excellent isolation, oriented luminance simple cell, has STA and PC1, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('M080217003.nex')); % excellent isolation, luminance simple cell, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('M080817002.nex')); % excellent isolation, center-surround color cell, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('M081417003.nex')); % excellent isolation, center-surround color cell, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('M082117002.nex')); % excellent isolation, luminance complex cell, Neurothresh completed for 1 TFR, STAvsPC1
% stro = nex2stro(findfile('M082217001.nex')); % excellent isolation, luminance complex cell, Neurothresh completed for 1 TFR, STAvsPC1
% stro = nex2stro(findfile('M082517004.nex')); % excellent isolation, oriented luminance simple cell, Neurothresh completed for 1 TFR, subunit, linear cell
% stro = nex2stro(findfile('M082917001.nex')); % excellent isolation, center-surround color cell, Neurothresh completed for 1 TFR, subunit, few points
% stro = nex2stro(findfile('P082917001.nex')); % excellent isolation, luminance simple cell, Neurothresh completed for 1 TFR, subunit, maybe nonlinear
% stro = nex2stro(findfile('M090117002.nex')); % good isolation, DO cell, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('P090117001.nex')); % excellent isolation, color cell, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('P090117002.nex')); % excellent isolation, color cell, Neurothresh completed for 1 TFR, subunit, non-linear cell
% stro = nex2stro(findfile('P090217001.nex')); % excellent isolation, color cell, Neurothresh completed for 1 TFR, subunit, less points
% stro = nex2stro(findfile('M090617002.nex')); % good/excellent isolation, color cell, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('M090717002.nex')); % excellent isolation, color cell, has PC1, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('M090817001.nex')); % good/excellent isolation, BY color cell, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('M091117002.nex')); % excellent isolation, color cell, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('M091817003.nex')); % excellent isolation, luminance simple cell, Neurothresh completed for 1 TFR, subunit, few points
% stro = nex2stro(findfile('P091817001.nex')); % excellent isolation, has STA and PC1, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('M092017005.nex')); % excellent isolation, kinda luminancy but need to double check, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('M092117001.nex')); % excellent isolation, high sf luminance cell, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('M092717003.nex')); % excellent isolation, luminance simple cell, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('M092917001.nex')); % excellent isolation, luminance simple cell, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('M092917002.nex')); % excellent isolation, luminance simple cell, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('M100417001.nex')); % excellent isolation, BY-DO cell, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('M101017002.nex')); % excellent isolation, color cell, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('M101617001.nex')); % good isolation, BY DO cell, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('M101617002.nex')); % excellent isolation, BY DO cell, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('M101717002.nex')); % good isolation, maybe multi-unit, DO color cell, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('P101717001.nex')); % good isolation, DO cell, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('M101817002.nex')); % good/excellent isolation, DO cell, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('M102017002.nex')); % excellent isolation, luminance simple cell, multi-unit, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('P102017001.nex')); % good/excellent isolation, color cell has PC1, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('P102017004.nex')); % excellent isolation, DO cell, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('M102517001.nex')); % excellent isolation, color cell, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('P102617002.nex')); % excellent isolation, RG DO cell, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('P102617006.nex')); % excellent isolation, color cell, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('P102717003.nex')); % excellent isolation, DO cell, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('P103017001.nex')); % excellent isolation, lum cell/color cell ?, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('P103017002.nex')); % good isolation, color cell, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('M110317001.nex')); % excellent isolation, color cell, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('P110717002.nex')); % excellent isolation, luminance simple cell, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('P111017002.nex')); % excellent isolation, luminance simple cell, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('M111417002.nex')); % excellent isolation, maybe color cell?, Neurothresh completed for 2 TFR, subunit, less points
% stro = nex2stro(findfile('M111517006.nex')); % excellent isolation, luminance simple cell, Neurothresh completed for 1 TFR
% stro = nex2stro(findfile('P112717002.nex')); % excellent isolation, color cell, Neurothresh completed for 1 TFR, subunit, linear cell
% stro = nex2stro(findfile('P112817002.nex')); % excellent isolation, luminance simple cell, 2 subunits + surround, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('P112817003.nex')); % excellent isolation, luminance simple cell, 2 subunits + surround, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('P120517001.nex')); % excellent isolation, nonlinear cell, maybe color cell, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('M120617003.nex')); % excellent isolation, complex cell, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('P121117006.nex')); % good isolation, color cell, Neurothresh completed for 1 TFR, subunit, nonlinear cell
% stro = nex2stro(findfile('P121117007.nex')); % excellent isolation,  DO cell, Neurothresh completed for 1 TFR, subunit, linear cell
% stro = nex2stro(findfile('P121217002.nex')); % excellent isolation, maybe DO cell, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('M121317002.nex')); % good isolation, color cell, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('M121317003.nex')); % excellent isolation, luminance simple cell,Neurothresh completed for 1 TFR,subunit
% stro = nex2stro(findfile('P121317003.nex')); % excellent isolation, DO cell, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('P121317004.nex')); % good isolation, color cell, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('M121417002.nex')); % excellent isolation, luminance simple cell, Neurothresh completed for 1 TFR, subunit, linear cell
% stro = nex2stro(findfile('P121417003.nex')); % excellent isolation, DO cell, Neurothresh completed for 1 TFR, nonlinear cell
% stro = nex2stro(findfile('P121917006.nex')); % good isolation, luminance simple cell, Neurothresh completed for 1 TFR, subunit
% stro = nex2stro(findfile('P091018005.nex')); % excellent isolation, DO cell, Neurothresh completed for 1 TFR (less points), subunit
% stro = nex2stro(findfile('P091418002.nex')); % excellent isolation, blue-yellow DO cell, Neurothresh completed for 1 TFR , subunit
% stro = nex2stro(findfile('P102918001.nex')); % excellent isolation, color cell, Neurothresh completed for 1 TFR , subunit
% stro = nex2stro(findfile('P102918002.nex')); % excellent isolation, DO cell, Neurothresh completed for 1 TFR , subunit
% stro = nex2stro(findfile('P110118003.nex')); % excellent isolation, color cell, Neurothresh completed for 1 TFR, subunit


% ********************************************************************************
% Analysis pending (Maui) Opto dlx Whitenoise files
% ********************************************************************************
% stro = nex2stro(findfile('M032518003.nex'));
% stro = nex2stro(findfile('M032518006.nex'));
% stro = nex2stro(findfile('M032518009.nex'));
% stro = nex2stro(findfile('M032518010.nex'));

% ********************************************************************************
%                Declaring the local and global variables here
% ********************************************************************************
% Declaring global variables
global spikename maskidx spikeidx neurothreshidx nstixperside ngammasteps seedidx nframesidx correctidx
global fponidx stimoffidx stimonidx muidxs sigmaidxs basisvecidx weightsidx fpacqidx targetspikerateidx basisvecdiridx latencyidx reversalflagidx parentverticesidx
global msperframe ntrials maxT xx yy M linepredtol stepsizescale stepsize nreversals oogscale
spikename = getSpikenum(stro);
maskidx = strcmp(stro.sum.rasterCells(1,:),'subunit_mask');
spikeidx = strcmp(stro.sum.rasterCells(1,:),spikename);
basisvecidx = strcmp(stro.sum.rasterCells(1,:),'basis_vec');
weightsidx = strcmp(stro.sum.rasterCells(1,:),'weights');
parentverticesidx = strcmp(stro.sum.rasterCells(1,:),'parentvertices');
nstixperside = stro.sum.exptParams.nstixperside;
ngammasteps = 2^16; % 65536
linepredtol = stro.sum.exptParams.linepredtol;
stepsizescale = stro.sum.exptParams.stepsizescale;
stepsize = stro.sum.exptParams.stepsize;
nreversals = stro.sum.exptParams.nreversals;
oogscale = stro.sum.exptParams.oogscale;
seedidx = strcmp(stro.sum.trialFields(1,:),'seed');
nframesidx = strcmp(stro.sum.trialFields(1,:),'num_frames');
stimonidx = strcmp(stro.sum.trialFields(1,:),'stim_on');
stimoffidx = strcmp(stro.sum.trialFields(1,:),'stim_off');
fponidx = strcmp(stro.sum.trialFields(1,:),'fp_on');
fpacqidx = strcmp(stro.sum.trialFields(1,:),'fpacq');
basisvecdiridx = strcmp(stro.sum.trialFields(1,:),'weights_idx');
neurothreshidx = strcmp(stro.sum.trialFields(1,:),'neurothresh'); % when exactly the neurothresh trials started
targetspikerateidx = strcmp(stro.sum.trialFields(1,:),'targetspikerate');
correctidx = strcmp(stro.sum.trialFields(1,:),'correct');
muidxs = [find(strcmp(stro.sum.trialFields(1,:),'mu1')), ...
    find(strcmp(stro.sum.trialFields(1,:),'mu2')), ...
    find(strcmp(stro.sum.trialFields(1,:),'mu3'))];
sigmaidxs = [find(strcmp(stro.sum.trialFields(1,:),'sigma1')), ...
    find(strcmp(stro.sum.trialFields(1,:),'sigma2')), ...
    find(strcmp(stro.sum.trialFields(1,:),'sigma3'))];
latencyidx = strcmp(stro.sum.trialFields(1,:),'latency');
reversalflagidx = strcmp(stro.sum.trialFields(1,:),'reversalflag');
msperframe = 1000/stro.sum.exptParams.framerate;
ntrials = size(stro.trial,1);
maxT = 15; % this represents the temporal part in the spatiotemporal receptive field
xx = linspace(stro.sum.exptParams.gauss_locut/1000, stro.sum.exptParams.gauss_hicut/1000,ngammasteps); % xx represents the probabilities. For more info, have a look at the MATLAB 'norminv' function.
yy = norminv(xx'); % defining norminv to extract the values for which the cdf values range between gauss_locut and gauss_hicut

% Obtaining the M matrix, code extracted from Greg, fitting a cubic spline
% using the command 'spline'. 'SplineRaw' only availabe through
% psychtoolbox which I currently don't have now.
fundamentals = stro.sum.exptParams.fundamentals; % CONE FUNDAMENTALS: L,M,S
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]); %1st column - L, 2nd- M, 3rd- S
mon_spd = stro.sum.exptParams.mon_spd; % MONITOR SPECTRAL DISTRIBUTION IN R,G,B
mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
mon_spd = spline([380:4:780], mon_spd', [380:5:780]); % fitting a cubic spline
M = fundamentals'*mon_spd'; % matrix that converts RGB phosphor intensites to L,M,S cone fundamentals
M = inv(M');

mask_changes = [2];
all_masks = stro.ras(:,maskidx);
Fx = @(xi) any(isnan(xi)); % function that finds 'NaN' in a cell array
inds = find(cellfun(Fx,stro.ras(:,basisvecidx))==0);
if isempty(inds)
    inds = size(stro.trial,1)-1;
end
last_wntrial =  inds(1)-1;
for k = 3:last_wntrial
    if isequal(all_masks{k}, all_masks{k-1}) %|| all(all_masks{k} == 0) && any(isnan(all_masks{k-1}))
        continue
    else
        mask_changes = [mask_changes k-1 k]; %#ok<AGROW>
    end
end
if mask_changes(end) == last_wntrial
    mask_changes(end) = [];
else
    mask_changes = [mask_changes  last_wntrial];
end
mask_changes = reshape(mask_changes , 2, []);

% Plotting the results before the subunit selection, white noise pixel stimulus
[plot_counter,basis_vec] = WNpixelplot(stro,plot_counter,mask_changes(:,1),'WNthresh_ST'); % Uses STCOVmex for calculating STC

% ******************************************************************************************
%  This is the meat of the analysis for generating the spike triggered covariance matrix.
%  Modes 1,2 and 3 code perform the analysis only for the the subunits. The modes 4,5 and 6
%  perform the same analysis but consider background as an additional subunit.
% ******************************************************************************************
use_STCOVmex_ST = 1; % 1 - Will use STCOVmex_ST, 0 - Will use STCOVmex
mode_steps = [1 2];
% 1 - calculate the STA and the PCs
% 2 - projects the subunit stimuli onto the STAs and the PCs
mask_changes_idx = [2 2];
% 1 - white noise pixel stimulus
% 2 - subunit white noise stimulus w/o background

prompt = 'Project along multiple directions?(Y/N)';
str = input(prompt,'s');
if ((strcmp('Y',str) || strcmp('y',str)))
    flag = 1;
else
    flag = 0;
end
% 0 - just project onto 'return_vec'
% 1 - project onto 'new_vec' (PC1 and STA)

for mode = mode_steps
    % Only takes trials which have subunits selected on them. The selected trials can be accessed in 'mask_changes'
    trial_span = mask_changes(:,mask_changes_idx(mode));
    for mask_span = trial_span
        if  (mod(mode,3) == 1)
            % Just store enough space to accomodate the 9 frames of the basis vector
            st_mask = stro.ras{trial_span(1),maskidx}; % subunit mask
            st_mask(st_mask == 0) = Inf;
            [stIdxs,~,~] = unique(st_mask); % now the Infs map to nsubunits+1
            num_subunits = length(stIdxs)-any(isinf(stIdxs)); % nsubunits, like subunits A and B
            if use_STCOVmex_ST
                STCOV_st('init', {num_subunits 3 maxT});
            else
                STCOVmex('init', {num_subunits 3 maxT});
            end
            
        elseif (mod(mode,3) == 2)
            % define basis vectors so that the stimuli can be projected onto it
            initargs = define_basis_vec(nrandnums_perchannel,stro,basis_vec,new_vec,flag);
            STPROJmod('init',initargs); % initialising the STPROJmod
        elseif (mod(mode,3) == 0)
            STPROJmod('init',initargs); % initialising the STPROJmod
        end
        
        for k = mask_span(1):mask_span(2)
            nframes = stro.trial(k,nframesidx);
            if (nframes == 0)
                continue;
            end
            seed = stro.trial(k,seedidx);
            mu = stro.trial(k,muidxs)/1000;
            sigma = stro.trial(k,sigmaidxs)/1000;
            
            % org_mask tells u if u have updated the mask or not. If org_mask is non-zero it means at this particular trial
            % u have selected the subunits and need to analyse its computation
            org_mask = stro.ras{k,maskidx};
            if any(org_mask)
                org_mask(org_mask == 0) = Inf;
                [subunitIdxs,~,mask] = unique(org_mask); % now the Infs map to nsubunits+1
                nrandnums_perchannel = length(subunitIdxs)-any(isinf(subunitIdxs)); % nsubunits, like subunits A and B
                mask = [mask; mask+max(mask); mask+2*max(mask)]; %#ok<AGROW>
            else
                nrandnums_perchannel = nstixperside^2; % In this case, it is the 100 pixels which are flickering when no subunits are selected
            end
            
            % assuming Gaussian gun noise only, random number generator
            % routine as a mexfile (getEJrandnums.mexw64)
            invnormcdf = bsxfun(@plus, bsxfun(@times, yy, sigma), mu);
            randnums = getEJrandnums(3*nrandnums_perchannel*nframes, seed); % random numbers for 300 x 9 pixels
            % This is the extracted colors for subunits/pixels using the seed number
            randnums = reshape(randnums, 3*nrandnums_perchannel, nframes);
            for gun = 1:3
                idxs = (1:nrandnums_perchannel)+nrandnums_perchannel*(gun-1);
                randnums(idxs,:) = reshape(invnormcdf(randnums(idxs,:)+1,gun),[length(idxs) nframes]);
            end
            
            rgbs = randnums;
            t_stimon = stro.trial(k, stimonidx);
            spiketimes = (stro.ras{k,spikeidx}-t_stimon)*1000; % observing spiketimes in milliseconds
            frametimes = linspace(0, nframes*msperframe, nframes)+(msperframe/2)';
            % Deleting the spikes taking place before the first 9 frames as I need to look at the 9 preceding frames
            spiketimes(spiketimes < maxT*msperframe) = [];
            % Deleting the spikes that take place after the stimulus was
            % removed as it would imply that these spikes do not result from
            % the stimulus shown on the screen
            spiketimes(spiketimes > frametimes(end)) = [];
            n = hist(spiketimes, frametimes);
            if (mod(mode,3) == 1)
                % need to make some modifications here
                if use_STCOVmex_ST
                    STCOV_st(rgbs(:),n);
                else
                    STCOVmex(rgbs(:),n);
                end
            else
                STPROJmod(rgbs(:),n);
            end
        end
        if (mod(mode,3) == 1)
            % Plot the STA and PCs if mode = 1
            if use_STCOVmex_ST
                [ch,new_vec, basis_vec, plot_counter,subunits,sub_rgb_idx] = compute_spatio_temp_basis_vec(plot_counter,nrandnums_perchannel,mask,flag,mode);
            else
                [basis_vec, plot_counter,subunits,sub_rgb_idx] = compute_basis_vec(plot_counter,nrandnums_perchannel,mask);
                ch = []; new_vec = []; flag = 0;
            end
        elseif (mod(mode,3) == 2)
            % Calculate the projection values of the stimuli on the white noise subunit basis vectors
            [plot_counter,nbins,nbins1] = compute_projection_val_thresh(plot_counter,mode,subunits,1,1,flag,ch);
        end
    end
end

% Estimating the latency from the White Noise or non-Neurothresh trials - new way of calculating t_offset
tmp_basis = basis_vec;
frametimes = [1:maxT]*msperframe;
basisvec_energy =  sum(tmp_basis.^2);
t_offset1 = frametimes(basisvec_energy == max(basisvec_energy(1,1:10))); % Time of the peak energy in milliseconds from 1 to 10 frames so that it matches the offline t_offset
centroid_t = (basisvec_energy/sum(basisvec_energy))*frametimes';
figure(plot_counter);
plot(frametimes,basisvec_energy,'ko-','Linewidth',2); xlabel('Frame Times'); ylabel('Energy'); title('Estimation of latency');
plot_counter = plot_counter + 1;
t_offset1 = t_offset1/1000; % Converting it into seconds
t_offset = stro.trial(end,latencyidx)/1000;
disp([ 'Offline: ' num2str(t_offset1*1000) ' ,Online: ' num2str(stro.trial(end,latencyidx))]);

% Determining when Neurothresh mode was active, plotting the basis vector, working correctly
neurothreshmode = stro.trial(:,neurothreshidx);
basisvec_dropidx = inds(end);
neurothresh_startidx = max(find(neurothreshmode==1,1),basisvec_dropidx+1); % There is possibility of pressing the Neurothresh button multiple times
num_targetspikerates = unique(stro.trial(neurothresh_startidx:end,targetspikerateidx));
vect = stro.ras{basisvec_dropidx,basisvecidx};
basisvec_size = nstixperside*nstixperside*3;
numvect = (numel(vect)/basisvec_size)-1;
basisvec = cell(1,numvect);
figure(plot_counter);
% Actual basis vec
bkgnd_monitor = reshape(vect(numvect*basisvec_size+1:basisvec_size*(numvect+1)),[nstixperside nstixperside 3]);
for ii = 1:numvect
    tmp_vec = vect((ii-1)*basisvec_size+1:basisvec_size*ii) + vect(numvect*basisvec_size+1:basisvec_size*(numvect+1));
    basisvec{ii} = reshape(tmp_vec,[nstixperside nstixperside 3]);
    subplot(1,numvect,ii),image(basisvec{ii}); set(gca,'Xtick',[],'YTick',[]); hold on;
end
hold off;
plot_counter = plot_counter + 1;
tmpmask = reshape(mask,[nstixperside nstixperside 3]);
[s1x,s1y] = find(tmpmask(:,:,1)==1,1);
[s2x,s2y] = find(tmpmask(:,:,1)==2,1);
vec1 = basisvec{1}; vec2 = basisvec{2};
S1RGB = squeeze(vec1(s1x,s1y,:) - bkgnd_monitor(s1x,s1y,:));
S2RGB = squeeze(vec2(s2x,s2y,:) - bkgnd_monitor(s2x,s2y,:));
clear vec1 vec2

% plotting the normalized basis vector
figure(plot_counter);
norm_basisvec = cell(1,numvect);
for ii = 1:numvect
    tmp_vec = basisvec{ii}-bkgnd_monitor;
    normfactor = 0.5/((max(abs(tmp_vec(:))))*1.05);
    norm_basisvec{ii} = normfactor*(tmp_vec)+0.5;
    subplot(1,numvect,ii),image(norm_basisvec{ii}); set(gca,'Xtick',[],'YTick',[]); hold on;
end
hold off;
plot_counter = plot_counter + 1;

% plot the images along 4 basis directions and diagonals to get a sense what do they look like
coeffs = [1 0; 1 1; 0 1;-1 1; -1 0; -1 -1; 0 -1; 1 -1];
tmp_vec1 = basisvec{1}-bkgnd_monitor;
tmp_vec2 = basisvec{2}-bkgnd_monitor;
figure(plot_counter);
for ii = 1:size(coeffs,1)
    tmp_vec = coeffs(ii,1)*tmp_vec1 + coeffs(ii,2)*tmp_vec2;
    normfactor = 0.5/((max(abs(tmp_vec(:))))*1.05);
    tmp_vec = normfactor*(tmp_vec)+0.5;
    subplot(2,4,ii),image(tmp_vec); set(gca,'Xtick',[],'YTick',[]); hold on;
end
hold off;
plot_counter = plot_counter + 1;

% This section works correctly
% plotting the weights corresponding to the weight directions - useful for
% checking if the direction indexes and the directions are aligned or not
weight_direction_mat = [];
for i=neurothresh_startidx:size(stro.trial,1)
    weight_direction_mat = [weight_direction_mat; stro.trial(i,basisvecdiridx) stro.ras{i,weightsidx}'/norm(stro.ras{i,weightsidx}')];
end
[~,idx] = sort(weight_direction_mat(:,1)); % sorting it according to the increasing order of weight direction indexes
weight_direction_mat1 = weight_direction_mat(idx,:);

% plotting the rasterplots for neurothresh trials
norms = cell(1,numel(num_targetspikerates));
completed_search_alongdir = cell(1,numel(num_targetspikerates));
for jj = 1: numel(num_targetspikerates)
    idxs = find(~isnan(stro.trial(:,correctidx)) & stro.trial(:,targetspikerateidx)==num_targetspikerates(jj));
    idxs(idxs<=neurothresh_startidx) = [];
    different_weights = unique(stro.trial(idxs,basisvecdiridx));
    tmp_norm = [];
    tmp_completed_search_alongdir = [];
    
    for kk = 1:numel(different_weights)
        idxs1 = find(stro.trial(:,basisvecdiridx) == different_weights(kk));
        idxs1(idxs1<neurothresh_startidx) = [];
        raster_data = stro.ras(idxs1,1);
        tmp_norm = [tmp_norm; stro.ras{idxs1(end),weightsidx}'];
        for ii = 1:size(raster_data,1)
            tmp = raster_data{ii} ;
            spikes = tmp(tmp>stro.trial(idxs1(ii),stimonidx)+t_offset & tmp < stro.trial(idxs1(ii),stimoffidx));
            spikes = spikes - stro.trial(idxs1(ii),stimonidx)-t_offset;
        end
        [flag, gamutViolation] = Check_ifcompletely_probed(stro,basisvec,bkgnd_monitor,idxs1,tmp_norm(end,:));
        % flag = 0, incompletely probed
        % flag = 1, completely probed
        % gamutViolation = 1, out of gamut point
        tmp_completed_search_alongdir = [tmp_completed_search_alongdir; flag gamutViolation];
    end
    norms{jj} = tmp_norm;
    completed_search_alongdir{jj} = tmp_completed_search_alongdir;
end

% Refer NTpreprocess written by Greg,
% Plotting the end norm or contrast values for each search direction and Converting the end norms into polar coordinates
% Need to write a small function to check for Gamut Violation based on the
% reversalflagidx

grid on;
color = ['g-';'k-';'r-';'b-';'m-'];
lo = -6.0; hi = 6.0;
for ii = 1:size(norms,2)
    tmp = norms{ii};
    completed_dir = completed_search_alongdir{ii};
    probed_dirs = logical(completed_dir(:,1)==1); % only including the directions that have been completely probed
    oog_idx = find(completed_dir(:,1)==1 & completed_dir(:,2)==1); % probed and out of gamut
    not_oog_idx = find(completed_dir(:,1)==1 & completed_dir(:,2)==0);
    fact = 0.5./sqrt(tmp(probed_dirs,1).^2 + tmp(probed_dirs,2).^2); % factor needed to extract unit vector
    figure(plot_counter), subplot(size(norms,2),2,2*ii-1);h = compass(fact.*tmp(probed_dirs,1),fact.*tmp(probed_dirs,2),color(ii,:)); title(['Target FR=' num2str(num_targetspikerates(ii))]);
    set(h,'LineWidth',1);
    [THETA,RHO] = cart2pol(tmp(:,1),tmp(:,2));
    ind = (1:numel(THETA))';
    r = fliplr(linspace(0,1,numel(ind)));
    b = fliplr(r);
    THETA = THETA * (180/pi);
    %Earlier points in time are blue in color and later points in time are red in color
    %     for jj = 1: numel(ind)
    %         figure(plot_counter), subplot(size(norms,2),2,2*ii); plot(tmp(ind(jj),1), tmp(ind(jj),2),'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[r(jj) 0 b(jj)]); hold on;
    %         figure(plot_counter+1), plot(THETA(ind(jj)), log(RHO(ind(jj))),'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[r(jj) 0 b(jj)],'MarkerEdgeColor',[r(jj) 0 b(jj)]); hold on;
    %     end
    figure(plot_counter+1), subplot(size(norms,2),2,2*ii-1); plot(tmp(not_oog_idx,1), tmp(not_oog_idx,2),'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor','b','MarkerEdgeColor','b'); hold on;
    set(gca,'Xlim',[lo,hi],'Ylim',[lo,hi]); colormap(gca,winter); grid on; axis square; xlabel('Basisvec 1'); ylabel('Basisvec 2');
    figure(plot_counter+2),subplot(size(norms,2),2,2*ii-1); plot(THETA(not_oog_idx), log(RHO(not_oog_idx)),'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor','b'); hold on;
    grid on; axis square; xlabel('Theta'); ylabel('Log R');
    if ~isempty(oog_idx)
        figure(plot_counter+1),plot(tmp(oog_idx,1), tmp(oog_idx,2),'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor','g');hold off;
        figure(plot_counter + 2), plot(THETA(oog_idx), log(RHO(oog_idx)),'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor','g'); set(gca,'Xlim',[-180,180]); hold off;
    else
        figure(plot_counter+1), hold off;
        figure(plot_counter+2), hold off;
    end
end
plot_counter = plot_counter + 3;

% Estimating the baseline firing rate for non-Neurothresh trials
% The number of spikes between the fpon and stimon codes
figure(plot_counter);
for i = 1:2
    % First plot the raster plots of the WhiteNoise trials and then the Neurothresh trials
    idxs = find(neurothreshmode==(i-1));
    raster_data = stro.ras(idxs,1);
    num_spikes =[];
    num_dur = [];
    spikerate = [];
    add_time = 1;
    stimofftime = [];
    for ii = 1:size(raster_data,1)
        tmp = raster_data{ii} ;
        spikes_disp = tmp(tmp<stro.trial(idxs(ii),stimonidx)+ add_time & tmp>stro.trial(idxs(ii),fpacqidx));
        spikes = tmp(tmp<stro.trial(idxs(ii),stimonidx) & tmp>stro.trial(idxs(ii),fpacqidx));
        subplot(2,3,1+3*(i-1));plot(spikes_disp-stro.trial(idxs(ii),fpacqidx),(ii-1)*ones(1,length(spikes_disp)),'k.'); hold on;
        num_spikes = [num_spikes; numel(spikes)];
        num_dur = [num_dur; (stro.trial(idxs(ii),stimonidx)- stro.trial(idxs(ii),fpacqidx))];
        stimofftime = [stimofftime; (stro.trial(idxs(ii),stimoffidx)- stro.trial(idxs(ii),fpacqidx))];
        firing_rate = numel(spikes)/(stro.trial(idxs(ii),stimonidx)- stro.trial(idxs(ii),fpacqidx)); % add time not added as I need am only counting spikes between fpacqidx and stimonidx to obtain the firing rate
        spikerate = [spikerate; firing_rate];
    end
    spikerate = num_spikes./num_dur; set(gca,'Xlim',[0-0.05 max(num_dur)+ add_time])
    xlabel('Time in s'); ylabel('Trials'); title('Baseline Rasterplots');
    line([0 0],[0 ii-1],'Color',[1 0 0])% Time when fpacq code was dropped
    line([min(num_dur) min(num_dur)],[0 ii-1],'Color',[1 0 0]) % Time when stimon code was droppped
    line([min(num_dur)+t_offset min(num_dur)+t_offset],[0 ii-1],'Color',[1 0 0]) % Time when stimon code was droppped
    line([mean(stimofftime) mean(stimofftime)],[0 ii-1],'Color',[1 0 0]) % Time when stimoff code was droppped
    hold off;
    figure(plot_counter),subplot(2,3,2+3*(i-1)), hist(stro.trial(idxs,stimonidx) - stro.trial(idxs,fpacqidx)); xlabel('Time'); ylabel('Trials'); title('Time durations');
    figure(plot_counter),subplot(2,3,3+3*(i-1)), hist(spikerate); xlabel('Spikerates in spikes/s'); ylabel('Trials'); title('Spikerate Histogram');
    BL_firingrate = mean(spikerate);
end
plot_counter = plot_counter + 1;

% Code to see if there are any OFF responses in the rasters from the Neurothresh trials
% The number of spikes between the stimoff and stimon codes
% First plot the raster plots of the WhiteNoise trials and then the
% Neurothresh trials
figure(plot_counter);
idxs = find(neurothreshmode==1);
raster_data = stro.ras(idxs,1);
num_spikes =[];
num_dur = [];
spikerate = [];
add_time_before = 0.1; % Want to see spikes from 100 ms before the stimulus was presented on the screen
add_time_after = 0.3; % Want to see spikes till 300 ms after the stimulus was erased from the screen
for ii = 1:size(raster_data,1)
    tmp = raster_data{ii} ;
    spikes = tmp(tmp<(stro.trial(idxs(ii),stimoffidx)+ add_time_after) & tmp>(stro.trial(idxs(ii),stimonidx)-add_time_before));
    plot(spikes-stro.trial(idxs(ii),stimonidx),(ii-1)*ones(1,length(spikes)),'k.'); hold on;
    num_spikes = [num_spikes; numel(spikes)];
    num_dur = [num_dur; stro.trial(idxs(ii),stimoffidx)- stro.trial(idxs(ii),stimonidx)];
end
xlabel('Time in s'); ylabel('Trials'); title('Are there OFF responses?');
line([0 0],[0 ii],'Color',[1 0 0]); % Time when fpacq code was dropped
line([t_offset t_offset],[0 ii],'Color',[1 0 0]); %
line([mean(num_dur) mean(num_dur)],[0 ii],'Color',[1 0 0]) % Time when stimon code was droppped
hold off;
plot_counter = plot_counter + 1;

% plotting Firing rate as a function of contrast along multiple directions
% and fitting a spline to those points
figure(plot_counter);
for ss = 1:numel(num_targetspikerates)
    tmp_n = [];
    tmp_wts = [];
    num_dur =[];
    firing_rate = [];
    for kk = 1:size(norms{ss},1)
        dir = kk;
        idxs1 = find(stro.trial(:,basisvecdiridx) == dir);
        for jj = 1:numel(idxs1)
            tmp_wts = [tmp_wts; stro.ras{idxs1(jj),weightsidx}'];
            tmp_n = [tmp_n; norm(stro.ras{idxs1(jj),weightsidx})];
        end
        raster_data = stro.ras(idxs1,1);
        for ii = 1:size(raster_data,1)
            tmp = raster_data{ii} ;
            % FR calculation
            spikes = tmp(tmp<stro.trial(idxs1(ii),stimoffidx) & tmp>stro.trial(idxs1(ii),stimonidx)+t_offset);
            num_dur = [num_dur; (stro.trial(idxs1(ii), stimoffidx)- stro.trial(idxs1(ii),stimonidx))-t_offset];
            firing_rate = [firing_rate; numel(spikes)/num_dur(end)];
        end
    end
    subplot(1,numel(num_targetspikerates),ss);
    plot3(tmp_wts(:,1),tmp_wts(:,2),firing_rate,'o','LineWidth',0.5,'MarkerFaceColor',[1 0 0],'MarkerSize',4); hold on;
    st = tpaps(tmp_wts',firing_rate');
    subplot(1,numel(num_targetspikerates),ss); fnplt(st); xlabel('Constrast S1'), ylabel('Contrast S2'), zlabel('FR'),hold off;
end
plot_counter = plot_counter + 1;

% Trying out some new analysis, I want to see if there is any meaningful relationship among the RGBs of the endpoints of the searches
% As a start, I am only concerned with just the one/first targetfiring rate
tmp_norm = norms{1};
allRGBs = cat(2,repmat(tmp_norm(:,1),[1 3]).*repmat(S1RGB',[size(tmp_norm,1) 1]),repmat(tmp_norm(:,2),[1 3]).*repmat(S2RGB',[size(tmp_norm,1) 1]));
coeff = mvregress(allRGBs,repmat(num_targetspikerates,[size(allRGBs,1) 1]));
% I need to think of a better analysis as ths rank of this vector is still
% 2.

% Trying out a new representation of looking at the data
figure(plot_counter);
for ii = 1:size(norms{1},1)
    dir = ii;
    for jj = 1: numel(num_targetspikerates)
        tmp_n = [];
        tmp_wts = [];
        tmp_parentvertices = [];
        idxs1 = find(stro.trial(:,basisvecdiridx) == dir);
        if ~isempty(idxs1)
            for jj = 1:numel(idxs1)
                tmp_parentvertices = [tmp_parentvertices; stro.ras{idxs1(jj),parentverticesidx}];
                tmp_wts = [tmp_wts; stro.ras{idxs1(jj),weightsidx}'];
                tmp_n = [tmp_n; norm(stro.ras{idxs1(jj),weightsidx})];
            end
            raster_data = stro.ras(idxs1,1);
            num_dur =[];
            firing_rate = [];
            if ~isempty(tmp_n)
                figure(plot_counter), subplot(221); plot(mean(tmp_wts(end-2:end,1)),mean(tmp_wts(end-2:end,2)),'o','MarkerFacecolor',[0 0 1]); hold on;
                subplot(222); plot(tmp_wts(end,1),tmp_wts(end,2),'o','MarkerFacecolor',[0 0 1]); hold on;
                subplot(223); plot(tmp_wts(:,1),tmp_wts(:,2),'o','MarkerFacecolor',[0 0 1]); hold on;
            end
        end
    end
end
figure(plot_counter); subplot(221); title('Mean last 3 points');hold off; 
subplot(222); title('Last point'); hold off; 
subplot(223); title('All points'); hold off;
plot_counter = plot_counter + 1;

%% Wanted to check if the staircase algorithm is working fine for all the directions
plot_counter = 1;
figure(plot_counter);
cum_flag = []; % to check if I am looking at all the probed directions
for kk = 1:size(norms{1},1)
    close(figure(plot_counter));
    t_offset1 = -0.2;
    figure(plot_counter);
    dir = kk;
    for jj = 1: numel(num_targetspikerates)
        tmp_n = [];
        tmp_wts = [];
        tmp_parentvertices = [];
        idxs1 = find(stro.trial(:,basisvecdiridx) == dir);
        disp([dir stro.trial(idxs1(end),reversalflagidx)]);
        for jj = 1:numel(idxs1)
            tmp_parentvertices = [tmp_parentvertices; stro.ras{idxs1(jj),parentverticesidx}];
            tmp_wts = [tmp_wts; stro.ras{idxs1(jj),weightsidx}'];
            tmp_n = [tmp_n; norm(stro.ras{idxs1(jj),weightsidx})];
        end
        [flag, gamutViolation] = Check_ifcompletely_probed(stro,basisvec,bkgnd_monitor,idxs1,tmp_wts(end,:));
        cum_flag = [cum_flag; flag];
        if gamutViolation == 1
            c = [0 1 0];
        else
            c = [0 0 1];
        end
        subplot(2,2,1),plot(tmp_n,'-o','LineWidth',2,'MarkerFaceColor',c); xlabel('Trials'), ylabel('Contrast');
        img = tmp_wts(end,1)*(basisvec{1}-bkgnd_monitor) + tmp_wts(end,2)*(basisvec{2}-bkgnd_monitor);
        subplot(2,2,2),imagesc(img+0.5); set(gca,'XTick',[],'YTick',[]); axis square;
        raster_data = stro.ras(idxs1,1);
        num_dur =[];
        firing_rate = [];
        subplot(2,2,4); cla;
        for ii = 1:size(raster_data,1)
            tmp = raster_data{ii} ;
            spikes1 = tmp(tmp<stro.trial(idxs1(ii),stimoffidx) & tmp>stro.trial(idxs1(ii),stimonidx)+t_offset1);
            plot(spikes1-stro.trial(idxs1(ii),stimonidx),(ii-1)*ones(1,length(spikes1)),'k.'); hold on;
            %FR calculation
            spikes = tmp(tmp<stro.trial(idxs1(ii),stimoffidx) & tmp>stro.trial(idxs1(ii),stimonidx)+t_offset);
            num_dur = [num_dur; (stro.trial(idxs1(ii), stimoffidx)- stro.trial(idxs1(ii),stimonidx))-t_offset];
            firing_rate = [firing_rate; numel(spikes)/num_dur(end)];
        end
        xlabel('Time in s'); ylabel('Trials');
        line([0 0],[0 size(raster_data,1)+1],'Color',[1 0 0])
        line([min(num_dur)+t_offset min(num_dur)+t_offset],[0 size(raster_data,1)+1],'Color',[1 0 0])
        set(gca,'Xlim',[0-0.05 max(num_dur)+t_offset+0.05])
        hold off;
        subplot(2,2,3),plot(firing_rate,'-o','LineWidth',2); xlabel('Trials'), ylabel('FR');
    end
    plot_counter = plot_counter + 1;
end

%% Trying something new with the firing rate map.
% Don't know why I wrote this code.
nbins1 =linspace(min(tmp_wts(:,1)),max(tmp_wts(:,1)),51);
nbins2 =linspace(min(tmp_wts(:,2)),max(tmp_wts(:,2)),51);
FRmap = zeros(numel(nbins1),numel(nbins2));
for ii = 1:numel(nbins1)
    for jj = 1:numel(nbins2)
        pt = [nbins1(ii) nbins2(jj)];
        dist = tmp_wts - repmat(pt,[size(tmp_wts,1) 1]);
        dist = sum(dist.*dist,2);
        L = dist<0.1;
        dist(L) = 0;
        FRmap(ii,jj) = (dist'*firing_rate)/(sum(dist));
    end
end
[A B] = meshgrid(nbins1,nbins2);
FRmap = triu(FRmap,0); % returning the upper triangular matrix
FRmap(FRmap<=0) = inf;
figure,subplot(121),surf(A,B,FRmap);
subplot(122),contour(FRmap);
%%

% This is in continuation to the previous analysis which is based on Fred's suggestion to check the existence of any suppressive effect
% Run this section of code if the subunit paradigm was used
figure(plot_counter);
add_time_before = 0.1;
add_time_after = 0.3;
first_half_dur = (0.3-t_offset)/2;
for ss = 1:numel(num_targetspikerates)
    num_dur =[];
    firing_rate = [];
    all_dirs = norms{ss};
    exc_supp_dir_idx = find(all_dirs(:,1).*all_dirs(:,2)<0); % picks up the directions in 2nd and 4th quadrant
    first_quad_dir_idx = find(all_dirs(:,1).*all_dirs(:,2)>=0);
    subplot(2,numel(num_targetspikerates),2*ss-1);
    count = 0;
    num_spikes_first_half = [];
    num_spikes_second_half = [];
    for kk = 1:numel(exc_supp_dir_idx)
        dir = exc_supp_dir_idx(kk);
        idxs1 = find(stro.trial(:,basisvecdiridx) == dir);
        raster_data = stro.ras(idxs1,1);
        for ii = 1:size(raster_data,1)
            tmp = raster_data{ii} ;
            spikes = tmp(tmp<stro.trial(idxs1(ii),stimoffidx)+add_time_after & tmp>stro.trial(idxs1(ii),stimonidx)-add_time_before);
            spikes_first_half = tmp(tmp<stro.trial(idxs1(ii),stimonidx)+t_offset+first_half_dur & tmp>stro.trial(idxs1(ii),stimonidx)+t_offset);
            num_spikes_first_half = [num_spikes_first_half; numel(spikes_first_half)];
            spikes_second_half = tmp(tmp>stro.trial(idxs1(ii),stimonidx)+t_offset+first_half_dur & tmp<stro.trial(idxs1(ii),stimoffidx));
            num_spikes_second_half = [num_spikes_second_half; numel(spikes_second_half)];
            plot(spikes-stro.trial(idxs1(ii),stimonidx),count*ones(1,length(spikes)),'k.'); hold on;
            count = count + 1;
            num_dur = [num_dur; (stro.trial(idxs1(ii), stimoffidx)- stro.trial(idxs1(ii),stimonidx))-t_offset];
        end
    end
    xlabel('Time'), ylabel('Trials');
    line([0 0],[0 numel(num_spikes_first_half)],'Color',[1 0 0]); % Time when fpacq code was dropped
    line([t_offset t_offset],[0 numel(num_spikes_first_half)],'Color',[1 0 0]); %
    line([mean(num_dur)+t_offset mean(num_dur)+t_offset],[0 numel(num_spikes_first_half)],'Color',[1 0 0]) % Time when stimon code was droppped
    title('2nd & 4th quad');
    hold off;
    subplot(2,numel(num_targetspikerates),2*ss), plot(num_spikes_first_half, num_spikes_second_half,'o','LineWidth',0.5,'MarkerFaceColor',[1 0 0],'MarkerSize',4);
    xlabel('first half'), ylabel('second half'), title('Spike Counts');
end
plot_counter = plot_counter + 1;


% Exact same analysis as before but now I am looking at the first quadrant directions
figure(plot_counter);
for ss = 1:numel(num_targetspikerates)
    num_dur =[];
    firing_rate = [];
    all_dirs = norms{ss};
    first_quad_dir_idx = find(all_dirs(:,1).*all_dirs(:,2)>=0); % Only looking at the first quadrant dirs
    subplot(2,numel(num_targetspikerates),2*ss-1);
    count = 0;
    num_spikes_first_half = [];
    num_spikes_second_half = [];
    for kk = 1:numel(first_quad_dir_idx)
        dir = first_quad_dir_idx(kk);
        idxs1 = find(stro.trial(:,basisvecdiridx) == dir);
        raster_data = stro.ras(idxs1,1);
        for ii = 1:size(raster_data,1)
            tmp = raster_data{ii} ;
            spikes = tmp(tmp<stro.trial(idxs1(ii),stimoffidx)+add_time_after & tmp>stro.trial(idxs1(ii),stimonidx)-add_time_before);
            spikes_first_half = tmp(tmp<stro.trial(idxs1(ii),stimonidx)+t_offset+first_half_dur & tmp>stro.trial(idxs1(ii),stimonidx)+t_offset);
            num_spikes_first_half = [num_spikes_first_half; numel(spikes_first_half)];
            spikes_second_half = tmp(tmp>stro.trial(idxs1(ii),stimonidx)+t_offset+first_half_dur & tmp<stro.trial(idxs1(ii),stimoffidx));
            num_spikes_second_half = [num_spikes_second_half; numel(spikes_second_half)];
            plot(spikes-stro.trial(idxs1(ii),stimonidx),count*ones(1,length(spikes)),'k.'); hold on;
            count = count + 1;
            num_dur = [num_dur; (stro.trial(idxs1(ii), stimoffidx)- stro.trial(idxs1(ii),stimonidx))-t_offset];
        end
    end
    xlabel('Time'), ylabel('Trials');
    line([0 0],[0 numel(num_spikes_first_half)],'Color',[1 0 0]); % Time when fpacq code was dropped
    line([t_offset t_offset],[0 numel(num_spikes_first_half)],'Color',[1 0 0]); %
    line([mean(num_dur)+t_offset mean(num_dur)+t_offset],[0 numel(num_spikes_first_half)],'Color',[1 0 0]) % Time when stimon code was droppped
    title('1st quad & 3rd quad');
    hold off;
    subplot(2,numel(num_targetspikerates),2*ss), plot(num_spikes_first_half, num_spikes_second_half,'o','LineWidth',0.5,'MarkerFaceColor',[1 0 0],'MarkerSize',4);
    xlabel('first half'), ylabel('second half'), title('Spike Counts');
end
plot_counter = plot_counter + 1;

%% See if the firing rate surface changes or not when the time window for counting the spikes is varied
% Still working on it
figure(plot_counter);
start_time = [0.0; 0.01];
end_time = [0.01; 0.02];
for mm = 1:numel(start_time)
    for ss = 1:numel(num_targetspikerates)
        tmp_n = [];
        tmp_wts = [];
        num_dur =[];
        firing_rate = [];
        for kk = 1:size(norms{ss},1)
            dir = kk;
            idxs1 = find(stro.trial(:,basisvecdiridx) == dir);
            for jj = 1:numel(idxs1)
                tmp_wts = [tmp_wts; stro.ras{idxs1(jj),weightsidx}'];
                tmp_n = [tmp_n; norm(stro.ras{idxs1(jj),weightsidx})];
            end
            raster_data = stro.ras(idxs1,1);
            for ii = 1:size(raster_data,1)
                tmp = raster_data{ii} ;
                % FR calculation
                spikes = tmp(tmp<stro.trial(idxs1(ii),stimonidx)+t_offset+end_time(mm) & tmp>stro.trial(idxs1(ii),stimonidx)+t_offset+start_time(mm));
                num_dur = [num_dur; end_time(mm) - start_time(mm)];
                firing_rate = [firing_rate; numel(spikes)/num_dur(end)];
            end
        end
        subplot(numel(start_time),numel(num_targetspikerates),mm);
        plot3(tmp_wts(:,1),tmp_wts(:,2),firing_rate,'o','LineWidth',0.5,'MarkerFaceColor',[1 0 0],'MarkerSize',4); hold on;
        st = tpaps(tmp_wts',firing_rate');
        subplot(numel(start_time),numel(num_targetspikerates),mm); fnplt(st); xlabel('Constrast S1'), ylabel('Contrast S2'), zlabel('FR'),hold off;
    end
end
plot_counter = plot_counter + 1;

%% Plotting the contrast response function in multiple directions but displaying them in a 2-D plot
% Not a useful analysis, would have been meaningful if there were multiple repeats for each stimulus
figure(plot_counter);
for ss = 1:numel(num_targetspikerates)
    tmp_n = [];
    tmp_wts = [];
    num_dur =[];
    firing_rate = [];
    for kk = 1:size(norms{ss},1)
        dir = kk;
        idxs1 = find(stro.trial(:,basisvecdiridx) == dir);
        for jj = 1:numel(idxs1)
            tmp_wts = [tmp_wts; stro.ras{idxs1(jj),weightsidx}'];
            tmp_n = [tmp_n; norm(stro.ras{idxs1(jj),weightsidx})];
        end
        raster_data = stro.ras(idxs1,1);
        for ii = 1:size(raster_data,1)
            tmp = raster_data{ii} ;
            % FR calculation
            spikes = tmp(tmp<stro.trial(idxs1(ii),stimoffidx) & tmp>stro.trial(idxs1(ii),stimonidx)+t_offset);
            num_dur = [num_dur; (stro.trial(idxs1(ii), stimoffidx)- stro.trial(idxs1(ii),stimonidx))-t_offset];
            firing_rate = [firing_rate; numel(spikes)/num_dur(end)];
        end
        [~,ind] = sort(tmp_n);
        subplot(1,numel(num_targetspikerates),ss);plot(tmp_n(ind),firing_rate(ind),'o'); hold on;
    end
    hold off;
end

plot_counter = plot_counter + 1;


%% Code that validates if the Online selection of new weights is happening correctly, working correctly
for jj = 1:numel(num_targetspikerates)
    idxs = find(neurothreshmode == 1 & stro.trial(:,targetspikerateidx)==num_targetspikerates(jj));
    figure(plot_counter); hold on; grid on;
    for ii = 1:numel(idxs)
        tmp = stro.ras{idxs(ii),weightsidx};
        if (stro.trial(idxs(ii),basisvecdiridx)<=5*numel(num_targetspikerates))
            c = 'k';
        elseif (stro.trial(idxs(ii),basisvecdiridx)<=7*numel(num_targetspikerates))
            c = 'r';
        else
            c = 'b';
        end
        subplot(numel(num_targetspikerates),2,2*jj-1);scatter(tmp(1),tmp(2),c); hold on; pause(0.05); grid on;
        subplot(numel(num_targetspikerates),2,2*jj);compass(tmp(1),tmp(2),c); hold on;
    end
end
hold off;
plot_counter = plot_counter + 1;



