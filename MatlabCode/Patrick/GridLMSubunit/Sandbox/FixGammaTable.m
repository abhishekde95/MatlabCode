 
% This script will correct the gamma table, and thereby the true 
% distribution, of stimuli for GridLMSubunit.  The idea is to convert
% everything back to voltages using the "expected" gammatable and monspd
% (those we used in the experiment), and convert them back using the 
% "actual" gammatable and monspd.  This will result in the real stimuli
% that were used.


% To begin, unpack the datafile
library = '/Users/jpatrickweller/Dropbox/Patrick/GLMS Data/nex files/';
datafile = 'M090415001.nex'; % -L-M
datafile = 'M012516001.nex'; % looks kinda rando
rawdata = nex2stro([char(library) datafile]);

% Next, unpack the real cals data
load('/Users/jpatrickweller/Documents/Matlab Working Directory/Slave/Monitor Calibration/Monitor data/Dell 7/Dell7.mat')
calsact = cals{end}; % This will change with date. Use last one for files after Jan 22th, 2016

% Unpack the expected gamma table and M matrix
fundamentals = rawdata.sum.exptParams.fundamentals;
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
monspd = rawdata.sum.exptParams.monspd;
monspd = reshape(monspd,[length(monspd)/3, 3]);
monspd_exp = SplineRaw([380:4:780]', monspd, [380:5:780]');
M_exp = fundamentals' * monspd_exp;
gammaTable = rawdata.sum.exptParams.gammatable;
gammaTable_exp = reshape(gammaTable, length(gammaTable)/3, 3);

% Unpack the actual gamma table and calculate new M matrix
monspd_act = SplineRaw([380:4:780]', calsact.P_device, [380:5:780]');
gammaTable_act = calsact.gammaTable;
M_act = fundamentals' * monspd_act;

% Convert background gun intensities to gun voltages
bkgndrgb_exp = rawdata.sum.exptParams.bkgndrgb;
invgamma_exp = InvertGamma(gammaTable_exp, 0);
bkgndrgbvols_exp = rgb2RGB(bkgndrgb_exp,invgamma_exp); % voltage is normalized

% Pull expected Lcc, Mcc, and Scc out of datafile
GLMP_L = rawdata.trial(:,strcmp(rawdata.sum.trialFields(1,:),'Epoch'))==2 |...
    rawdata.trial(:,strcmp(rawdata.sum.trialFields(1,:),'Epoch'))==3;
Lcc_exp = rawdata.trial(GLMP_L,strcmp(rawdata.sum.trialFields(1,:),'Lcc'));
Mcc_exp = rawdata.trial(GLMP_L,strcmp(rawdata.sum.trialFields(1,:),'Mcc'));
Scc_exp = rawdata.trial(GLMP_L,strcmp(rawdata.sum.trialFields(1,:),'Scc'));

% Convert stimulus cc to cone excitations
bkgndlms_exp = M_exp * bkgndrgb_exp;
Lexcite_exp = Lcc_exp * bkgndlms_exp(1) + bkgndlms_exp(1);
Mexcite_exp = Mcc_exp * bkgndlms_exp(2) + bkgndlms_exp(2);
Sexcite_exp = Scc_exp * bkgndlms_exp(3) + bkgndlms_exp(3);
LMSexcite_exp = [Lexcite_exp Mexcite_exp Sexcite_exp];

% Convert stimulus cone excitations to gun intensities
rgbintens_exp = [inv(M_exp) * LMSexcite_exp']';

% Convert gun intensities to gun voltages
stimvolts_exp = rgb2RGB(rgbintens_exp,invgamma_exp);

% Now that we have background and stimulus modulations in terms of normlized volts,
% we can use the "actual" gamma table and monspd to calculate the LMScc's
% that were really used.

% convert background using actual gamma table
bkgndidx = round(bkgndrgbvols_exp * 255)+1;
bkgndRintens_act = gammaTable_act(bkgndidx(1),1);
bkgndGintens_act = gammaTable_act(bkgndidx(2),2);
bkgndBintens_act = gammaTable_act(bkgndidx(3),3);
bkgndrgb_act = [bkgndRintens_act; bkgndGintens_act; bkgndBintens_act];

% Convert stimuli using actual gamma table
stimidx = round(stimvolts_exp * 255) + 1;
stimRintens_act = gammaTable_act(stimidx(:,1),1);
stimGintens_act = gammaTable_act(stimidx(:,2),2);
stimBintens_act = gammaTable_act(stimidx(:,3),3);
stimrgb_act = [stimRintens_act stimGintens_act stimBintens_act];

% Convert background and stimuli from rgb to lms
bkgndlms_act = M_act * bkgndrgb_act;
stimlms_act = [M_act * stimrgb_act']';

% Calculate new Lcc, Mcc, Scc based on new lms values
Lcc_act = (stimlms_act(:,1) - bkgndlms_act(1)) ./ bkgndlms_act(1);
Mcc_act = (stimlms_act(:,2) - bkgndlms_act(2)) ./ bkgndlms_act(2);
Scc_act = (stimlms_act(:,3) - bkgndlms_act(3)) ./ bkgndlms_act(3);

% Plot the new results
figure(1); clf; hold on; grid on;
plot(Lcc_exp,Mcc_exp,'ko')
plot(Lcc_act,Mcc_act,'r*')
xlabel('Lcc'); ylabel('Mcc')














