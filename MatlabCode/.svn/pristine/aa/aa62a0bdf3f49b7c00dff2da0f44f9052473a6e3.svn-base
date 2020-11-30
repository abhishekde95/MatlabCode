% Author - Abhishek De, 5/18
% storing files for DToneloc in 3 categories: opto, behavior and behavior with laser control
close all; clearvars;
%************************************* MAUI files *******************************************
% opto files
filenameopto1 = {'M041218004.nex';'M041218006.nex';'M041218007.nex';'M041218008.nex';'M041218009.nex'};
RWmultiplieropto1 = [1.0;1.0;1.0;1.0;1.0;];
laserdialopto1 = [2.9;2.0;2.0;2.0;2.0];
filenameopto2 = {'M042518002.nex';'M042518004.nex';'M042518005.nex';'M042518008.nex';'M042518006.nex';'M042518007.nex'}; % 4/25
RWmultiplieropto2 =[1.0;1.0;1.0;1.0;1.0;1.0];
laserdialopto2 = [2.5;2.5;2.5;2.5;2.5;2.5];
filenameopto3 = {'M042618003.nex';'M042618004.nex';'M042618005.nex';'M042618006.nex';'M042618007.nex'}; % 4/26, 200 stim, 300 laser
RWmultiplieropto3 =[1.0;1.0;1.0;1.0;1.0];
laserdialopto3 = [2.5;2.5;2.5;2.5;2.5];
filenameopto4 = {'M042718006.nex';'M042718007.nex';'M042718008.nex';'M042718009.nex';'M042718010.nex'}; % 4/27, 200 stim, 300 laser
RWmultiplieropto4 =[1.0;1.0;1.0;1.0;1.0];
laserdialopto4 = [2.7;2.7;2.7;2.7;3.2];
filenameopto5 = {'M042718012.nex';'M042718013.nex';'M042718014.nex'}; % 4/27, 100 stim, 150 laser
RWmultiplieropto5 =[1.0;1.0;1.0];
laserdialopto5 = [3.2;3.2;3.2];
filenameopto6 = {'M051318008.nex';'M051318009.nex';'M051318010.nex';'M051318011.nex'}; %5/13, 200 stim, 300 laser 
RWmultiplieropto6 =[1.0;0.70;0.60;1.0];
laserdialopto6 = 4.0*ones(size(RWmultiplieropto6));
filenameopto7 = {'M051418006.nex';'M051418008.nex';'M051418009.nex';'M051418010.nex';'M051418011.nex';'M051418013.nex';'M051418014.nex';'M051418015.nex';'M051418016.nex'}; %5/14, 200 stim, 300 laser 
RWmultiplieropto7 =[1.0;1.0;1.0;0.7;0.7;1.0;1.0;0.7;0.7];
laserdialopto7 = 4.0*ones(size(RWmultiplieropto7));
filenameopto8 = {'M051518003.nex';'M051518004.nex';'M051518005.nex';'M051518006.nex';'M051518007.nex';'M051518008.nex';'M051518009.nex';'M051518010.nex';'M051518011.nex';}; %5/15, 200 stim, 300 laser 
RWmultiplieropto8 =[1.0;0.7;0.5;0.5;0.3;0.6;1.0;0.4;0.7];
laserdialopto8 = 4.0*ones(size(RWmultiplieropto8));
filenameopto9 = {'M051618003.nex';'M051618004.nex';'M051618005.nex';'M051618006.nex';'M051618007.nex';'M051618008.nex';'M051618009.nex';'M051618010.nex';'M051618015.nex';'M051618016.nex';'M051618011.nex';'M051618012.nex'}; %5/16, 200 stim, 300 laser
RWmultiplieropto9 =[1.0;0.7;0.5;0.4;0.3;0.45;0.55;0.60;1.0;1.0;0.6;1.0];
laserdialopto9 = 4.0*ones(size(RWmultiplieropto9));
filenameopto10 = {'M051718018.nex';'M051718019.nex'}; %5/17, 200 stim, 300 laser
RWmultiplieropto10 = [1.0;1.0];
laserdialopto10 = [1.0;1.0];
filenameopto11 = {'M051818002.nex';'M051818003.nex';'M051818015.nex';'M051818006.nex';'M051818007.nex';'M051818008.nex';'M051818009.nex';'M051818012.nex'};
RWmultiplieropto11 = [1.0;1.0;1.0;1.0;1.0;1.0;1.0;1.0];
laserdialopto11 = [2.0;2.0;2.0;2.0;2.0;2.0;2.0;2.0];
filenameopto12 = {'M052518006.nex'};
RWmultiplieropto12 = [1.0];
laserdialopto12 = [3.0];
filenameopto = [filenameopto1;filenameopto2;filenameopto3;filenameopto4;filenameopto5;filenameopto6;filenameopto7;filenameopto8;filenameopto9;filenameopto10;filenameopto11;filenameopto12];
RWmultiplieropto = [RWmultiplieropto1;RWmultiplieropto2;RWmultiplieropto3;RWmultiplieropto4;RWmultiplieropto5;RWmultiplieropto6;RWmultiplieropto7;RWmultiplieropto8;RWmultiplieropto9;RWmultiplieropto10;RWmultiplieropto11;RWmultiplieropto12];
laserdialopto = [laserdialopto1;laserdialopto2;laserdialopto3;laserdialopto4;laserdialopto5;laserdialopto6;laserdialopto7;laserdialopto8;laserdialopto9;laserdialopto10;laserdialopto11;laserdialopto12];
save filenameoptoM filenameopto
save RWmultiplieroptoM RWmultiplieropto
save laserdialoptoM laserdialopto

% behavior files - Maui
filenamebeh1 = {'M050418002.nex';'M050418003.nex';'M050418004.nex';'M050418005.nex';'M050418006.nex';'M050418007.nex'}; % 5/4
RWmultiplierbeh1 = [0.55;0.4;0.4;0.55;0.65;0.65];
filenamebeh2 = {'M050518001.nex';'M050518002.nex';'M050518003.nex';'M050518004.nex';'M050518005.nex';'M050518006.nex';'M050518007.nex';'M050518008.nex';'M050518009.nex';}; % 5/5
RWmultiplierbeh2 = [0.85;0.7;0.65;0.65;0.5;0.5;0.4;0.4;0.55];
filenamebeh3 = {'M050618001.nex';'M050618002.nex';'M050618003.nex';'M050618004.nex';'M050618005.nex';'M050618006.nex'}; % 5/6
RWmultiplierbeh3 = [0.55;0.4;0.4;0.3;0.3;0.3];
filenamebeh4 = {'M050718001.nex';'M050718002.nex';'M050718003.nex';'M050718004.nex'}; % 5/7
RWmultiplierbeh4 = [0.85;0.5;0.2;0.35];
filenamebeh5 = {'M050818001.nex';'M050818002.nex';'M050818003.nex';'M050818004.nex';'M050818006.nex';'M050818007.nex';'M050818008.nex';'M050818009.nex'}; % 5/8
RWmultiplierbeh5 = [1.0;0.8;0.6;0.4;0.5;0.2;0.6;1.0];
filenamebeh6 = {'M050918001.nex';'M050918002.nex';'M050918003.nex';'M050918004.nex';'M050918005.nex';'M050918006.nex';'M050918007.nex';'M050918008.nex';'M050918009.nex';'M050918010.nex'}; % 5/9
RWmultiplierbeh6 = [1.0;0.75;0.50;0.25;0.50;0.30;0.40;0.40;1.0;1.0];
filenamebeh7 = {'M062818001.nex';'M062818002.nex';'M062818003.nex';'M062818004.nex';'M062818005.nex'};
RWmultiplierbeh7 = [0.5;0.3;0.2;0.2;0.1];
filenamebeh8 = {'M062918001.nex';'M062918002.nex';'M062918003.nex';'M062918004.nex';'M062918005.nex';'M062918006.nex';'M062918007.nex';'M062918008.nex'};
RWmultiplierbeh8 = [1.0;1.0;1.0;1.0;1.0;1.0;1.0;1.0];
filenamebeh9 = {'M070218001.nex';'M070218002.nex';'M070218003.nex';'M070218004.nex';'M070218005.nex';'M070218006.nex';'M070218007.nex';'M070218008.nex';'M070218009.nex'};
RWmultiplierbeh9 = [1.0;1.0;1.0;1.0;1.0;1.0;1.0;1.0;1.0];
filenamebeh10 = {'M070318001.nex';'M070318002.nex';'M070318003.nex';'M070318004.nex';'M070318005.nex';'M070318006.nex';'M070318007.nex';'M070318008.nex';'M070318009.nex';'M070318010.nex'};
RWmultiplierbeh10 = [1.0;1.0;0.5;0.5;0.3;0.3;0.2;0.2;0.2;0.15;0.15];
filenamebeh11 = {'M070518001.nex';'M070518002.nex';'M070518003.nex';'M070518004.nex';'M070518005.nex';'M070518006.nex';'M070518007.nex';'M070518008.nex';'M070518009.nex';'M070518010.nex';'M070518011.nex';'M070518012.nex'};
RWmultiplierbeh11 = [1.0;1.0;1.0;0.5;0.3;0.3;0.3;0.3;0.2;0.2;0.2;0.2];
filenamebeh12 = {'M070618001.nex';'M070618003.nex';'M070618005.nex';'M070618007.nex';'M070618008.nex';'M070618009.nex';'M070618010.nex'};
RWmultiplierbeh12 = [1.0;1.0;1.0;1.0;1.0;1.0;1.0];
filenamebeh13 = {'M070918001.nex';'M070918002.nex';'M070918003.nex';'M070918004.nex';'M070918005.nex';'M070918006.nex';'M070918007.nex';'M070918008.nex';'M070918009.nex';'M070918010.nex';'M070918011.nex';'M070918012.nex';'M070918013.nex'};
RWmultiplierbeh13 = [1.0;1.0;0.5;0.5;0.4;1.0;1.0;0.4;0.3;0.3;0.3;0.3;0.3];
filenamebeh14 = {'M071018001.nex';'M071018002.nex';'M071018003.nex';'M071018004.nex';'M071018005.nex';'M071018006.nex'};
RWmultiplierbeh14 = [1.0;1.0;1.0;1.0;1.0;1.0];
filenamebeh15 = {'M071118001.nex';'M071118002.nex';'M071118003.nex';'M071118004.nex';'M071118005.nex';'M071118006.nex';'M071118007.nex';'M071118008.nex'};
RWmultiplierbeh15 = [1.0;1.0;0.5;0.3;0.15;0.10;1.0;1.0];
filenamebeh16 = {'M071218001.nex';'M071218002.nex';'M071218003.nex';'M071218004.nex';'M071218005.nex';'M071218006.nex'};
RWmultiplierbeh16 = [1.0;1.0;1.0;1.0;1.0;1.0];
filenamebeh17 = {'M071318001.nex';'M071318002.nex';'M071318003.nex';'M071318004.nex';'M071318005.nex';'M071318006.nex';'M071318007.nex';'M071318008.nex';'M071318009.nex';'M071318010.nex';'M071318011.nex'};
RWmultiplierbeh17 = [1.0;1.0;1.0;1.0;1.0;1.0;1.0;1.0;1.0;1.0;1.0];
filenamebeh18 = {'M071618001.nex';'M071618002.nex';'M071618003.nex';'M071618004.nex';'M071618005.nex';'M071618006.nex';'M071618007.nex';'M071618008.nex'};
RWmultiplierbeh18 = [1.0;1.0;1.0;1.0;1.0;1.0;1.0;1.0];
filenamebeh19 = {'M071718001.nex';'M071718002.nex';'M071718003.nex';'M071718004.nex';'M071718005.nex'};
RWmultiplierbeh19 = [1.0;1.0;1.0;1.0;1.0];
filenamebeh20 = {'M071818001.nex';'M071818002.nex';'M071818003.nex';'M071818004.nex';'M071818005.nex';'M071818006.nex';'M071818007.nex'};
RWmultiplierbeh20 = [1.0;1.0;1.0;1.0;1.0;1.0;1.0];
filenamebeh21 = {'M071918001.nex';'M071918002.nex';'M071918003.nex';'M071918004.nex'};
RWmultiplierbeh21 = [1.0;1.0;1.0;1.0];
% filenamebeh = [filenamebeh1;filenamebeh2;filenamebeh3;filenamebeh4;filenamebeh5;filenamebeh6;filenamebeh7;filenamebeh8];
% RWmultiplierbeh = [RWmultiplierbeh1;RWmultiplierbeh2;RWmultiplierbeh3;RWmultiplierbeh4;RWmultiplierbeh5;RWmultiplierbeh6;RWmultiplierbeh7;RWmultiplierbeh8];
filenamebeh = [filenamebeh9;filenamebeh12;filenamebeh14;filenamebeh16;filenamebeh17;filenamebeh18;filenamebeh19;filenamebeh20];
RWmultiplierbeh = [RWmultiplierbeh9;RWmultiplierbeh12;RWmultiplierbeh14;RWmultiplierbeh16;RWmultiplierbeh17;RWmultiplierbeh18;RWmultiplierbeh19;RWmultiplierbeh20];
laserdialbeh = zeros(size(RWmultiplierbeh));
save filenamebehM filenamebeh
save RWmultiplierbehM RWmultiplierbeh
save laserdialbehM laserdialbeh

% laser control behavior files obtained during behavior
filenamelc1 = {'M051018006.nex';'M051018007.nex';'M051018008.nex';'M051018009.nex';'M051018010.nex';'M051018011.nex';'M051018012.nex'}; % 5/10, laser control trials
RWmultiplierlc1 =[1.0;0.75;0.50;0.25;0.40;0.60;0.80];
laserdiallc1 = [2.0;2.0;2.0;2.0;2.0;2.0;2.0];
filenamelc2 = {'M051118005.nex';'M051118006.nex';'M051118007.nex';'M051118008.nex';'M051118009.nex';'M051118010.nex';'M051118011.nex';'M051118012.nex';'M051118013.nex';}; % 5/11, laser control trials
RWmultiplierlc2 =[1.0;0.70;0.40;0.60;0.50;0.30;0.65;0.80;0.45];
laserdiallc2 = [2.0;2.0;2.0;2.0;2.0;2.0;2.0;2.0;2.0];
filenamelc = [filenamelc1;filenamelc2];
RWmultiplierlc = [RWmultiplierlc1;RWmultiplierlc2];
laserdiallc = [laserdiallc1;laserdiallc2];
save filenamelcM filenamelc
save RWmultiplierlcM RWmultiplierlc
save laserdiallcM laserdiallc

% position control opto files 
filenamepcopto1 = {'M051618011.nex';'M051618012.nex'}; % 5/16, RF=47,-67, 1 degree below the SCstimcue derived scotoma
RWmultiplierpcopto1 = [0.6;1.0];
laserdialpcopto1 = [4.0;4.0];
filenamepcopto2 = {'M051618013.nex';'M051618014.nex'}; % 5/16, RF=10,-47, far from the SCstimcue derived scotoma
RWmultiplierpcopto2 = [1.0;1.0];
laserdialpcopto2 = [4.0;4.0];
filenamepcopto3 = {'M051718021.nex'}; % 5/17, RF=47,-57, far from the SCstimcue derived scotoma
RWmultiplierpcopto3 = [1.0];
laserdialpcopto3 = [2.0];
filenamepcopto4 = {'M051718022.nex'}; % 5/17, RF=47,0 far from the SCstimcue derived scotoma
RWmultiplierpcopto4 = [1.0];
laserdialpcopto4 = [2.0];
filenamepcopto5 = {'M051818004.nex';'M051818005.nex'}; % 5/18, RF=47,0 far from the SCstimcue derived scotoma
RWmultiplierpcopto5 = [1.0;1.0];
laserdialpcopto5 = [2.0; 2.0];
filenamepcopto6 = {'M051818006.nex';'M051818007.nex'}; % 5/18, RF=47,-57 far from the SCstimcue derived scotoma
RWmultiplierpcopto6 = [1.0; 1.0];
laserdialpcopto6 = [2.0; 2.0];
filenamepcopto7 = {'M051818008.nex';'M051818009.nex'}; % 5/18, RF=27,-60 far from the SCstimcue derived scotoma
RWmultiplierpcopto7 = [1.0; 1.0];
laserdialpcopto7 = [2.0; 2.0];
filenamepcopto8 = {'M051818010.nex';'M051818011.nex'}; % 5/18, RF=20,-30 far from the SCstimcue derived scotoma
RWmultiplierpcopto8 = [1.0; 1.0];
laserdialpcopto8 = [2.0; 2.0];
filenamepcopto9 = {'M051818012.nex'}; % 5/18, RF=27,-57 far from the SCstimcue derived scotoma
RWmultiplierpcopto9 = [1.0];
laserdialpcopto9 = [2.0];
filenamepcopto10 = {'M051818013.nex';'M051818014.nex'}; % 5/18, RF=27,-80 far from the SCstimcue derived scotoma
RWmultiplierpcopto10 = [2.0; 1.0];
laserdialpcopto10 = [2.0; 1.0];
filenamepcopto = [filenamepcopto1;filenamepcopto2;filenamepcopto3;filenamepcopto4;filenamepcopto5;filenamepcopto6;filenamepcopto7;filenamepcopto8;filenamepcopto9;filenamepcopto10];
RWmultiplierpcopto = [RWmultiplierpcopto1;RWmultiplierpcopto2;RWmultiplierpcopto3;RWmultiplierpcopto4;RWmultiplierpcopto5;RWmultiplierpcopto6;RWmultiplierpcopto7;RWmultiplierpcopto8;RWmultiplierpcopto9;RWmultiplierpcopto10];
laserdialpcopto = [laserdialpcopto1;laserdialpcopto2;laserdialpcopto3;laserdialpcopto4;laserdialpcopto5;laserdialpcopto6;laserdialpcopto7;laserdialpcopto8;laserdialpcopto9;laserdialpcopto10];
save filenamepcoptoM filenamepcopto 
save RWmultiplierpcoptoM RWmultiplierpcopto
save laserdialpcoptoM laserdialpcopto


% randomly interleaved multiple locations opto experiments
filenamemultiloc1 = {'M052318001.nex';'M052318002.nex';'M052318003.nex';'M052318004.nex'}; % 5/23, opto session, tested 3 spatial locations
RWmultipliermultiloc1 = [1.0;1.0;1.0;1.0];
laserdialmultiloc1 = [2.0;2.0;2.0;2.0];
filenamemultiloc2 = {'M052518007.nex';'M052518009.nex'}; % 5/25, opto session, tested 3 spatial locations
RWmultipliermultiloc2 = [1.0;1.0];
laserdialmultiloc2 = [2.0;2.0];
filenamemultiloc = [filenamemultiloc1; filenamemultiloc2];
RWmultipliermultiloc = [RWmultipliermultiloc1; RWmultipliermultiloc2];
laserdialmultiloc = [laserdialmultiloc1; laserdialmultiloc2];
save filenamemultilocM filenamemultiloc 
save RWmultipliermultilocM RWmultipliermultiloc
save laserdialmultilocM laserdialmultiloc

%************************************* APOLLO files *******************************************
filenameopto1 = {'A011419007.nex';'A011419008.nex';'A011419009.nex';'A011419010.nex'};
RWmultiplieropto1 = [0.7;0.7;0.5;0.5];
laserdialopto1 = [2.5;2.5;2.5;2.8];
filenameopto2 = {'A011519005.nex';'A011519006.nex';'A011519007.nex';'A011519008.nex'};
RWmultiplieropto2 = [1.0;0.7;0.5;0.5];
laserdialopto2 = [2.7;2.2;2.2;2.2];
filenameopto3 = {'A011719009.nex';'A011719010.nex';'A011719011.nex';'A011719012.nex'};
RWmultiplieropto3 = [1.0;1.0;1.0;1.0];
laserdialopto3 = [2.4;2.4;2.4;2.4];
filenameopto4 = {'A012219006.nex';'A012219009.nex';'A012219010.nex';'A012219011.nex';'A012219012.nex';'A012219013.nex'};
RWmultiplieropto4 = [1.0;0.8;0.8;0.8;0.8;0.8;];
laserdialopto4 = [2.5;2.5;2.7;2.7;2.7;2.7];
filenameopto5 = {'A012319004.nex';'A012319007.nex';'A012319008.nex';'A012319009.nex';'A012319010.nex';'A012319011.nex'};
RWmultiplieropto5 = [1.0;0.7;0.7;0.7;0.7;0.7];
laserdialopto5 = [2.0;1.5;1.5;1.5;1.5;1.5];
filenameopto6 = {'A013019003.nex';'A013019004.nex';'A013019005.nex';'A013019006.nex';'A013019007.nex';'A013019008.nex';'A013019009.nex';'A013019012.nex'};
RWmultiplieropto6 = [1.0;0.7;0.7;0.7;0.7;0.5;0.5;0.7];
laserdialopto6 = [2.0;2.0;2.0;2.2;2.2;2.2;2.2;2.2];
filenameopto7 = {'A013119001.nex';'A013119002.nex';'A013119003.nex';'A013119004.nex';'A013119005.nex';'A013119006.nex';'A013119007.nex';'A013119008.nex'};
RWmultiplieropto7 = [0.7;0.7;0.7;0.7;0.6;0.6;0.7;1.0];
laserdialopto7 = [1.0;1.0;1.5;1.5;1.0;1.0;1.0;1.0];
filenameopto8 = {'A020119001.nex';'A020119002.nex';'A020119003.nex';'A020119004.nex';'A020119005.nex';'A020119006.nex';'A020119009.nex';'A020119010.nex';'A020119011.nex';'A020119012.nex';'A020119013.nex'};
RWmultiplieropto8 = [0.8;0.8;0.8;0.8;1.0;0.7;0.8;0.8;0.7;0.8;0.7];
laserdialopto8 = [2.0;2.0;2.0;2.0;2.0;2.0;2.0;2.0;2.0;2.0;0.2];
filenameopto9 = {'A020319001.nex';'A020319002.nex';'A020319003.nex';'A020319006.nex';'A020319007.nex';'A020319008.nex';'A020319009.nex'};
RWmultiplieropto9 = [0.8;0.8;0.8;0.8;0.8;0.8;0.8;0.8];
laserdialopto9 = [1.3;1.8;2.3;2.0;2.0;2.5;2.0];
filenameopto10 = {'A020519001.nex';'A020519002.nex';'A020519003.nex';'A020519004.nex';'A020519005.nex';'A020519006.nex';'A020519007.nex';'A020519008.nex';'A020519009.nex'};
RWmultiplieropto10 = [0.8;0.7;0.7;0.8;1.0;0.8;0.8;0.8;0.8];
laserdialopto10 = [1.5;1.5;1.5;1.5;1.5;1.3;1.3;1.3;1.5];
filenameopto11 = {'A021419003.nex';'A021419004.nex';'A021419005.nex';'A021419007.nex';'A021419008.nex';'A021419009.nex'};
RWmultiplieropto11 = [1.0;0.8;0.7;0.8;0.7;0.6];
laserdialopto11 = [2.4;2.4;2.8;2.4;2.8;1.8];
filenameopto12 = {'A021519003.nex';'A021519004.nex';'A021519005.nex';'A021519006.nex';'A021519007.nex';'A021519008.nex';'A021519009.nex';'A021519010.nex'};
RWmultiplieropto12 = 0.8*ones(numel(filenameopto12),1);
laserdialopto12 = [2.0;1.0;0.5;1.5;0.75;2.0;0.25;2.5];
filenameopto = [filenameopto1; filenameopto2; filenameopto3; filenameopto4; filenameopto5; filenameopto6; filenameopto7; filenameopto8;filenameopto9;filenameopto10;filenameopto11;filenameopto12];
RWmultiplieropto = [RWmultiplieropto1; RWmultiplieropto2; RWmultiplieropto3; RWmultiplieropto4; RWmultiplieropto5; RWmultiplieropto6; RWmultiplieropto7; RWmultiplieropto8; RWmultiplieropto9; RWmultiplieropto10; RWmultiplieropto11;RWmultiplieropto12];
laserdialopto = [laserdialopto1; laserdialopto2; laserdialopto3; laserdialopto4; laserdialopto5; laserdialopto6; laserdialopto7; laserdialopto8;laserdialopto9;laserdialopto10;laserdialopto11;laserdialopto12];
save filenameoptoA filenameopto
save RWmultiplieroptoA RWmultiplieropto
save laserdialoptoA laserdialopto

% behavior files - Apollo
filenamebeh1 = {'A020619001.nex';'A020619002.nex';'A020619003.nex';'A020619004.nex';'A020619005.nex';'A020619006.nex';'A020619007.nex';'A020619008.nex';'A020619009.nex';'A020619010.nex';'A020619011.nex';'A020619012.nex';'A020619013.nex';'A020619014.nex';'A020619015.nex';'A020619016.nex';'A020619017.nex';'A020619018.nex';'A020619019.nex';'A020619020.nex'};
RWmultiplierbeh1 = ones(size(filenamebeh1));
laserdialbeh1 = zeros(size(filenamebeh1));
filenamebeh2 = {'A020719001.nex';'A020719002.nex';'A020719006.nex';'A020719003.nex';'A020719004.nex';'A020719005.nex';'A020719007.nex';'A020719008.nex';'A020719009.nex';'A020719010.nex';'A020719011.nex';'A020719012.nex';'A020719013.nex';'A020719014.nex';'A020719015.nex';'A020719016.nex';'A020719017.nex';'A020719018.nex'};
RWmultiplierbeh2 = [1.0;0.7;0.5;1.0;0.7;0.5;1.0;0.7;0.5;1.0;0.7;0.5;1.0;0.7;0.5;1.0;0.7;0.5]; 
laserdialbeh2 = zeros(size(filenamebeh1));
filenamebeh = [filenamebeh1];
RWmultiplierbeh = [RWmultiplierbeh1];
laserdialbeh = [laserdialbeh2];
save filenamebehA filenamebeh
save RWmultiplierbehA RWmultiplierbeh
save laserdialbehA laserdialbeh

%************************************* PANGU files *******************************************
filenameopto1 = {'P102418010.nex';'P102418011.nex';'P102418012.nex';'P102418013.nex';'P102418014.nex';'P102418015.nex';'P102418016.nex';'P102418017.nex';'P102418018.nex'};%10/24
RWmultiplieropto1 =[1.0;1.0;1.0;1.0;1.0;1.0;1.0;1.0;1.0];
laserdialopto1 = [2.0;2.0;2.0;2.0;2.0;2.0;2.0;2.0;2.0];
filenameopto = [filenameopto1];
RWmultiplieropto = [RWmultiplieropto1];
laserdialopto = [laserdialopto1];
save filenameoptoP filenameopto
save RWmultiplieroptoP RWmultiplieropto
save laserdialoptoP laserdialopto

% behavior files - Pangu
filename1 = {'P060418001.nex';'P060418002.nex';'P060418003.nex';'P060418004.nex';'P060518001.nex';'P060518002.nex';'P060518003.nex';'P060518004.nex';'P060518005.nex';'P060518006.nex';'P060518007.nex';'P060518008.nex';'P060518009.nex';'P060518010.nex';'P060518011.nex';'P060518012.nex';'P060518013.nex'};
RWmultiplier1 = [1.0;1.0;1.0;1.0;1.0;0.7;0.7;0.4;0.4;0.2;0.2;0.3;0.5;0.5;0.8;0.8;1.0];
filename2 = {'P060618001.nex';'P060618002.nex';'P060618003.nex';'P060618004.nex';'P060618005.nex';'P060618006.nex';'P060618007.nex';'P060618008.nex';'P060618009.nex';'P060618010.nex';'P060618011.nex';'P060618012.nex';'P060618013.nex'};
RWmultiplier2 = [1.0;0.5;0.5;0.0;0.25;0.25;0.25;0.60;0.60;0.35;0.35;0.10;0.10];
filename3 = {'P060718001.nex';'P060718002.nex';'P060718003.nex';'P060718004.nex';'P060718005.nex'};
RWmultiplier3 = [0.3;0.3;0.3;0.2;0.2];
filename4 = {'P060818001.nex';'P060818002.nex';'P060818003.nex';'P060818004.nex';'P060818005.nex';'P060818006.nex';'P060818007.nex';'P060818008.nex';'P060818009.nex';'P060818010.nex';'P060818011.nex';'P060818012.nex';'P060818013.nex'};
RWmultiplier4 = [1.0;0.6;0.4;0.2;0.8;0.1;0.0;0.5;0.3;0.7;0.9;1.0;1.0];
filename5 = {'P061118001.nex';'P061118002.nex';'P061118003.nex';'P061118004.nex';'P061118005.nex';'P061118006.nex';'P061118007.nex'};
RWmultiplier5 = [0.5;0.8;0.3;0.6;0.1;0.9;0.05]; 
filename6 = {'P061218001.nex';'P061218002.nex';'P061218003.nex';'P061218004.nex';'P061218005.nex';'P061218006.nex';'P061218007.nex';'P061218008.nex';'P061218009.nex';'P061218010.nex'};
RWmultiplier6 = [1.0;0.5;0.5;0.3;0.3;0.3;0.2;0.2;0.2;0.2];
% Before the nblocks and trials/blocks were changed, Now its just 1 block, previously it used to be 30 blocks
filename7 = {'P061318001.nex';'P061318002.nex';'P061318003.nex';'P061318004.nex';'P061318005.nex';'P061318006.nex';'P061318007.nex';'P061318008.nex';'P061318009.nex';'P061318010.nex';'P061318011.nex'};
RWmultiplier7 = [0.8;0.6;0.6;0.6;0.4;0.4;0.4;0.2;0.2;0.2;0.2];
filename8 = {'P061418002.nex';'P061418003.nex';'P061418004.nex';'P061418005.nex';'P061418006.nex';'P061418008.nex';'P061418009.nex';'P061418010.nex';'P061418011.nex';'P061418012.nex'};
RWmultiplier8 = [0.8;0.6;0.6;0.4;0.4;0.4;0.3;0.3;0.3;0.2;0.2];
filename9 = {'P061518001.nex';'P061518002.nex';'P061518003.nex';'P061518004.nex';'P061518005.nex'};
RWmultiplier9 = [1.0;1.0;1.0;1.0;1.0];
filename10 = {'P061818002.nex';'P061818003.nex';'P061818004.nex';'P061818005.nex';'P061818006.nex';'P061818007.nex'};
RWmultiplier10 = [1.0;1.0;1.0;1.0;1.0;1.0];
filename11 = {'P061918001.nex';'P061918002.nex';'P061918003.nex';'P061918004.nex';'P061918005.nex';'P061918006.nex'};
RWmultiplier11 = [1.0;1.0;1.0;1.0;1.0;1.0];
filename12 = {'P062018001.nex';'P062018002.nex';'P062018003.nex';'P062018004.nex';'P062018005.nex'};
RWmultiplier12 = [1.0;1.0;1.0;1.0;1.0];
filename13 = {'P062118001.nex';'P062118002.nex';'P062118003.nex';'P062118004.nex';'P062118005.nex'};
RWmultiplier13 = [1.0;1.0;1.0;1.0;1.0];
filename14 = {'P062518001.nex';'P062518002.nex';'P062518003.nex'};
RWmultiplier14 = [1.0;1.0;1.0];
filename15 = {'P062918001.nex';'P062918002.nex';'P062918003.nex';'P062918004.nex';'P062918005.nex';'P062918006.nex';'P062918007.nex';'P062918008.nex'};
RWmultiplier15 = [1.0;1.0;1.0;1.0;1.0;1.0;1.0;1.0];
filename16 = {'P070218001.nex';'P070218002.nex';'P070218003.nex';'P070218004.nex';'P070218005.nex';'P070218006.nex';'P070218007.nex';'P070218008.nex';'P070218009.nex';'P070218010.nex'};
RWmultiplier16 = [1.0;1.0;1.0;1.0;1.0;1.0;0.4;0.25;0.20;0.15];
filename17 = {'P070318001.nex';'P070318002.nex';'P070318003.nex';'P070318004.nex';'P070318005.nex';'P070318006.nex';'P070318007.nex';'P070318008.nex';'P070318009.nex';'P070318010.nex';'P070318011.nex';'P070318012.nex'};
RWmultiplier17 = [1.0;1.0;0.5;0.3;0.3;0.2;0.2;0.15;0.15;0.1;0.1;0.07];
filename18 = {'P070518002.nex';'P070518003.nex';'P070518004.nex';'P070518005.nex';'P070518006.nex';'P070518007.nex';'P070518008.nex'};
RWmultiplier18 = [1.0;0.5;0.5;0.5;0.5;0.3;0.2];
filename19 = {'P070618001.nex';'P070618002.nex';'P070618003.nex';'P070618004.nex';'P070618006.nex';'P070618008.nex';'P0706180010.nex'};
RWmultiplier19 = [1.0;1.0;0.5;0.5;0.3;1.0;1.0;1.0];
filename20 = {'P070918001.nex';'P070918003.nex';'P070918004.nex';'P070918006.nex';'P070918007.nex';'P070918009.nex';'P070918010.nex';'P070918011.nex'};
RWmultiplier20 = [1.0;1.0;0.5;0.5;0.5;0.25;0.25;0.25];
filename21 = {'P071018001.nex';'P071018002.nex';'P071018003.nex';'P071018004.nex';'P071018005.nex'};
RWmultiplier21 = [1.0;1.0;0.3;0.2;1.0];
filename22 = {'P071118001.nex';'P071118002.nex';'P071118003.nex'};
RWmultiplier22 = [1.0;1.0;1.0];
filename23 = {'P071218001.nex';'P071218002.nex';'P071218003.nex';'P071218004.nex';'P071218005.nex'};
RWmultiplier23 = [1.0;1.0;1.0;1.0;1.0];
filename24 = {'P071618001.nex';'P071618002.nex';'P071618003.nex';'P071618004.nex';'P071618005.nex';'P071618006.nex';'P071618007.nex';'P071618008.nex';'P071618009.nex';'P071618010.nex'};
RWmultiplier24 = [1.0;1.0;1.0;1.0;1.0;1.0;1.0;1.0;1.0;1.0];
filename25 = {'P071718001.nex';'P071718002.nex';'P071718003.nex';'P071718004.nex';'P071718005.nex';'P071718006.nex';'P071718007.nex';'P071718008.nex'};
RWmultiplier25 = [1.0;1.0;1.0;1.0;1.0;1.0;1.0;1.0];
filename26 = {'P071918001.nex';'P071918002.nex';'P071918003.nex';'P071918004.nex';'P071918005.nex';'P071918006.nex'};
RWmultiplier26 = [1.0;1.0;1.0;1.0;1.0;1.0];
filename27 = {'P072018001.nex';'P072018002.nex';'P072018003.nex';'P072018004.nex';'P072018005.nex';'P072018006.nex';'P072018007.nex';'P072018008.nex'};
RWmultiplier27 = [1.0;1.0;1.0;1.0;1.0;1.0;1.0;0.2];
% filenamebeh = [filename1; filename2; filename3; filename4; filename5; filename6; filename7; filename8; filename9; filename10; filename11; filename12; filename13; filename14];
% RWmultiplierbeh = [RWmultiplier1; RWmultiplier2; RWmultiplier3; RWmultiplier4; RWmultiplier5; RWmultiplier6; RWmultiplier7; RWmultiplier8; RWmultiplier9; RWmultiplier10; RWmultiplier11; RWmultiplier12; RWmultiplier13; RWmultiplier14];
filenamebeh = [filename22; filename23; filename24; filename25; filename26;  filename27];
RWmultiplierbeh = [RWmultiplier22; RWmultiplier23; RWmultiplier24; RWmultiplier25; RWmultiplier26; RWmultiplier27];
% filenamebeh = [filename9; filename10; filename11;filename12;filename13;filename14];
% RWmultiplierbeh = [RWmultiplier9; RWmultiplier10; RWmultiplier11; RWmultiplier12; RWmultiplier13; RWmultiplier14];
laserdialbeh = zeros(size(RWmultiplierbeh));
save filenamebehP filenamebeh
save RWmultiplierbehP RWmultiplierbeh
save laserdialbehP laserdialbeh




