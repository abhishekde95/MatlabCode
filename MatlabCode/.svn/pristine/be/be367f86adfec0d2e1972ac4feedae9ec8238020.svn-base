% This code sets up the analysis GUI for GridLMSubunit datafiles.
% 3/22/13       Created.        JPW
global GLMP DN

if ismac
    library = '/Users/jpatrickweller/Dropbox/Patrick/GLMS Data/nex files/';
elseif ispc
    %library = 'C:\Users\jpweller\Dropbox\Patrick\GLMS Data\';
    library = 'C:\Documents and Settings\JPatrickWeller\My Documents\Dropbox\Patrick\GLMS Data\nex files\';
    %library = 'C:\PlexonData\Old nex files\Old plx files\';
    %library = 'N:\PlexonData\Patrick\'
end

% GLMS v1.0
datafile = 'N060613003.nex'; % +L-M
%datafile = 'N060913002.nex'; % +L-M 
%datafile = 'N061113001.nex'; % +L-M (sub2)
%datafile = 'N061113003.nex'; % +L-M  
%datafile = 'N061213002.nex'; % +L-M (sub1)
%datafile = 'N061313001.nex'; % +L-M (sub 1)
%datafile = 'N061513002.nex'; % -L+M Responsive (#34)
%datafile = 'N061513001.nex'; % -L+M Responsive (#34)
%datafile = 'N062313002.nex'; % RG Chromaic (center and surround, both chromatic, may be one subunit) (sub3 #35)
%datafile = 'N062713004.nex'; % RG Chromatic (both +L-M and-L+M) (#36)
%datafile = 'N062813002.nex'; % +L-M Chromatic (tested 2 possible subunits, only 1 looks real) (sub2 #37)
%datafile = 'N062913003.nex'; % -L+M Chromatic (#38)
%datafile = 'N062913004.nex'; % -L+M Responsive (#39)
%datafile = 'N070213002.nex'; % +L-M Chromatic (tested 2 possible subunits, only 1 looks real) (sub2 #40)
%datafile = 'N070313003.nex'; % -L+M Chromatic (#41)
%datafile = 'N061513002.nex'; % -L+M 
%datafile = 'N061813001.nex'; % +L-M  
%datafile = 'N062313002.nex'; % +L-M (sub1)
%datafile = 'N062713004.nex'; % +L-M and -L+M (not in pop bc no RF in DN)
%datafile = 'N062813002.nex'; % +L-M and -L+M (sub1 and sub2)
%datafile = 'N062913003.nex'; % -L+M
%datafile = 'N062913004.nex'; % -L+M
%datafile = 'N070213002.nex'; % +L-M (sub 2)
%datafile = 'N070313001.nex'; % +L-M (sub 1) (DN looks weak)
%datafile = 'N070313003.nex'; % -L+M

% GridLMSubunit_StimControl Datasets GLMS v1.1
%datafile = 'N100313002.nex'; % +L-M (sub 3)
%datafile = 'N100713001.nex'; % -L+M 
%datafile = 'N102413002.nex'; % +L-M
%datafile = 'N111713001.nex'; % +L-M
%datafile = 'N112113002.nex'; % +L-M
%datafile = 'N112513002.nex' % +L-M (used high contrast lum WN, so this is an L-cell)
%datafile = 'N112613001.nex'; % +L-M
%datafile = 'N112713002.nex'; % +L-M
%datafile = 'N112713003.nex'; % +L-M
%datafile = 'N112713005.nex'; % +L-M

% Added special high luminance case GLMS v1.2
%datafile = 'N120513003.nex'; % Horseshoe
%datafile = 'N120613002.nex'; % Horseshoe
%datafile = 'N120713001.nex'; % +L+M and -L-M
%datafile = 'N121613001.nex'; % +L-M

% GLMSv2.0 (Incorporated higher luminance WN and GLMP)
%datafile = 'N031814001.nex'; % horseshoe? show greg.
%datafile = 'N031914001.nex'; % +L-M
%datafile = 'N032114001.nex'; % -L-M
%datafile = 'N032514001.nex'; % -L-M and +L+M
%datafile = 'N032614001.nex'; % +/- luminance (#59)

% Nut's second chamber (GLMSv2.0 continued)
%datafile = 'N102114003.nex'; % +L+M and -L-M (sub1 and sub2)
%datafile = 'N102214002.nex'; % +L+M and -L-M
%datafile = 'N102614001.nex'; % +L+M and -L-M 
%datafile = 'N102714002.nex'; % +L+M and -L-M
%datafile = 'N103014006.nex'; % +L+M and -L-M
%datafile = 'N103114002.nex'; % +L+M and -L-M (maybe chuck?)
%datafile = 'N031914001.nex'; % +L-M
%datafile = 'N110514001.nex'; % +L+M and -L-M
%datafile = 'N110514002.nex'; % add to 001
%datafile = 'N111214001.nex'; % merged with 002 into 003
%datafile = 'N111214002.nex'; % merged with 001 into 003
%datafile = 'N111214003.nex'; % Horseshoe (-L+M)
%datafile = 'N111414001.nex'; % +L+M and -L-M
%datafile = 'N111414002.nex'; % 2 luminance subunits, responses look identical (#2 in pop)
%datafile = 'N111714002.nex'; % 2 luminance subunits (using both in pop)
%datafile = 'N112014002.nex'; % bipolar luminance
%datafile = 'N112014005.nex'; % bipolar luminance
%datafile = 'N112514001.nex'; % Horseshoe (-L+M)
%datafile = 'N112614001.nex'; % horseshoe? pan color?
%datafile = 'N120314003.nex'; % ** DS cell ** (not in pop)
%datafile = 'N120414001.nex'; % ** DS cell ** (not in pop)
%datafile = 'N120414004.nex'; % ** DS cell ** (not in pop)
%datafile = 'N120414005.nex'; % +L+M and -L-M (sub1 and sub2)
%datafile = 'N120514002.nex'; % -L-M and +L+M
%datafile = 'N120714004.nex'; % +L+M and -L-M (sub1 and sub2)
%datafile = 'N120914001.nex'; % +L+M and -L-M (non stationarity!)
%datafile = 'N121714001.nex'; % +L+M and -L-M
%datafile = 'N121914001.nex'; % Horseshoe cell (+L-M)
%datafile = 'N121914002.nex'; % +L+M and -L-M
%datafile = 'N122114001.nex'
%datafile = 'N122314002.nex'; % +L+M and -L-M (non stationarity!)
%datafile = 'N010515001.nex'; %-L-M (and a bit of +L-M)
%datafile = 'N010515003.nex'; % Horseshoe/ L-cell (-L-M and -L+M)
%datafile = 'N010615001.nex'; % -L-M (bad DN and small responses)
%datafile = 'N010715002.nex'; % -L-M
%datafile = 'N010715003.nex'; % same cell as 002. 2 subunits, but only 2nd is correctly aligned (not in pop)
%datafile = 'N010815002.nex'; % crazy cell. WN: -L+M; GLMPL: suppression to high lum contrast
%datafile = 'N010915002.nex'; % Bipolar luminance
%datafile = 'N011115001.nex'; % Bipolar luminance
%datafile = 'N011115002.nex'; % Bipolar luminance (weak WN)
%datafile = 'N011215001.nex'; % bipolar luminance
%datafile = 'N011215002.nex'; % Horseshoe cell (+L+M)
%datafile = 'N011315001.nex'; % bipolar luminance (unselected surround, or 2nd subunit)
%datafile = 'N012115001.nex'; % -L+M and gabor
%datafile = 'N020215001.nex'; % 2 lum subunits, also Gabor
%datafile = 'N020215002.nex'; % 2 lum subunits, very little data. Also gabor. (sub1 #67. sub2 is strange fit)
%datafile = 'N020515001.nex'; % 2 lum subunits
%datafile = 'N020515002.nex'; % L-cell (-L+M and +L+M). Also Gabor.
%datafile = 'N021115001.nex'; % High baseline (not in pop) ** alternate fit candidate **
%datafile = 'N021115003.nex'; % High baseline (not in pop) ** alternate fit candidate **
%datafile = 'N021915002.nex'; % High baseline (not in pop) ** alternate fit candidate **
%datafile = 'N022015003.nex'; % -L-M cell
%datafile = 'N022015004.nex'; % 2 luminance subunits
%datafile = 'N022315002.nex'; % pan color? sloppy iso
%datafile = 'N022315004.nex'; % -L-M
%datafile = 'N022315005.nex'; % Simple cell, 2 lum subunits. (#73 & # 74)
datafile = 'N022415001.nex'; % doubel opponent cell
%datafile = 'N022415002.nex'; % doubel opponent cell *** missing? ****
%datafile = 'N022715001.nex'; % +L+M. Also gabor. (#77)
%datafile = 'N030215001.nex'; % +L-M (#78) ** good latency file
%datafile = 'N030415001.nex'; % re-examine later, DS cell?
%datafile = 'N032415001.nex'; % Bipolar luminance, +Gabor *good latency file
%datafile = 'N032415002.nex'; % Bipolar lum
%datafile = 'N032815001.nex'; % Bipolar luminance
%datafile = 'N032915001.nex'; % High baseline (not in pop) ** alternate fit candidate **
%datafile = 'N033015001.nex'; % Horseshoe (but barely)
%datafile = 'N033115001.nex'; % pan color
%datafile = 'N033115002.nex'; % L-cell
%datafile = 'N040215001.nex'; % Horseshoe (-L+M)
%datafile = 'N040215002.nex'; % L-cell
%datafile = 'N040315002.nex'; % bipolar luminance? Maybe M-cone responsive.
%datafile = 'N040515001.nex'; % bipolar luminance
%datafile = 'N040615003.nex'; % bipolar lum
%datafile = 'N040815003.nex'; % +L+M
%datafile = 'N040815004.nex'; % +L+M
%datafile = 'N040915001.nex'; % +L-M
%datafile = 'N040915002.nex'; % +lum and +L+M (#91)
%datafile = 'N041315002.nex'; % -L+M (#92)
%datafile = 'N041315004.nex'; % pan color? Horseshoe? (#93)
%datafile = 'N041415001.nex'; % simple cell, 2 subunits! (#94 & #95) Also Gabor.
%datafile = 'N041415002.nex'; % Simple cell, 2 subunits, but too little data (not thru rnd 1). Also Gabor.
%datafile = 'N041615001.nex'; % Small RF. Two lum subunits, but almost overlapping.
%datafile = 'N041715002.nex'; % Sub 1 looks good. Lum.
%datafile = 'N042215001.nex'; % lum? maybe pan color?
%datafile = 'N042615001.nex'; % Pan color
%datafile = 'N042715001.nex'; % too many headers?
%datafile = 'N051515002.nex'; % Lum (maybe pan color?) + gabor
%datafile = 'N051915001.nex'; % looks random? Take out of pop
%datafile = 'N052115001.nex'; % bipolar lum + gabor
%datafile = 'N052415001.nex'; % lum (maybe some chrom - check iso) +gabor
%datafile = 'N052715001.nex'; % simple lum cell + gabor
%datafile = 'N052915001.nex'; % symmetric simple lum cell! + gabor
%datafile = 'N052915002.nex'; % pan color?
%datafile = 'N060415002.nex'; % bipolar lum
%datafile = 'N060915001.nex'; % bipolar lum

% Maui Datafiles (these first ones are on the bad monitor)
%datafile = 'M061515001.nex'; % luminance! nice off resp!
%datafile = 'M061615001.nex'; % luminance - STA is unclear.
%datafile = 'M061715001.nex'; % maybe DS, no structure, large (but unclear) RF
%datafile = 'M062015002.nex'; % WN looks unclear, GLMP -L+M horseshoe
%datafile = 'M062015003.nex'; % -L-M excitatory, +L-M inhibitory
%datafile = 'M062115001.nex'; % unclear WN, lum GLMP
%datafile = 'M062115002.nex'; % -L-M
%datafile = 'M062215001.nex'; % 2 subunits, lum + chrom? Great iso.
%datafile = 'M062315003.nex'; % Nice lum STA, some chrom resps in GLMP

% ------ New monitor ------ % (Probably the datasets above are bad-ish)
% **** These datafiles are using the wrong gamma table *****
%datafile = 'M062715002.nex'; % Horseshoe!
%datafile = 'M070515002.nex'; % Luminance! -L-M Only got one subunit (of 2)
%datafile = 'M070715001.nex'; % +L-M chromatic!
%datafile = 'M071415001.nex'; % -L-M 
%datafile = 'M071615002.nex'; % -L-M, but 3.3 sp max, should it go in pop? (not in pop)
%datafile = 'M071915001.nex'; % Strange - looks like a nonlinear off cell (not in pop)
%datafile = 'M072215001.nex'; % Small -L-M resp
%datafile = 'M072315003.nex'; % Two subunits. one horseshoe, one pan color...
%datafile = 'M072715001.nex'; % lum (maybe pan color)
%datafile = 'M072815002.nex'; % pancolor
%datafile = 'M072915001.nex'; % bipolar lum
%datafile = 'M073015001.nex'; % pan color
%datafile = 'M073015002.nex'; % (+L-M) Combine with 003
%datafile = 'M073015003.nex'; % Combine with 002
%datafile = 'M073115001.nex'; % bipolar lum
%datafile = 'M073115002.nex'; % lum, maybe pan color
%datafile = 'M080115001.nex'; % -L-M lum
%datafile = 'M080315001.nex'; % +L-M chromatic
%datafile = 'M080415001.nex'; % pan color
%datafile = 'M080415002.nex'; % lum and some -L+M
%datafile = 'M080415003.nex'; % lum? lots of nice spikes, but very little signal
%datafile = 'M080415004.nex'; % lum (small dataset)
%datafile = 'M080515002.nex'; % lum
%datafile = 'M080715001.nex'; % lum
%datafile = 'M080815001.nex'; % Horseshoe
%datafile = 'M081315001.nex'; % +L+M
%datafile = 'M081715003.nex'; % pan color
%datafile = 'M081815002.nex'; % +L+M
%datafile = 'M081815003.nex'; % -L-M
%datafile = 'M081915001.nex'; % -L-M
%datafile = 'M081915002.nex'; % double opponent! 2 subs
%datafile = 'M082415001.nex'; % simple lum! w subs
%datafile = 'M082515001.nex'; % bipolar lum
%datafile = 'M082615002.nex'; % pan color
%datafile = 'M082615003.nex'; % not great
%datafile = 'M090315002.nex'; % Bracket cell? -L-M and +L-M
%datafile = 'M090415001.nex'; % -L-M

% *** Dell 7, correct gamma table
% datafile = 'M012516001.nex'; % looks kinda rando
% datafile = 'M012716004.nex'; % also kinda rando
% datafile = 'M012716005.nex'; % pan color (mostly lum) + gabor
% datafile = 'M012816002.nex'; % simple cell DN, horseshoe GLMP? use both subs??
% datafile = 'M013016002.nex'; % weird
% datafile = 'M013116001.nex'; % bipolar lum
% datafile = 'M013116002.nex'; % bipolar lum, maybe horseshoe
% datafile = 'M020316001.nex'; % bipolar lum
% datafile = 'M020616003.nex'; % -L-M
 datafile = 'M020816001.nex'; % pan color
datafile = 'M020816002.nex'; % pan color? maybe chuck
datafile = 'M020916001.nex'; % -L-M
% datafile = 'M020916002.nex'; % bipolar lum?
% datafile = 'M021016001.nex' % simple lum cell. 2 subs
% datafile = 'M021116001.nex'; % High BL, suppression to +lum (not in pop)
% datafile = 'M030816002.nex'; % high BL, no modulation
% datafile = 'M031016001.nex'; % Pan color
% datafile = 'M041816003.nex'; % -L+M chromatic! Finally!
% datafile = 'M042716001.nex'; % bipolar lum + Gabor
% datafile = 'M060716001.nex'; % Two simple cell lum subunits. Add to pop?
% datafile = 'M061516003.nex'; % bad isolation. Weird responses. Not in pop.
% datafile = 'M062116001.nex'; % +L-M
% datafile = 'M070616001.nex'; % Bracket cell (-lum and red chrom) 1,2,3 are all same cell, probably.
% datafile = 'M070616002.nex'; % Bracket cell (-lum and red chrom)
% datafile = 'M070616003.nex'; % Bracket Cell (-lum and red chrom)
% datafile = 'M071516001.nex'; % Bipolar lum 
% datafile = 'M072916003.nex'; % Horshoe cell (-L+M)
% datafile = 'M080716001.nex'; % +L-M
% datafile = 'M081816003.nex'; % -L+M
% datafile = 'M082616001.nex'; 
% datafile = 'M090516001.nex'; % +L-M
% datafile = 'M090616002.nex'; 

try
%     rawdata1 = nex2stro([char(library) datafile1]);    
%     rawdata2 = nex2stro([char(library) datafile2]);
%     rawdata3 = nex2stro([char(library) datafile3]);
%     rawdata = strocat(rawdata1,rawdata2,rawdata3);
    rawdata = nex2stro([char(library) datafile]);  
catch
    rawdata = nex2stro;
end

disp(['Now processing ' datafile])

% Organize Rawdata
[GLMP,DN] = OrganizeRawGLMSData(rawdata);


% Analysis GUI (interactive data display and analysis)
GLMSGUI_Overview();
GLMS_AnalysisGUI(GLMP,DN);

load gong
sound(y,Fs)
clear y
clear Fs

