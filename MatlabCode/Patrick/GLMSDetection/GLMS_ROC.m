% Ideal Observer Analysis for GLMS Datasets
function params = GLMS_ROC(glmp,sub)
global GLMP
GLMP = glmp;

%datafile = 'N060613003.nex'; % +L-M Responsive (#29)
%datafile = 'N060913002.nex'; % +L-M Chromatic, with small -L+M response (Very large dataset) (#30) *
%datafile = 'N061113001.nex'; % +L-M Chromatic, with smaller -L+M response (tested 2 possible subunits, only 1 looks real) (sub2 #31)
%datafile = 'N061113003.nex'; % +L-M responsive (#32)
%datafile = 'N061213002.nex'; % RG Chromatic (tested 2 possible subunits, only 1 looks real) (sub1 #33)
%datafile = 'N061513002.nex'; % -L+M Responsive (#34) *
%datafile = 'N062313002.nex'; % RG Chromaic (center and surround, both chromatic, may be one subunit) (sub3 #35)
%datafile = 'N062713004.nex'; % RG Chromatic (both +L-M and-L+M) (#36) *
%datafile = 'N062813002.nex'; % +L-M Chromatic (tested 2 possible subunits, only 1 looks real) (sub2 #37)
%datafile = 'N062913003.nex'; % -L+M Chromatic (#38)
%datafile = 'N062913004.nex'; % -L+M Responsive (#39) *
%datafile = 'N070213002.nex'; % +L-M Chromatic (tested 2 possible subunits, only 1 looks real) (sub2 #40) *
%datafile = 'N070313003.nex'; % -L+M Chromatic (#41) *

% GridLMSubunit_StimControl Datasets GLMS v1.1
%datafile = 'N100313002.nex'; %+L-M (add with 001) (sub1+2 is #42) ** trim 001 header
%datafile = 'F100513002.nex'; % very nice -L-M DN, ambiguous GLMP
%datafile = 'F100513003.nex'; % Same clear subunit! Cannot combine bc GLMP RF's are different
%datafile = 'N101113001.nex' - trim header, look at tuning

%datafile = 'N102413002.nex'; % +L-M Chromatic (#43) *
%datafile = 'N111313003.nex'; % DN lacking, but clear +L-M GLMP (#44) - maybe chuck
%datafile = 'N111713001.nex'; % +L-M Chromatic (#45) *
%datafile = 'N111713003.nex'; % -L-M Chromatic (#46)
%datafile = 'N112113002.nex'; % +L-M Chromatic (noisy) (#47)
%datafile = 'N112513001.nex' % +L-M Chromatic (# assign number)
%datafile = 'N112513002.nex' % +L-M Chromatic (# assign number) - same as 001?

%datafile = 'N112613001.nex'; % +L-M Chromatic (#48)
%datafile = 'N112713002.nex'; % +L-M Chromatic (#49) *
%datafile = 'N112713003.nex'; % Pan Color (DN: -L-M; GLMP: +L-M) (#50)
%datafile = 'N112713005.nex'; % +L-M Chromatic (#51) *

% Added special high luminance case GLMS v1.2
%datafile = 'N120613002.nex'; % "Horseshoe Cell" (+/- lum and +L-M) (#52)
%datafile = 'N120713001.nex'; % Simple luminance cell! 2 subunits! (#53 & #54)
%datafile = 'N121613001.nex'; % +L-M Chromatic (#55)

% GLMSv2.0 (Incorporated higher luminance WN and GLMP)
%datafile = 'N031814001.nex'; % +/- luminance (#56) *
%datafile = 'N031914001.nex'; % +L-M Chromatic (#57) *
%datafile = 'N032514001.nex'; % +/- luminance (very little data) (#58)
%datafile = 'N032614001.nex'; % +/- luminance (#59)

if ismac
    library = '/Users/jpatrickweller/Dropbox/Patrick/GLMS Data/';
elseif ispc
    library = 'C:\Users\jpweller\Dropbox\Patrick\GLMS Data\';
end

% datafile = [GLMP.datafile '.nex'];
% rawdata = nex2stro([char(library) datafile]);
% [GLMP,~] = OrganizeRawGLMSData(rawdata);

%% ROC

AUC = nan(numel(GLMP.subunit{sub}.uniqueLcc),1);

for s = 1:numel(GLMP.subunit{sub}.uniqueLcc) 
    
    % Define terms
%     spHist = [];
%     blHist = [];

%     for t = 1:numel(GLMP.subunit{sub}.spiketimes_col(s,:))
%         if ~isempty(GLMP.subunit{sub}.spiketimes_col(s,t))
%             spHist = [spHist numel(GLMP.subunit{sub}.spiketimes_col{s,t})];
%         end
%     end
    
    idx = GLMP.subunit{sub}.uniqueIdx{s};
    spHist = GLMP.subunit{sub}.fr(idx)';
    %blHist = GLMP.subunit{sub}.blfr(idx)';
    blHist = GLMP.subunit{sub}.blfr';

    %blHist = randi([20 30],size(spHist));
    %blIdx = randi([1 numel(GLMP.fpon)],size(spHist));
    %blIdx = 
    %blHist = GLMP.blnspikes(blIdx)';
    %blHist = GLMP.blfr(blIdx)';
    thresh = unique([spHist blHist]);
    TPR = nan(1,length(thresh));
    FPR = TPR;

    for n=1:length(thresh)
        TPR(n) = sum(spHist >= thresh(n))/length(spHist);
        FPR(n) = sum(blHist >= thresh(n))/length(blHist);
    end
    
    AUC(s) = trapz(fliplr(FPR),fliplr(TPR));
    %keyboard
    
end

x = GLMP.subunit{sub}.uniqueLcc;
y = GLMP.subunit{sub}.uniqueMcc;
z =  AUC;
F = TriScatteredInterp(x,y,z);
[qx,qy] = meshgrid(min(x):.01:max(x),min(y):.01:max(y));


figure(4); clf; hold on;
plot3(x,y,AUC,'k*')
surfc(qx,qy,F(qx,qy))

%% 1D fit

% Set up some variables
vlb = [.0001  .005  .005    2   .5  -pi 0];
vub = [1       50    50     5   .5   pi 0];
options = optimset('MaxFunEvals',300,'MaxIter',300,'TolFun',10^-6,'TolX',10^-6,'Display','off','Algorithm','active-set');
topfr = 1;
sigmaguess = 0.1;
expguess = 2;
baselineguess = .5;
angguess = linspace(-pi,pi,9);

resp = AUC;
LMpts = [GLMP.subunit{sub}.uniqueLcc GLMP.subunit{sub}.uniqueMcc];
f1 = nan(numel(angguess),numel(vlb));
fval = nan(numel(angguess),1);

for n = 1:numel(angguess)
       
    paramsGuess = [topfr, sigmaguess, sigmaguess,...
        expguess,baselineguess,angguess(n) 0];
    
    [f1(n,:),fval(n)] = fmincon('FitNakaRushtonFunJPW',paramsGuess,[],[],[],[],vlb,vub,[],options,...
        LMpts,resp,'surface7','bernoulli');
    
end

[~,best] = max(fval);
params = f1(best,:);
fval = fval(best);



