% Statistics

global GLMSPopData

% Grab saved population data
if ismac
    library = '/Users/jpatrickweller/Dropbox/Patrick/GLMS Data/';
elseif ispc
    library = 'C:\Users\jpweller\Dropbox\Patrick\GLMS Data\';
end
load([library 'GLMSPopData.mat'])
datatypes = GLMSPopData(1,:);






%% Average number of stimulus presentations


ParamsFig = get(552,'userdata');
conpanel = get(ParamsFig.conpanel,'userdata');
poppanel = get(ParamsFig.poppanel,'userdata');

popidx = find(~poppanel.excludeL)+1;
npres = nan(numel(popidx),1);

for n = 1:numel(popidx)
   
    GLMP = GLMSPopData{popidx(n),strcmp(datatypes,'GLMP')};
    sub = GLMSPopData{popidx(n),strcmp(datatypes,'Subunit')};
    
    npres(n) = size(GLMP.subunit{sub}.spiketimes_col,2);
    
end


%% RF size

%GLMSPopGUI_Params

% Load figure and pop variables
ParamsFig = get(552,'userdata');
conpanel = get(ParamsFig.conpanel,'userdata');
poppanel = get(ParamsFig.poppanel,'userdata');

L = poppanel.oneDL | poppanel.twoDL;

rfs = cat(1,GLMSPopData{2:end,strcmp(datatypes,'RF Params')});
rfrhos = [rfs.rfrho]';

% [maxrho,maxidx] = max(rfrhos(~poppanel.excludeL));
% [minrho,minidx] = min(rfrhos(~poppanel.excludeL));
% meanrho = mean(rfrhos(~poppanel.excludeL));

[maxrho,maxidx] = max(rfrhos(poppanel.oneDL));
[minrho,minidx] = min(rfrhos(poppanel.oneDL));
meanrho = mean(rfrhos);

disp(['Smallest RF rho = ' num2str(minrho), ', idx = ' num2str(minidx)])
disp(['Largest RF rho = ' num2str(maxrho), ', idx = ' num2str(maxidx)])
disp(['Average RF rho = ' num2str(meanrho)])



%% Regression of c50 on stim size and eccentricity


%GLMSPopGUI_Params
ParamsFig = get(552,'userdata');
conpanel = get(ParamsFig.conpanel,'userdata');
poppanel = get(ParamsFig.poppanel,'userdata');
ndL = poppanel.oneDL;
angL = conpanel.fitcard.oneD.angL;

% % Regression of c50 on rfsize and ecc
L = ndL & any(angL(:,1:2),2); % lum
lumpdsigs =  1./poppanel.params(L,2);
lumapdsigs = 1./poppanel.params(L,3);
rfparams = [GLMSPopData{2:end,strcmp(datatypes,'RF Params')}];
lumrfsize = [rfparams(L).sqdva]';
lumrfecc = [rfparams(L).rfrho]';
intercept = ones(size(lumrfecc));
[b,bint,~,~,stats] = regress(lumpdsigs,cat(2,lumrfsize,lumrfecc,intercept),.05)

L = ndL & any(angL(:,3:4),2); % all pts not in wedges
colpdsigs =  1./poppanel.params(L,2); % col
colapdsigs = 1./poppanel.params(L,3);
%rfparams = [GLMSPopData{2:end,strcmp(datatypes,'RF Params')}]';
colrfsize = [rfparams(L).sqdva]';
colrfecc = [rfparams(L).rfrho]';
intercept = ones(size(colrfecc));
[b,bint,~,~,stats] = regress(colpdsigs,cat(2,colrfsize,colrfecc,intercept),.05)

colparatio = colpdsigs./colapdsigs;
lumparatio = lumpdsigs./lumapdsigs;

ranksum(lumparatio,colparatio)

ranksum(lumrfsize,colrfsize)

ranksum(lumrfecc,colrfecc)

%figure(1); clf; hold on; grid on;
%plot(lumrfecc,lumrfsize,'ko')
%plot(colrfecc,colrfsize,'ro')


% Fit chomatic and luminance (together) residuals with a spline
smoothval = .999;
allrfsize = cat(1,lumrfsize,colrfsize);
allpdsigs = cat(1,lumpdsigs,colpdsigs);

allpp = csaps(allrfsize,allpdsigs,smoothval);
lumpp = csaps(lumrfsize,lumpdsigs,smoothval);
colpp = csaps(colrfsize,colpdsigs,smoothval);

figure(1); clf; hold on; grid on; box on;
plot(lumrfsize,lumpdsigs,'ko')
xlabel('stim size')
plot(colrfsize,colpdsigs,'ro')
ylabel('c50 size')
fnplt(allpp,'b')
fnplt(lumpp,'k--')
fnplt(colpp,'r--')

allpred = csaps(allrfsize,allpdsigs,smoothval,allrfsize);
lumpred = csaps(lumrfsize,lumpdsigs,smoothval,lumrfsize);
colpred = csaps(colrfsize,colpdsigs,smoothval,colrfsize);

allsse = sum((allpdsigs - allpred).^2);
lumsse = sum((lumpdsigs - lumpred).^2);
colsse = sum((colpdsigs - colpred).^2);

realsseratio = allsse/(lumsse+colsse);

%% Permute and recalculate ratio
%idvect = zeros(size(allpred));
%idvect(1:numel(lumpred)) = 1;

nperms = 10000;
permsseratio = nan(nperms,1);

for n = 1:nperms
    
    permvect = randperm(numel(allpred));
    permlumrfsize = allrfsize(permvect(1:numel(lumrfsize)));
    permcolrfsize = allrfsize(permvect(numel(lumrfsize)+1:end));
    permlumc50 = allpdsigs(permvect(1:numel(lumrfsize)));
    permcolc50 = allpdsigs(permvect(numel(lumrfsize)+1:end));
    
    permlumpred = csaps(permlumrfsize,permlumc50,smoothval,permlumrfsize);
    permcolpred = csaps(permcolrfsize,permcolc50,smoothval,permcolrfsize);
    
    permlumsse = sum((permlumc50 - permlumpred).^2);
    permcolsse = sum((permcolc50 - permcolpred).^2);
    
    permsseratio(n) = allsse/(permlumsse+permcolsse);
    
end


figure(2); clf; hold on; box on;
hist(permsseratio,100)

pval = sum(permsseratio >= realsseratio)./ nperms


%% Compare pan color sensitivity to red vs green and +lum vs -lum

datatypes = GLMSPopData(1,:);

GLMSPopGUI_Params
ParamsFig = get(552,'userdata');
conpanel = get(ParamsFig.conpanel,'userdata');
poppanel = get(ParamsFig.poppanel,'userdata');
ndL = poppanel.twoDL;
whichL = poppanel.twoD.bichrom.eliL;
L = ndL & whichL;




%% For greg, comparing 1d lum bipolar to 2d lum sigs

%GLMSPopGUI_Params
ParamsFig = get(552,'userdata');
conpanel = get(ParamsFig.conpanel,'userdata');
poppanel = get(ParamsFig.poppanel,'userdata');
ndL = poppanel.oneDL;
allangL = conpanel.fitcard.oneD.angL;
angL = any(allangL(:,1:2),2);
bipolarL = poppanel.oneD.bichrom.L;
L = ndL & angL & bipolarL; % 1D lum bipolar


mean(mean(1./poppanel.params(L,2:3)))


%% Subunits per monkey

L = poppanel.oneDL | poppanel.twoDL;
fns = cat(1,GLMSPopData{find(L)+1,1});

sum(fns(:,1)=='N')
sum(fns(:,1)=='M')



%% Compare 1D and 2D parameters

ParamsFig = get(552,'userdata');
conpanel = get(ParamsFig.conpanel,'userdata');
poppanel = get(ParamsFig.poppanel,'userdata');

params1d = poppanel.params(poppanel.oneDL,:);
params2d = poppanel.params(poppanel.twoDL,:);

avals1d = params1d(:,1);
avals2d = params2d(:,1);
blvals1d = params1d(:,6);
blvals2d = params2d(:,6);
expvals1d = params1d(:,5);
expvals2d = params2d(:,5);
kvals1d = params1d(:,end);
kvals2d = params2d(:,end);

apval = ranksum(avals1d,avals2d)
blpval = ranksum(blvals1d,blvals2d)
expval = ranksum(expvals1d,expvals2d)
kpval = ranksum(kvals1d,kvals2d)

%% Compare lum and col sigmas to pc sigmas

% Compare distributions using statistical tests
angL = conpanel.fitcard.oneD.angL;
L1 = poppanel.oneDL & angL(:,1) & ~poppanel.excludeL;
L2 = poppanel.oneDL & angL(:,2) & ~poppanel.excludeL;
L3 = poppanel.oneDL & angL(:,3) & ~poppanel.excludeL;
L4 = poppanel.oneDL & angL(:,4) & ~poppanel.excludeL;

pcL = poppanel.twoDL & poppanel.twoD.bichrom.eliL;




%% Chromatic vs Lum comparisons

ndL = poppanel.oneDL;
angL = conpanel.fitcard.oneD.angL;


% Compare distributions using statistical tests
L1 = ndL & angL(:,1) & ~poppanel.excludeL;
L2 = ndL & angL(:,2) & ~poppanel.excludeL;
L3 = ndL & angL(:,3) & ~poppanel.excludeL;
L4 = ndL & angL(:,4) & ~poppanel.excludeL;


%%% Compare polarities

upperApLpM = poppanel.params(L1,1);
upperAmLmM = poppanel.params(L2,1);
upperApLmM = poppanel.params(L3,1);
upperAmLpM = poppanel.params(L4,1);

sigmaspLpM = 1./poppanel.params(L1,2);
sigmasmLmM = 1./poppanel.params(L2,2);
sigmaspLmM = 1./poppanel.params(L3,2);
sigmasmLpM = 1./poppanel.params(L4,2);

exppLpM = poppanel.params(L1,5);
expmLmM = poppanel.params(L2,5);
exppLmM = poppanel.params(L3,5);
expmLpM = poppanel.params(L4,5);

blpLpM = poppanel.params(L1,6);
blmLmM = poppanel.params(L2,6);
blpLmM = poppanel.params(L3,6);
blmLpM = poppanel.params(L4,6);

kappapLpM = poppanel.params(L1,end);
kappamLmM = poppanel.params(L2,end);
kappapLmM = poppanel.params(L3,end);
kappamLpM = poppanel.params(L4,end);

% luminance axes
pvalAlum = ranksum(upperApLpM,upperAmLmM);
disp(['Lum upperAs: ' num2str(mean(upperApLpM)) ', ' num2str(mean(upperAmLmM)) ', '...
    'mean = ' num2str(rndofferr(mean(cat(1,upperApLpM,upperAmLmM)),2)) ', '...
    'diff = ' num2str(abs(mean(upperApLpM)-mean(upperAmLmM))) ', '...
    'p = ' num2str(pvalAlum)])

pvalsiglum = ranksum(sigmaspLpM,sigmasmLmM);
disp(['Lum c50s: ' num2str(mean(sigmaspLpM)) ', ' num2str(mean(sigmasmLmM)) ', '...
    'mean = ' num2str(rndofferr(mean(cat(1,sigmaspLpM,sigmasmLmM)),2)) ', '...
    'diff = ' num2str(abs(mean(sigmaspLpM)-mean(sigmasmLmM))) ', '...
    'p = ' num2str(pvalsiglum)])

pvalexplum = ranksum(exppLpM,expmLmM);
disp(['Lum exps:  ' num2str(mean(exppLpM)) ', ' num2str(mean(expmLmM)) ', '...
    'mean = ' num2str(rndofferr(mean(cat(1,exppLpM,expmLmM)),2)) ', '...
    'diff = ' num2str(abs(mean(exppLpM)-mean(expmLmM))) ', '...
    'p = ' num2str(pvalexplum)])

pvalbllum = ranksum(blpLpM,blmLmM);
disp(['Lum bls:  ' num2str(mean(blpLpM)) ', ' num2str(mean(blmLmM)) ', '...
    'mean = ' num2str(rndofferr(mean(cat(1,blpLpM,blmLmM)),2)) ', '...
    'diff = ' num2str(abs(mean(blpLpM)-mean(blmLmM))) ', '...
    'p = ' num2str(pvalbllum)])

pvalkappalum = ranksum(kappapLpM,kappamLmM);
disp(['Lum kappas:  ' num2str(mean(kappapLpM)) ', ' num2str(mean(kappamLmM)) ', '...
    'mean = ' num2str(rndofferr(mean(cat(1,kappapLpM,kappamLmM)),2)) ', '...
    'diff = ' num2str(abs(mean(kappapLpM)-mean(kappamLmM))) ', '...
    'p = ' num2str(pvalkappalum)])

% Chromatic axes
pvalAcol = ranksum(upperApLmM,upperAmLpM);
disp(['Col upperAs: ' num2str(mean(upperApLmM)) ', ' num2str(mean(upperAmLpM)) ', '...
    'mean = ' num2str(rndofferr(mean(cat(1,upperApLmM,upperAmLpM)),2)) ', '...
    'diff = ' num2str(abs(mean(upperApLmM)-mean(upperAmLpM))) ', '...
    'p = ' num2str(pvalAcol)])

pvalsigcol = ranksum(sigmaspLmM,sigmasmLpM);
disp(['Col c50s: ' num2str(mean(sigmaspLmM)) ', ' num2str(mean(sigmasmLpM)) ', '...
    'mean = ' num2str(rndofferr(mean(cat(1,sigmaspLmM,sigmasmLpM)),2)) ', '...
    'diff = ' num2str(abs(mean(sigmaspLmM)-mean(sigmasmLpM))) ', '...
    'p = ' num2str(pvalsigcol)])

pvalexpcol = ranksum(exppLmM,expmLpM);
disp(['Col exps: ' num2str(mean(exppLmM)) ', ' num2str(mean(expmLpM)) ', '...
    'mean = ' num2str(rndofferr(mean(cat(1,exppLmM,expmLpM)),2)) ', '...
    'diff = ' num2str(abs(mean(exppLmM)-mean(expmLpM))) ', '...
    'p = ' num2str(pvalexpcol)])

pvalblcol = ranksum(blpLmM,blmLpM);
disp(['Col bls: ' num2str(mean(blpLmM)) ', ' num2str(mean(blmLpM)) ', '...
    'mean = ' num2str(rndofferr(mean(cat(1,blpLmM,blmLpM)),2)) ', '...
    'diff = ' num2str(abs(mean(blpLmM)-mean(blmLpM))) ', '...
    'p = ' num2str(pvalblcol)])

pvalkappacol = ranksum(kappapLmM,kappamLpM);
disp(['Col kappas: ' num2str(mean(kappapLmM)) ', ' num2str(mean(kappamLpM)) ', '...
    'mean = ' num2str(rndofferr(mean(cat(1,kappapLmM,kappamLpM)),2)) ', '...
    'diff = ' num2str(abs(mean(kappapLmM)-mean(kappamLpM))) ', '...
    'p = ' num2str(pvalkappacol)])


%%% Luminace vs Color

% Lum vs col A
if pvalAlum > .05 && pvalAcol > .05
    lumAs = cat(1,upperApLpM,upperAmLmM);
    colAs = cat(1,upperApLmM,upperAmLpM);
    pvalAcollum = ranksum(lumAs,colAs);
    disp(['Lum vs col upperA: '...
        'means = ' num2str(rndofferr(mean(lumAs),2)) ', '... 
        num2str(rndofferr(mean(colAs),2)) ', '...
        'diff = ' num2str(rndofferr(mean(lumAs)-mean(colAs),2)) ', '...
        'pval = ' num2str(pvalAcollum)])
else
    disp('Cannot compare lum and col upperAs, one or more polarities differ significantly.')
end

% Lum vs col sig
if pvalsiglum > .05 && pvalsigcol > .05
    lumsigs = cat(1,sigmaspLpM,sigmasmLmM);
    colsigs = cat(1,sigmaspLmM,sigmasmLpM);
    pvalsigcollum = ranksum(lumsigs,colsigs);
    disp(['Lum vs col sig: '...
        'means = ' num2str(rndofferr(mean(lumsigs),2)) ', '... 
        num2str(rndofferr(mean(colsigs),2)) ', '...
        'diff = ' num2str(rndofferr(mean(lumsigs)-mean(colsigs),2)) ', '...
        'pval = ' num2str(pvalsigcollum)])
else
    disp('Cannot compare lum and col sigs, one or more polarities differ significantly.')
end

% Lum vs col exp
if pvalexplum > .05 && pvalexpcol > .05
    lumexps = cat(1,exppLpM,expmLmM);
    colexps = cat(1,exppLmM,expmLpM);
    pvalexpcollum = ranksum(lumexps,colexps);
    disp(['Lum vs col exp: '...
        'means = ' num2str(rndofferr(mean(lumexps),2)) ', '... 
        num2str(rndofferr(mean(colexps),2)) ', '...
        'diff = ' num2str(rndofferr(mean(lumexps)-mean(colexps),2)) ', '...
        'pval = ' num2str(pvalexpcollum)])
else
    disp('Cannot compare lum and col exps, one or more polarities differ significantly.')
end

% Lum vs col bl
if pvalbllum > .05 && pvalblcol > .05
    lumbls = cat(1,blpLpM,blmLmM);
    colbls = cat(1,blpLmM,blmLpM);
    pvalblcollum = ranksum(lumbls,colbls);
    disp(['Lum vs col bl: '...
        'means = ' num2str(rndofferr(mean(lumbls),2)) ', '... 
        num2str(rndofferr(mean(colbls),2)) ', '...
        'diff = ' num2str(rndofferr(mean(lumbls)-mean(colbls),2)) ', '...
        'pval = ' num2str(pvalblcollum)])
else
    disp('Cannot compare lum and col bls, one or more polarities differ significantly.')
end

% Lum vs col kappa
if pvalkappalum > .05 && pvalkappacol > .05
    lumkappas = cat(1,kappapLpM,kappamLmM);
    colkappas = cat(1,kappapLmM,kappamLpM);
    pvalkappacollum = ranksum(lumkappas,colkappas);
    disp(['Lum vs col kappa: '...
        'means = ' num2str(rndofferr(mean(lumkappas),2)) ', '... 
        num2str(rndofferr(mean(colkappas),2)) ', '...
        'diff = ' num2str(rndofferr(mean(lumkappas)-mean(colkappas),2)) ', '...
        'pval = ' num2str(pvalkappacollum)])
else
    disp('Cannot compare lum and col kappas, one or more polarities differ significantly.')
end

%%% Compare all 4 directions

upperAs = cat(1,upperApLpM,upperAmLmM,upperApLmM,upperAmLpM);
gr = repmat({'pLpM'},numel(upperApLpM),1);
gr = cat(1,gr,repmat({'mLmM'},numel(upperAmLmM),1));
gr = cat(1,gr,repmat({'pLmM'},numel(upperApLmM),1));
gr = cat(1,gr,repmat({'mLpM'},numel(upperAmLpM),1));
p = kruskalwallis(upperAs,gr,'off');
disp(['All upperAs from same dist: pval = ' num2str(p)])
disp(['Mean upperA = ' num2str(mean(upperAs))])

sigmas = cat(1,sigmaspLpM,sigmasmLmM,sigmaspLmM,sigmasmLpM);
gr = repmat({'pLpM'},numel(sigmaspLpM),1);
gr = cat(1,gr,repmat({'mLmM'},numel(sigmasmLmM),1));
gr = cat(1,gr,repmat({'pLmM'},numel(sigmaspLmM),1));
gr = cat(1,gr,repmat({'mLpM'},numel(sigmasmLpM),1));
p = kruskalwallis(sigmas,gr,'off');
disp(['All sigmas from same dist: pval = ' num2str(p)])
disp(['Mean sigmas = ' num2str(mean(sigmas))])

exps = cat(1,exppLpM,expmLmM,exppLmM,expmLpM);
gr = repmat({'pLpM'},numel(exppLpM),1);
gr = cat(1,gr,repmat({'mLmM'},numel(expmLmM),1));
gr = cat(1,gr,repmat({'pLmM'},numel(exppLmM),1));
gr = cat(1,gr,repmat({'mLpM'},numel(expmLpM),1));
p = kruskalwallis(exps,gr,'off');
disp(['All exps from same dist: pval = ' num2str(p)])
disp(['Mean exp = ' num2str(mean(exps))])


bls = cat(1,blpLpM,blmLmM,blpLmM,blmLpM);
gr = repmat({'pLpM'},numel(blpLpM),1);
gr = cat(1,gr,repmat({'mLmM'},numel(blmLmM),1));
gr = cat(1,gr,repmat({'pLmM'},numel(blpLmM),1));
gr = cat(1,gr,repmat({'mLpM'},numel(blmLpM),1));
p = kruskalwallis(bls,gr,'off');
disp(['All bls from same dist: pval = ' num2str(p)])
disp(['Mean bl = ' num2str(mean(bls))])


kappas = cat(1,kappapLpM,kappamLmM,kappapLmM,kappamLpM);
gr = repmat({'pLpM'},numel(kappapLpM),1);
gr = cat(1,gr,repmat({'mLmM'},numel(kappamLmM),1));
gr = cat(1,gr,repmat({'pLmM'},numel(kappapLmM),1));
gr = cat(1,gr,repmat({'mLpM'},numel(kappamLpM),1));
p = kruskalwallis(kappas,gr,'off');
disp(['All kappas from same dist: pval = ' num2str(p)])
disp(['Mean kappa = ' num2str(mean(kappas))])



%% Figure 4E %%

ndL = poppanel.oneDL;
whichsurfL = poppanel.oneD.bichrom.L;
angL = conpanel.fitcard.oneD.angL;

for n = 1:4
    
    L = ndL & whichsurfL & angL(:,n) & ~poppanel.excludeL;
    angs = poppanel.tuning(L);
    sigmas = 1./poppanel.params(L,2);
    apdsigmas = 1./poppanel.params(L,3);

    disp([cols{n} ' n = ' num2str(sum(L)) ...
        ' Mean PD sig = ' num2str(rndofferr(mean(sigmas),3))... 
        ' Std PD sig = ' num2str(rndofferr(std(sigmas),3))...
        ' Mean APD sig = ' num2str(rndofferr(mean(apdsigmas),3))...
        ' Std APD sig = ' num2str(rndofferr(std(apdsigmas),3))...
        ' n = ' num2str(sum(L))])
end

% Compare distributions using statistical tests
L1 = ndL & whichsurfL & angL(:,1) & ~poppanel.excludeL;
L2 = ndL & whichsurfL & angL(:,2) & ~poppanel.excludeL;
L3 = ndL & whichsurfL & angL(:,3) & ~poppanel.excludeL;
L4 = ndL & whichsurfL & angL(:,4) & ~poppanel.excludeL;
sigmaspLpM = 1./poppanel.params(L1,2);
sigmasmLmM = 1./poppanel.params(L2,2);
sigmaspLmM = 1./poppanel.params(L3,2);
sigmasmLpM = 1./poppanel.params(L4,2);
apdsigmaspLpM = 1./poppanel.params(L1,3);
apdsigmasmLmM = 1./poppanel.params(L2,3);
apdsigmaspLmM = 1./poppanel.params(L3,3);
apdsigmasmLpM = 1./poppanel.params(L4,3);

%%% Luminance %%%
a = L1 | L2;
b = ndL & any(angL(:,1:2),2) & ~poppanel.excludeL;
lumperc = (sum(a) / sum(b))*100;
disp([num2str(rndofferr(lumperc,1)) '% (' num2str(sum(a)) '/' num2str(sum(b)) ') of luminance tuned responses are bipolar.'])

% Luminance confs
% pvalconflum = ranksum(confspLpM,confsmLmM);
% disp(['Lum confs: diff = ' num2str(abs(mean(confspLpM)-mean(confsmLmM))) ' p = ' num2str(pvalconflum)])

% luminance PD sigmas
pvalsiglum = ranksum(sigmaspLpM,sigmasmLmM);
disp(['Lum PD sig: diff = ' num2str(abs(mean(sigmaspLpM)-mean(sigmasmLmM))) ' p = ' num2str(pvalsiglum)])

% Luminance APD sigmas
pvalapdsiglum = ranksum(apdsigmaspLpM,apdsigmasmLmM);
disp(['Lum APD sig: diff = ' num2str(abs(mean(apdsigmaspLpM)-mean(apdsigmasmLmM))) ' p = ' num2str(pvalapdsiglum)])

% PD vs APD
pval = ranksum(sigmaspLpM,apdsigmaspLpM);
disp(['+Lum PD vs APD sig: diff = ' num2str(abs(mean(sigmaspLpM)-mean(apdsigmaspLpM))) ' p = ' num2str(pval)]);
pval = ranksum(sigmasmLmM,apdsigmasmLmM);
disp(['-Lum PD vs APD sig: diff = ' num2str(abs(mean(sigmasmLmM)-mean(apdsigmasmLmM))) ' p = ' num2str(pval)]);


%%% Chromatic %%%
a = L3 | L4;
b = ndL & any(angL(:,3:4),2) & ~poppanel.excludeL;
colperc = sum(a) / sum(b)*100;
disp([num2str(rndofferr(colperc,1)) '% (' num2str(sum(a)) '/' num2str(sum(b))...
    ' of chromatic tuned responses are bipolar.'])

% denomL = ndL & any(angL(:,3:4),2) & ~poppanel.excludeL;
% colperc = sum(L3 | L4) / sum(denomL(31:end))*100;
% disp([num2str(rndofferr(colperc,1)) '% of *adjusted* chromatic tuned responses are bipolar, n = ' num2str(sum(L3|L4))])

% Col confs
% pvalconfcol = ranksum(confspLmM,confsmLpM);
% disp(['Col confs: diff = ' num2str(abs(mean(confspLmM)-mean(confsmLpM))) ' p = ' num2str(pvalconfcol)])

% chrom PD sigmas
pvalsigcol = ranksum(sigmaspLmM,sigmasmLpM);
disp(['Col PD sig: diff = ' num2str(abs(mean(sigmaspLmM)-mean(sigmasmLpM))) ' p = ' num2str(pvalsigcol)])

% Chrom APD sigmas
pvalapdsigcol = ranksum(apdsigmaspLmM,apdsigmasmLpM);
disp(['Col APD sig: diff = ' num2str(abs(mean(apdsigmaspLmM)-mean(apdsigmasmLpM))) ' p = ' num2str(pvalapdsigcol)])

% PD vs APD
pval = ranksum(sigmaspLmM,apdsigmaspLmM);
disp(['red PD vs APD sig: diff = ' num2str(abs(mean(sigmaspLpM)-mean(apdsigmaspLpM))) ' p = ' num2str(pval)]);
pval = ranksum(sigmasmLpM,apdsigmasmLpM);
disp(['green PD vs APD sig: diff = ' num2str(abs(mean(sigmasmLmM)-mean(apdsigmasmLmM))) ' p = ' num2str(pval)]);


%%% Lum vs Col %%%

% % Lum vs Col APD conf
% if pvalconflum > .05 * pvalconfcol > .05
%     lumconfs = cat(1,confspLpM,confsmLmM);
%     colconfs = cat(1,confspLmM,confsmLpM);
%     pvalconfcollum = ranksum(lumconfs,colconfs);
%     disp(['Lum vs col conf: diff = ' num2str(abs(mean(lumconfs)-mean(colconfs))) ' pval = ' num2str(pvalconfcollum)])
% else
%     disp('Cannot compare lum and col axes bc one or more is significantly different...')
% end

% Lum vs Col PD c50
if pvalsiglum > .05 && pvalsigcol > .05
    lumsigs = cat(1,sigmaspLpM,sigmasmLmM);
    colsigs = cat(1,sigmaspLmM,sigmasmLpM);
    pvalsigcollum = ranksum(lumsigs,colsigs);
    disp(['Lum vs Col PD sig: diff = ' num2str(abs(mean(lumsigs)-mean(colsigs))) ' pval = ' num2str(pvalsigcollum)])
end

% Lum vs Col APD c50
if pvalapdsiglum > .05 && pvalapdsigcol > .05
    lumapdsigs = cat(1,apdsigmaspLpM,apdsigmasmLmM);
    colapdsigs = cat(1,apdsigmaspLmM,apdsigmasmLpM);
    pvalsigcollum = ranksum(lumapdsigs,colapdsigs);
    disp(['Lum vs Col APD sig: diff = ' num2str(abs(mean(lumapdsigs)-mean(colapdsigs))) ' pval = ' num2str(pvalsigcollum)])
end

pval = ranksum(lumsigs,lumapdsigs);
disp(['Lum PD vs APD sig: diff = ' num2str(abs(mean(lumsigs)-mean(lumapdsigs))) ' p = ' num2str(pval)]);
pval = ranksum(colsigs,colapdsigs);
disp(['Col PD vs APD sig: diff = ' num2str(abs(mean(colsigs)-mean(colapdsigs))) ' p = ' num2str(pval)]);



%% Are pancolor cells better fit by 2 symmetric axes?
global GLMP

GLMSPopGUI_Params()

% Load figure and pop variables
ParamsFig = get(552,'userdata');
conpanel = get(ParamsFig.conpanel,'userdata');
poppanel = get(ParamsFig.poppanel,'userdata');

% Pull out params
ndL = poppanel.twoDL;
pcL = ndL & poppanel.twoD.bichrom.eliL & ~poppanel.excludeL;
pcidx = find(pcL);

% Fitting stuff
options = optimset('Algorithm','interior-point','MaxFunEvals',5000,...
    'MaxIter',5000,'FinDiffType','central','Hessian','bfgs','display','off',...
    'TolFun',0,'TolCon',0,'FunValCheck','on','AlwaysHonorConstraints','bounds');
surftype = 'conicsection_xy';
errortype = 'NegativeBinomial';
thresh = .075; % For comparing types of fits
params = [];
a = [-1 0 0 0 0 1 0 0];
b = 0;
Aeq = [0 1 -1 0 0 0 0 0];
Beq = 0;
sympdparams = nan(numel(pcidx),numel(ub));
sympdnLL = nan(numel(pcidx),1);
sympdnormLL = nan(numel(pcidx),1);

% Run through each dataset, fit a pan color cell with a symmetric preferred
% axis. Record normalized LL... somewhere. How fast is this?
for n = 1:numel(pcidx)

    GLMP = GLMSPopData{pcidx(n)+1,strcmp(datatypes,'GLMP')};
    sub = GLMSPopData{pcidx(n)+1,strcmp(datatypes,'Subunit')};
    Lcc = cat(1,GLMP.subunit{sub}.Lcc,zeros(size(GLMP.subunit{sub}.blnspikes)));
    Mcc = cat(1,GLMP.subunit{sub}.Mcc,zeros(size(GLMP.subunit{sub}.blnspikes)));
    nsp = cat(1,GLMP.subunit{sub}.nspikes,GLMP.subunit{sub}.blnspikes);
    ub = [max(GLMP.nspikes)             300             300  300 10 max(GLMP.nspikes)  pi 5];
    lb = [min(GLMP.nspikes) 1/max(GLMP.rho) 1/max(GLMP.rho) -300  1              .001 -pi 0];
    
    % Recalculate min and max LL
    kappaguess = poppanel.params(pcidx(n),end);
    
    % Find minimum LL (using Negative Binomial)
    mu = repmat(mean(nsp),size(nsp));
    x = nsp;
    [kappa,f] = fmincon('minimizeKappa',kappaguess,[],[],[],[],0,[],[],options,mu,x);
    lbLL = -f;
    
    % Find maximum LL (using negative binomial)
    nsp = [];
    prednsp = [];
    for t = 1:numel(GLMP.subunit{sub}.uniqueIdx)
        nsp = cat(1,nsp,GLMP.subunit{sub}.nspikes(GLMP.subunit{sub}.uniqueIdx{t}));
        prednsp = cat(1,prednsp,repmat(GLMP.subunit{sub}.meannspikes(t),numel(GLMP.subunit{sub}.uniqueIdx{t}),1));
    end
    meanbl = repmat(mean(GLMP.subunit{sub}.blnspikes),numel(GLMP.subunit{sub}.blnspikes),1);
    prednsp = cat(1,prednsp,meanbl);
    nsp = cat(1,nsp,GLMP.subunit{sub}.blnspikes);
    mu = prednsp;
    x = nsp;
    [kappa,f] = fmincon('minimizeKappa',kappaguess,[],[],[],[],0,[],[],fitspanel.options,mu,x);
    ubLL = -f;    

    % Generating an initial guess
    paramsGuess = poppanel.params(pcidx(n),:);
    paramsGuess(2:3) = mean(paramsGuess(2:3));
    
    % Using paramsGuess as initial guess into fit
    warning off MATLAB:nearlySingularMatrix
    warning off MATLAB:illConditionedMatrix
    [sympdparams(n,:),sympdnLL(n)] = fmincon('FitNakaRushtonFunJPW',paramsGuess,...
        a,b,Aeq,Beq,lb,ub,[],options,...
        [Lcc Mcc],nsp,surftype,errortype);
    
%     if abs(ubLL-fitspanel.normLL.ubLL) > .005 | abs(lbLL - fitspanel.normLL.lbLL) > .005
%         keyboard
%     end
    sympdnLL(n) = FitNakaRushtonFunJPW(paramsGuess,[Lcc Mcc],nsp,'conicsection_xy','negativebinomial');
    sympdparams(n,:) = paramsGuess;
    
%     if -sympdnLL(n) < lbLL || -sympdnLL(n) > ubLL
%         keyboard
%     end
    
    sympdnormLL(n) = (-sympdnLL(n)-lbLL)/(ubLL-lbLL);
    
%     colormap('cool')
%     x = linspace(-max(GLMP.subunit{sub}.Lcc),max(GLMP.subunit{sub}.Mcc),50);
%     [xx,yy] = meshgrid(x,x);
%     surface = ComputeNakaRushtonJPW(sympdparams(n,:),[xx(:) yy(:)],surftype);
%     surface = reshape(surface,size(xx));
%     figure(666); cla; hold on; grid on;
%     ticks = rndofferr(linspace(min(Lcc),max(Lcc),5),2);
%     set(gca,'XTick',ticks,'YTick',ticks,'xlim',[min(xx(:)) max(xx(:))],'ylim',[min(yy(:)) max(yy(:))]);
%     p = surf(gca,xx,yy,surface);
%     set(p,'edgecolor','none')
%     alpha(.3);
%     contour3(xx,yy,surface,'linewidth',2);
%     uniqLcc = GLMP.subunit{sub}.uniqueLcc;
%     uniqMcc = GLMP.subunit{sub}.uniqueMcc;
%     meannsp = GLMP.subunit{sub}.meannspikes;
%     for i = 1:numel(uniqLcc)
%         mn = meannsp(i)/max(meannsp)*10+3;
%         h = plot3(uniqLcc(i),uniqMcc(i),meannsp(i),'ko');
%         set(h,'MarkerFaceColor','black','MarkerSize',mn,'MarkerEdgeColor','white')
%     end
%     xlabel('Lcc');
%     ylabel('Mcc');
%     zlabel('# of spikes')
%     title('Symmetric PD Fit')
        
end

unidnormLL = poppanel.twoD.unichrom.normLLs(pcidx);
bidnormLL = poppanel.twoD.bichrom.normLLs(pcidx);

disp([num2str(sum(bidnormLL - sympdnormLL < .075)) ' are well fit with a symmetric PD']);


