function GLMSD = OrganizeRawGLMSDData(rawdata)
%This code takes a stro structure of GridLMDetection data (rawdata), and
%organizes it into a GLMSD data structure, complete with original GLMS
%structure.

% Organize Dense Noise Data
GLMSD.datafile = rawdata.sum.fileName(end-13:end-4);
GLMSD.neuralFile = char(rawdata.sum.exptParams.filename);
if size(GLMSD.neuralFile,1) > size(GLMSD.neuralFile,2)
    GLMSD.neuralFile = GLMSD.neuralFile';
end
GLMSD.framerate = rawdata.sum.exptParams.framerate;
GLMSD.rf_x = rawdata.sum.exptParams.rf_x;
GLMSD.rf_y = rawdata.sum.exptParams.rf_y;
GLMSD.fpon = rawdata.trial(:,strcmp(rawdata.sum.trialFields(1,:),'fpon_t'));
GLMSD.fpacq = rawdata.trial(:,strcmp(rawdata.sum.trialFields(1,:),'fpacq_t'));
GLMSD.stimon = rawdata.trial(:,strcmp(rawdata.sum.trialFields(1,:),'stimon_t'));
GLMSD.stimoff = rawdata.trial(:,strcmp(rawdata.sum.trialFields(1,:),'stimoff_t'));
GLMSD.stimDur = GLMSD.stimoff - GLMSD.stimon;
GLMSD.DVAPerStix = rawdata.trial(:,strcmp(rawdata.sum.trialFields(1,:),'DVAPerStix'));
GLMSD.NStixGrid = rawdata.trial(:,strcmp(rawdata.sum.trialFields(1,:),'NStixGrid'));
GLMSD.Lcc = rawdata.trial(:,strcmp(rawdata.sum.trialFields(1,:),'Lcc'));
GLMSD.Mcc = rawdata.trial(:,strcmp(rawdata.sum.trialFields(1,:),'Mcc'));
GLMSD.Scc = rawdata.trial(:,strcmp(rawdata.sum.trialFields(1,:),'Scc'));
[GLMSD.theta,GLMSD.rho] = cart2pol(GLMSD.Lcc,GLMSD.Mcc);
GLMSD.RFCorrect = rawdata.trial(:,strcmp(rawdata.sum.trialFields(1,:),'RFCorrect'));
GLMSD.AnsCorrect = rawdata.trial(:,strcmp(rawdata.sum.trialFields(1,:),'CorrectAns'));
GLMSD.subunit = rawdata.sum.exptParams.subunit;

% Correcting thetas (roundoff error)
nrads = 64;
divs = -pi-pi/(nrads/2):pi/(nrads/2):pi+pi/(nrads/2);
rads = divs(2:2:end);
edges = divs(1:2:end);
[~,binIdx] = histc(GLMSD.theta,edges);
GLMSD.theta = rads(binIdx)';
if any(GLMSD.theta == -pi)
    GLMSD.theta(GLMSD.theta == -pi) = pi;
end

% Organize by RF
for n = 1:2
    str = ['RF' num2str(n)];
    L = GLMSD.RFCorrect == mod(n,2);
    GLMSD.(str).fpon = GLMSD.fpon(L);
    GLMSD.(str).fpacq = GLMSD.fpacq(L);
    GLMSD.(str).stimon = GLMSD.stimon(L);
    GLMSD.(str).stimoff = GLMSD.stimoff(L);
    GLMSD.(str).stimDur = GLMSD.stimDur(L);
    GLMSD.(str).Lcc = GLMSD.Lcc(L);
    GLMSD.(str).Mcc = GLMSD.Mcc(L);
    GLMSD.(str).Scc = GLMSD.Scc(L);
    GLMSD.(str).theta = GLMSD.theta(L);
    GLMSD.(str).rho = GLMSD.rho(L);
    GLMSD.(str).AnsCorrect = GLMSD.AnsCorrect(L);
end

% Collapse across unique stimuli
allStimCond = cat(2,GLMSD.Lcc,GLMSD.Mcc);
uniqueStim = unique(allStimCond,'rows');
GLMSD.uniqueLcc = uniqueStim(:,1);
GLMSD.uniqueMcc = uniqueStim(:,2);
[GLMSD.uniqueTheta,GLMSD.uniqueRho] = cart2pol(GLMSD.uniqueLcc,GLMSD.uniqueMcc);
GLMSD.uniqueAns = cell(size(uniqueStim,1),1);
GLMSD.pCorrect = nan(size(uniqueStim,1),1);
GLMSD.nStim = nan(size(uniqueStim,1),1);
GLMSD.nCorrect = nan(size(uniqueStim,1),1);
GLMSD.uniqueIdx = cell(size(uniqueStim,1),1);
GLMSD.uniqueAnsCorr = cell(size(uniqueStim,1),1);

% Correcting thetas (roundoff error)
[~,binIdx] = histc(GLMSD.uniqueTheta,edges);
GLMSD.uniqueTheta = rads(binIdx)';
if any(GLMSD.uniqueTheta == -pi)
    GLMSD.uniqueTheta(GLMSD.uniqueTheta == -pi) = pi;
end

for n = 1:size(uniqueStim,1)    
    L = ismember(allStimCond,uniqueStim(n,:),'rows');
    GLMSD.uniqueAnsCorr{n} = GLMSD.AnsCorrect(L)';
    GLMSD.nCorrect(n) = sum(GLMSD.AnsCorrect(L));
    GLMSD.nStim(n) = sum(L);
    GLMSD.pCorrect(n) = mean(GLMSD.AnsCorrect(L));
    GLMSD.uniqueIdx{n} = find(L);
end

% Collapse across unique stimuli by RF
for u = 1:2
    str = ['RF' num2str(u)];
    allStimCond = cat(2,GLMSD.(str).Lcc,GLMSD.(str).Mcc);
    uniqueStim = unique(allStimCond,'rows');
    GLMSD.(str).uniqueLcc = uniqueStim(:,1);
    GLMSD.(str).uniqueMcc = uniqueStim(:,2);
    [GLMSD.(str).uniqueTheta,GLMSD.(str).uniqueRho] = cart2pol(uniqueStim(:,1),uniqueStim(:,2));
    
    
    GLMSD.(str).uniqueAns = cell(size(uniqueStim,1),1);
    GLMSD.(str).pCorrect = nan(size(uniqueStim,1),1);
    GLMSD.(str).uniqueIdx = cell(size(uniqueStim,1),1);
    
    for n = 1:size(uniqueStim,1)        
        L = ismember(allStimCond,uniqueStim(n,:),'rows');
        GLMSD.(str).uniqueAns{n} = GLMSD.(str).AnsCorrect(L)';
        GLMSD.(str).pCorrect(n) = mean(GLMSD.(str).AnsCorrect(L));
        GLMSD.(str).uniqueIdx{n} = find(L);
    end
end


% Add original GLMS datafile
if ismac
    library = '/Users/jpatrickweller/Dropbox/Patrick/GLMS Data/nex files/';
elseif ispc
    library = 'C:\Users\jpweller\Dropbox\Patrick\GLMS Data\nex files\';
end
    

GLMSdata = nex2stro(char([library GLMSD.neuralFile '.nex']));
[GLMSD.GLMSdata,~] = OrganizeRawGLMSData(GLMSdata);

% % Calculate neurometric data
% GLMSD.AUC = nan(numel(GLMSD.GLMSdata.subunit{GLMSD.subunit}.uniqueLcc),1);
% for s = 1:numel(GLMSD.GLMSdata.subunit{GLMSD.subunit}.uniqueLcc)
%     idx = GLMSD.GLMSdata.subunit{GLMSD.subunit}.uniqueIdx{s};
%     spHist = GLMSD.GLMSdata.subunit{GLMSD.subunit}.fr(idx)';
%     blHist = GLMSD.GLMSdata.subunit{GLMSD.subunit}.blfr';
%     thresh = unique([spHist blHist]);
%     TPR = nan(1,length(thresh));
%     FPR = TPR;
%     for n=1:length(thresh)
%         TPR(n) = sum(spHist >= thresh(n))/length(spHist);
%         FPR(n) = sum(blHist >= thresh(n))/length(blHist);
%     end
%     GLMSD.AUC(s) = trapz(fliplr(FPR),fliplr(TPR));
% end

% Looking for rightward/leftward bias
Lans = sum(logical(GLMSD.RFCorrect) & logical(GLMSD.AnsCorrect));
Lans = Lans + sum(~GLMSD.RFCorrect & ~GLMSD.AnsCorrect);
Rans = numel(GLMSD.RFCorrect) - Lans;

figure(2); clf;
set(gcf,'position',[1150 100 250 400],'numbertitle','off',...
    'name','Right/Left Choices','toolbar','none');
bar([Lans Rans'])
xlim([0 3]);
labels{1} = 'Left';
labels{2} = 'Right';
labels = labels';
set(gca,'XTickLabelMode','manual','XTickLabel',labels,'ygrid','on')


end





