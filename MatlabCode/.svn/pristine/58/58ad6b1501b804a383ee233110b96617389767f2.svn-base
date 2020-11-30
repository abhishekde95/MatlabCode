function [GLMP,DN] = OrganizeRawGLMSData(rawdata)
%This code takes a stro structure of GridLMSubunitDN data (rawdata), and
%organizes it into two structures: Dense Noise (DN) and GridLMPlane (GLMP)
%data.
disp('Organizing Raw GLMS Data...')

DN = [];
GLMP = [];

DN = OrganizeDNData(rawdata);
GLMP = OrganizeGLMPData(rawdata);

disp('Raw GLMS Data Organized.')


end


function DN = OrganizeDNData(rawdata)

DN_L = rawdata.trial(:,strcmp(rawdata.sum.trialFields(1,:),'Epoch'))==1;

% Organize Dense Noise Data
DN.datafile = rawdata.sum.fileName(end-13:end-4);
DN.fpon = rawdata.trial(DN_L,strcmp(rawdata.sum.trialFields(1,:),'fpon_t'));
DN.fpacq = rawdata.trial(DN_L,strcmp(rawdata.sum.trialFields(1,:),'fpacq_t'));
DN.stimon = rawdata.trial(DN_L,strcmp(rawdata.sum.trialFields(1,:),'stimon_t'));
DN.stimoff = rawdata.trial(DN_L,strcmp(rawdata.sum.trialFields(1,:),'stimoff_t'));
DN.lumCC = rawdata.trial(DN_L,strcmp(rawdata.sum.trialFields(1,:),'Lcc'));
DN.colCC = rawdata.trial(DN_L,strcmp(rawdata.sum.trialFields(1,:),'Mcc'));
DN.seed = rawdata.trial(DN_L,strcmp(rawdata.sum.trialFields(1,:),'Seed'));
DN.nFrames = rawdata.trial(DN_L,strcmp(rawdata.sum.trialFields(1,:),'NFrames'));
DN.DVAPerStix = rawdata.trial(DN_L,strcmp(rawdata.sum.trialFields(1,:),'DVAPerStix'));
DN.NStixGrid = rawdata.trial(DN_L,strcmp(rawdata.sum.trialFields(1,:),'NStixGrid'));
DN.spiketimes = rawdata.ras(DN_L,strcmp(rawdata.sum.rasterCells(1,:),'sig001a'));
DN.framerate = rawdata.sum.exptParams.framerate;
DN.gammatable = rawdata.sum.exptParams.gammatable;
DN.bkgndrgb = rawdata.sum.exptParams.bkgndrgb;

% Calculate M matrix
funds = reshape(rawdata.sum.exptParams.fundamentals,numel(rawdata.sum.exptParams.fundamentals)/3,3); % have to adjust the spacing in wls
fundswls = linspace(380,780,size(funds,1))';
DN.monspds = reshape(rawdata.sum.exptParams.monspd,numel(rawdata.sum.exptParams.monspd)/3,3);
DN.wls = linspace(380,780,size(DN.monspds,1))';
DN.fundamentals = SplineSpd(fundswls,funds,DN.wls);
DN.M = DN.monspds' * DN.fundamentals;

% Set up time and spectral parameters
nframesback = 10;
msperframe = 1000/rawdata.sum.exptParams.framerate;
colordirlms = [1 1 0; -1 -1 0; 1 -1 0; -1 1 0];
idxnames = {'pLpML','mLmML','pLmML','mLpML'};
fieldnames = {'pLpM','mLmM','pLmM','mLpM'};

for f = 1:numel(fieldnames)
    DN.stats.(fieldnames{f}).STS = zeros(DN.NStixGrid(1).^2*3,nframesback);
    DN.stats.(fieldnames{f}).nspikes = 0;
end

for tr = 1:numel(DN.lumCC)
    
    randnums = getEJrandnums(DN.NStixGrid(tr)^2*DN.nFrames(tr), DN.seed(tr));
    randnums = reshape(randnums, [DN.NStixGrid(tr)^2, DN.nFrames(tr)]);
    randnums = mod(randnums, size(colordirlms,1))+1;
    randnums_orig = colordirlms(randnums,:);
    
    %Seperate +/- lum/chrom conditions (also called idxnames)
    pLpML = find(randnums_orig(:,1) > 0 & randnums_orig(:,2) > 0);
    mLmML = find(randnums_orig(:,1) < 0 & randnums_orig(:,2) < 0);
    pLmML = find(randnums_orig(:,1) > 0 & randnums_orig(:,2) < 0);
    mLpML = find(randnums_orig(:,1) < 0 & randnums_orig(:,2) > 0);
    
    for f = 1:numel(fieldnames)
        
        randnums = zeros(size(randnums_orig));
        randnums(eval(idxnames{f}),:) = randnums_orig(eval(idxnames{f}),:);
        randnums = reshape(randnums, [DN.NStixGrid(tr).^2, DN.nFrames(tr), 3]);
        randnums = permute(randnums,[1 3 2]);
        randnums = reshape(randnums, [DN.NStixGrid(tr)^2*3, DN.nFrames(tr)]);
        % Each column should be Ls followed by Ms followed by Ss for each frame.
        
        if any(DN.spiketimes{tr}) && DN.nFrames(tr) ~= 0
            frametimes = linspace(0, DN.nFrames(tr) * msperframe, DN.nFrames(tr))+(msperframe/2)';
            sptimes = (DN.spiketimes{tr}-DN.stimon(tr))*1000;
            sptimes(sptimes<0) = [];
            sptimes(sptimes>frametimes(end)) = [];
            [n,~] = hist(sptimes, frametimes);
            STCOVmex('init',{DN.NStixGrid(tr).^2,3,nframesback})
            STCOVmex(randnums(:), n);
            out = STCOVmex('return');
            DN.stats.(fieldnames{f}).STS = DN.stats.(fieldnames{f}).STS+out{1};
            DN.stats.(fieldnames{f}).nspikes = DN.stats.(fieldnames{f}).nspikes+out{3};
            DN.stats.(fieldnames{f}).STA = DN.stats.(fieldnames{f}).STS ./ DN.stats.(fieldnames{f}).nspikes;
        end
    end
end

end


function GLMP = OrganizeGLMPData(rawdata)

GLMP_L = rawdata.trial(:,strcmp(rawdata.sum.trialFields(1,:),'Epoch'))==2 |...
    rawdata.trial(:,strcmp(rawdata.sum.trialFields(1,:),'Epoch'))==3;
GLMP.datafile = rawdata.sum.fileName(end-13:end-4);
GLMP.rf_x = rawdata.sum.exptParams.rf_x;
GLMP.rf_y = rawdata.sum.exptParams.rf_y;
GLMP.fpon = rawdata.trial(GLMP_L,strcmp(rawdata.sum.trialFields(1,:),'fpon_t'));
GLMP.fpacq = rawdata.trial(GLMP_L,strcmp(rawdata.sum.trialFields(1,:),'fpacq_t'));
GLMP.stimon = rawdata.trial(GLMP_L,strcmp(rawdata.sum.trialFields(1,:),'stimon_t'));
GLMP.stimoff = rawdata.trial(GLMP_L,strcmp(rawdata.sum.trialFields(1,:),'stimoff_t'));
GLMP.Lcc = rawdata.trial(GLMP_L,strcmp(rawdata.sum.trialFields(1,:),'Lcc'));
GLMP.Mcc = rawdata.trial(GLMP_L,strcmp(rawdata.sum.trialFields(1,:),'Mcc'));
GLMP.Scc = rawdata.trial(GLMP_L,strcmp(rawdata.sum.trialFields(1,:),'Scc'));
[GLMP.theta,GLMP.rho] = cart2pol(GLMP.Lcc,GLMP.Mcc);
GLMP.gridX = rawdata.ras(GLMP_L,strcmp(rawdata.sum.rasterCells(1,:),'GridX'));
GLMP.gridY = rawdata.ras(GLMP_L,strcmp(rawdata.sum.rasterCells(1,:),'GridY'));
GLMP.stimDur = GLMP.stimoff - GLMP.stimon;
GLMP.spiketimes = rawdata.ras(GLMP_L,strcmp(rawdata.sum.rasterCells(1,:),'sig001a'));
GLMP.epoch = rawdata.trial(GLMP_L,strcmp(rawdata.sum.trialFields(1,:),'Epoch'));
GLMP.gaborL = GLMP.epoch == 3;
GLMP.gaboridx = find(GLMP.gaborL);
GLMP.flashL = GLMP.epoch == 2;
GLMP.flashidx = find(GLMP.flashL);

% Calculate M matrix
funds = reshape(rawdata.sum.exptParams.fundamentals,numel(rawdata.sum.exptParams.fundamentals)/3,3); % have to adjust the spacing in wls
fundswls = linspace(380,780,size(funds,1))';
GLMP.monspds = reshape(rawdata.sum.exptParams.monspd,numel(rawdata.sum.exptParams.monspd)/3,3);
GLMP.wls = linspace(380,780,size(GLMP.monspds,1))';
GLMP.fundamentals = SplineSpd(fundswls,funds,GLMP.wls);
GLMP.M = GLMP.monspds' * GLMP.fundamentals;

% The following conditions only apply if there is an actual GLMP dataset.
% If the datafile includes only DN data, bounce out.
if size(GLMP.fpon,1) == 0
    return
end

% Correcting thetas (roundoff error)
nrads = 64;
divs = -pi-pi/(nrads/2):pi/(nrads/2):pi+pi/(nrads/2);
rads = divs(2:2:end);
edges = divs(1:2:end);
[~,binIdx] = histc(GLMP.theta,edges);
GLMP.theta = rads(binIdx)';
if any(GLMP.theta == -pi)
    GLMP.theta(GLMP.theta == -pi) = pi;
end


%% Pull out spiking information

% Set up some variables
GLMP.nspikes = nan(size(GLMP.fpon));
GLMP.fr = nan(size(GLMP.fpon));
GLMP.normspiketimes = cell(size(GLMP.fpon));
GLMP.blnspikes = nan(size(GLMP.fpon));
GLMP.blfr = nan(size(GLMP.fpon));

% Normalized spike times
for k = 1:size(GLMP.fpon,1)
    GLMP.normspiketimes{k} = GLMP.spiketimes{k} - GLMP.stimon(k);
end

% Find baseline spikerate
for i = 1:size(GLMP.fpon,1)
    bldt = min(GLMP.stimDur); % = stimdur to make spike count equivalent between stim and baseline. 
    nsp = sum(GLMP.normspiketimes{i} < 0 & GLMP.normspiketimes{i} > -bldt);
    GLMP.blnspikes(i) = nsp;
    GLMP.blfr(i) = nsp./bldt;
end

% Countingthspikes after first significant deviation from baseline
binsize = .001; % in seconds
psthbins = -.2:binsize:.5;
blbins = -.2:binsize:0;
PSTH = histc(cat(1,GLMP.normspiketimes{GLMP.flashL}),psthbins);
PSTH = PSTH./sum(GLMP.flashL)./binsize;

% Bin and smooth baseline 
gaussSize = .01;
gaussSize = ceil(gaussSize/binsize);
gaussfilter = gausswin(gaussSize,3);
gaussfilter = gaussfilter./sum(gaussfilter); %normalizing
blpsth = histc(cat(1,GLMP.normspiketimes{:}),blbins);
blpsth = blpsth./numel(GLMP.normspiketimes)./binsize;
smoothblPSTH = conv(blpsth,gaussfilter,'same');
smoothPSTH = conv(PSTH,gaussfilter,'same');

% Define threshold for signal vs baseline noise
threshblfr = max(smoothblPSTH) * 1.2;
temp = find(smoothPSTH' > threshblfr & psthbins > .025 & psthbins < .2);
if isempty(temp) % if modulation never increases beyond BL threshold,
    [~,temp] = min(abs(psthbins-0.08)); % take fixed window of 80-280ms
    %[~,temp] = max(smoothPSTH)- 100;
end
GLMP.countingwin = [psthbins(temp(1)) psthbins(temp(1))+mean(GLMP.stimDur(GLMP.flashL))];
GLMP.blfrthresh = threshblfr;

% Count flash spikes 
for k = 1:sum(GLMP.flashL)
    idx = GLMP.flashidx(k);
    GLMP.nspikes(idx) = sum(GLMP.normspiketimes{idx} >= GLMP.countingwin(1) & GLMP.normspiketimes{idx} <= GLMP.countingwin(2));
end
GLMP.fr(GLMP.flashL) = GLMP.nspikes(GLMP.flashL)./GLMP.stimDur(GLMP.flashL);

% Count gabor spikes 
for k = 1:sum(GLMP.gaborL)
    idx = GLMP.gaboridx(k);
    GLMP.nspikes(idx) = sum(GLMP.normspiketimes{idx} > 0 & GLMP.normspiketimes{idx} < GLMP.stimDur(idx)+.2);
end
GLMP.fr(GLMP.gaborL) = GLMP.nspikes(GLMP.gaborL)./GLMP.stimDur(GLMP.gaborL);


%% Organize Subunits
% Alter gridX and gridY for gabors to treat them as a distinct subunit
%GLMP.gaborL = GLMP.epoch == 3;
[GLMP.gridX{GLMP.gaborL}] = deal(100);
[GLMP.gridY{GLMP.gaborL}] = deal(100);
temp = ones(size(GLMP.fpon));

% Matching subnits by gridX and gridY (Gabors = 100, see above)
for s = 1:3
    matchIdx = zeros(size(GLMP.fpon));
    for n = 1:size(GLMP.fpon,1)
        if ~isempty(find(temp,1))
            if numel(GLMP.gridX{find(temp,1)}) == numel(GLMP.gridX{n})
                matchIdx(n) = all(GLMP.gridX{find(temp,1)} == GLMP.gridX{n})...
                    && all(GLMP.gridY{find(temp,1)} == GLMP.gridY{n});
            else
                matchIdx(n) = 0;
            end
        end
    end
    matchIdx = logical(matchIdx);
    
    % Now assign quantities from the trial index into subunits
    if sum(matchIdx) ~= 0
        GLMP.subunit{s}.datafile = GLMP.datafile;
        GLMP.subunit{s}.GLMPIdx = find(matchIdx);
        GLMP.subunit{s}.fpon = GLMP.fpon(matchIdx);
        GLMP.subunit{s}.fpacq = GLMP.fpacq(matchIdx);
        GLMP.subunit{s}.stimon = GLMP.stimon(matchIdx);
        GLMP.subunit{s}.stimoff = GLMP.stimoff(matchIdx);
        GLMP.subunit{s}.stimDur = GLMP.stimDur(matchIdx);
        GLMP.subunit{s}.Lcc = GLMP.Lcc(matchIdx);
        GLMP.subunit{s}.Mcc = GLMP.Mcc(matchIdx);
        GLMP.subunit{s}.Scc = GLMP.Scc(matchIdx);
        GLMP.subunit{s}.theta = GLMP.theta(matchIdx);
        GLMP.subunit{s}.rho = GLMP.rho(matchIdx);
        GLMP.subunit{s}.gridX = cat(1,GLMP.gridX(matchIdx));
        GLMP.subunit{s}.gridY = cat(1,GLMP.gridY(matchIdx));
        GLMP.subunit{s}.spiketimes = GLMP.spiketimes(matchIdx);
        GLMP.subunit{s}.normspiketimes = GLMP.normspiketimes(matchIdx);
        GLMP.subunit{s}.nspikes = GLMP.nspikes(matchIdx);
        GLMP.subunit{s}.fr = GLMP.fr(matchIdx);
        GLMP.subunit{s}.blnspikes = GLMP.blnspikes(matchIdx);
        GLMP.subunit{s}.blfr = GLMP.blfr(matchIdx);
        
        % Collapsing across like stimuli and computing mean spike rates
        x = GLMP.subunit{s}.Lcc;
        y = GLMP.subunit{s}.Mcc;
        uniqueLM = unique([x y],'rows');
        GLMP.subunit{s}.meannspikes = nan(size(uniqueLM,1),1);
        GLMP.subunit{s}.varnspikes = nan(size(uniqueLM,1),1);
        GLMP.subunit{s}.stderr = nan(size(uniqueLM,1),1);
        GLMP.subunit{s}.meanfr = nan(size(uniqueLM,1),1);
        GLMP.subunit{s}.spiketimes_cat = cell(size(uniqueLM,1),1);
        GLMP.subunit{s}.uniqueIdx = cell(size(uniqueLM,1),1);
        GLMP.subunit{s}.uniqueLcc = uniqueLM(:,1);
        GLMP.subunit{s}.uniqueMcc = uniqueLM(:,2);
        for t = 1:size(uniqueLM,1)
            L = (x == uniqueLM(t,1)) & (y == uniqueLM(t,2));
            trlidx = find(L);
            GLMP.subunit{s}.uniqueIdx{t} = trlidx;
            nspikes = GLMP.subunit{s}.nspikes(trlidx);
            for j = 1:numel(trlidx)
                GLMP.subunit{s}.spiketimes_col{t,j} = GLMP.subunit{s}.normspiketimes{trlidx(j)};
            end
            GLMP.subunit{s}.spiketimes_cat{t} = cat(1,GLMP.subunit{s}.spiketimes_col{t,:});
            GLMP.subunit{s}.meannspikes(t) = mean(nspikes);
            GLMP.subunit{s}.varnspikes(t) = var(nspikes);
            GLMP.subunit{s}.stderr(t) = stderr(nspikes);
            GLMP.subunit{s}.meanfr(t) = mean(nspikes./GLMP.subunit{s}.stimDur(trlidx));
        end
        temp(matchIdx) = 0;
        
        %Unique thetas and rhos
        [GLMP.subunit{s}.uniquetheta,GLMP.subunit{s}.uniquerho] = cart2pol(GLMP.subunit{s}.uniqueLcc,GLMP.subunit{s}.uniqueMcc);
        
        % Correcting thetas (roundoff error)
        nrads = 64;
        divs = -pi-pi/(nrads/2):pi/(nrads/2):pi+pi/(nrads/2);
        rads = divs(2:2:end);
        edges = divs(1:2:end);
        [~,binIdx] = histc(GLMP.subunit{s}.uniquetheta,edges);
        GLMP.subunit{s}.uniquetheta = rads(binIdx)';
        if any(GLMP.subunit{s}.uniquetheta == -pi)
            GLMP.subunit{s}.uniquetheta(GLMP.subunit{s}.uniquetheta == -pi) = pi;
        end
        
        % Calculating neurometric data
        GLMP.subunit{s}.AUC = nan(numel(GLMP.subunit{s}.uniqueLcc),1);
        for t = 1:numel(GLMP.subunit{s}.AUC)
            idx = GLMP.subunit{s}.uniqueIdx{t};
            spHist = GLMP.subunit{s}.fr(idx)';
            blHist = GLMP.subunit{s}.blfr';
            thresh = unique([spHist blHist]);
            TPR = nan(1,length(thresh));
            FPR = TPR;
            for n=1:length(thresh)
                TPR(n) = sum(spHist >= thresh(n))/length(spHist);
                FPR(n) = sum(blHist >= thresh(n))/length(blHist);
            end
            GLMP.subunit{s}.AUC(t) = trapz(fliplr(FPR),fliplr(TPR));
        end
        
    else
        GLMP.subunit{s} = [];
    end
end

% Arrange subunits in correct order (sub1, sub2, subunit1+subunit2/Gabor)
if GLMP.subunit{1}.gridX{1} == 100
    tmpsub3 = GLMP.subunit{1};
    GLMP.subunit{1} = GLMP.subunit{3};
    GLMP.subunit{3} = tmpsub3;
elseif ~isempty(GLMP.subunit{2}) && any(GLMP.subunit{2}.gridX{1} == 100)
    tmpsub3 = GLMP.subunit{2};
    GLMP.subunit{2} = GLMP.subunit{3};
    GLMP.subunit{3} = tmpsub3;
elseif ~isempty(GLMP.subunit{3}) && all(GLMP.subunit{3}.gridX{1} ~= 100)
    [~,whichismax] = max([numel(GLMP.subunit{1}.gridX{1}),numel(GLMP.subunit{2}.gridX{1}),numel(GLMP.subunit{3}.gridX{1})]);
    tmpsub3 = GLMP.subunit{whichismax};
    GLMP.subunit{whichismax} = GLMP.subunit{3};
    GLMP.subunit{3} = tmpsub3;
end
if isempty(GLMP.subunit{1}) && ~isempty(GLMP.subunit{2})
    GLMP.subunit{1} = GLMP.subunit{2};
    GLMP.subunit{2} = [];
end


end

