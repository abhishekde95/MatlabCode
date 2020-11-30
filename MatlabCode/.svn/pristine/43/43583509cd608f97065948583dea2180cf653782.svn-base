% CellScreen.m ( this is a nonsense change)
%
% Contents:
% ---------
% Section 1) Going through the directory tree (down one layer only) and looking at
% all the grating data files that meet a specific set of criteria.  They have to
% have color tuning data, at least three trials per condition.  I'm
% wondering 
% if I'll be able to see relationships between color responses and spatial frequency 
% tuning (due to, say, chromatic aberration).
%
% Section 1.1) Finding all white noise file that have both gun and cone
% noise (with reasonable numbers of trials and reasonable spike rates for
% each).
%
% Section 2) Finding DTspot data files with QUEST contrast updating and consistent
% stimulus parameters (for analysis of fixational saccades).
%
% Section 3) Finding a set of cells studied with white noise whose STAs that are
% significantly L-M opponent. Also, single- and double-opponent cells.
%
% Section 4) Selecting NeuroThresh files that meet a few basic criteria.  Code adapted 
% from code written by Zack Lindbloom-Brown
%
% Section 5) Selecting NeuroThresh files that are paired (same cell) and
% finding pairs that have differ (mostly) in threshold.  Some difference in
% latency might be tolerated.  Also initial color directions.
%
% Section 6) Finding NeuroThresh files that were tested with "alternative"
% initial color directions.  Also pulling other header information.
%
% Section 7) Find NeuroThresh planar cells with L+M cone weights.
%
% ------------------------------
% DTSPOT
% Section 8) Finding DTspot files (for the DTMacPig project) that have
% aberrant thresholds
%
% Section 8.1) As above but for DTscotopic files.
%
% Section 9) Finding DTspot files that were done at 25 (or 15 Hz).
% Also getting target positions.
%
% ------------------------------
% SMurray size constancy
% Section 10) Going down one level of the directory tree and finding all
% the gray and corridor files that were collected on the same days.
%
% Section 10.1) Getting Some key numbers from size illusion experiments.
%
% ------------------------------
% DTNT
% Section 11) How many trials/threshold measurement?
%
% ------------------------------
% LMTF
% Section 12) How long was the stimulus on the screen? Until 4/16/16 LMTF.d
% requested the stimulus of the slave by passing nrampframes and
% nplatframes. These numbers were calculated on the basis of a 75 Hz
% refresh, but when I noticed this issue, the slave monitor was refreshing
% at 85 Hz. (13 frames ramping up, 24 plateau frames, 13 frames ramping down)
%
% Section 13) Geting RF locations for a bunch of LMTF nex files. (Also
% saccade latencies. Can I use this as a proxy for good behavior? Also
% confirming LMTF (as opposed to IsoSamp).)
%
% ------------------------------
% IsoSamp
% Section 14) Getting a list of all the RF locations and chucking any that
% are outside of some window (2-14°?)
%
% Section 15) Getting fixation window fixation point sizes. 
%
% Section 16) Getting the size and number of white noise stixels.
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
% SECTION 1
% Going through the directory tree (down one layer only) and looking at all
% the grating data files that meet a specific set of criteria.  They have to
% have color tuning data and at least three trials per condition.
%%%%%%%%%%%%%%%%%%%%%%%%%%
PARADIGMID = 150;
datapath = 'N:\NexFiles\Greg';
eval(['cd ',datapath]);
startdate = '110210';  % month, day, year (last two digits)
startdatenum = datenum(str2double(['20',startdate([5 6])]), str2double(startdate([1 2])), str2double(startdate([3 4])));
% First making a gigantic dirstruct
dirstruct = dir;
for i = 1:size(dirstruct,1)
    if (dirstruct(i).isdir && ~strncmp(dirstruct(i).name, '.',1))
        cd(dirstruct(i).name);
        dirstruct = [dirstruct; dir];
        cd('../')
    end
end
counter = 1;
filenames = {};
spikenumbers = [];
for i = 1:size(dirstruct,1)
    if (length(dirstruct(i).name) < 4)
        continue;
    end
    [num2str(i),' ',dirstruct(i).name]
    if (strcmp(dirstruct(i).name(end-3:end),'.nex'))
        if (length(dirstruct(i).name) >= 13)  % minimum # of characters for file named "<prefix><date><number>.nex"
            % Check the date
            filedate = dirstruct(i).name(end-12:end-7);
            filedatenum = datenum(str2double(['20',filedate([5 6])]), str2double(filedate([1 2])), str2double(filedate([3 4])));
            if (filedatenum <= startdatenum)
                continue
            end
        end
        try
            stro = nex2stro(findfile(dirstruct(i).name));
        catch
            disp(['skipping a file: ',dirstruct(i).name]);
            continue;
        end
        if ~isfield(stro,'sum')
             disp(['skipping a file: ',dirstruct(i).name]);
            continue;
        end
        if ~(stro.sum.paradigmID == PARADIGMID)
             disp(['skipping a file: ',dirstruct(i).name]);
            continue;
        end
        ntrials = size(stro.trial,1);
        orients = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'orient'));
        sfs = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'sf'));
        diams = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'diam'));
        Lcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'lcont'));
        Mcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'mcont'));
        Scc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'scont'));
        protocols = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'protocol'));
        if (isempty(protocols) || ~any(protocols == 4))
            continue;
        end
        
        % This part is having problems.
        % Need to include a check on minmum number of 
        % trials per cond, but this isn't doing it.
        %trialspecs = [orients sfs diams Lcc Mcc Scc];
        %trialtypes = unique(trialspecs,'rows');
        %n = [];
        %for j = 1:length(trialtypes)
        %    n(j) = sum(all(repmat(trialtypes(j,:),ntrials,1) == trialspecs,2));
        %end
        %if (any(n < 3))
        %     disp(['skipping a file: ',dirstruct(i).name]);
        %     keyboard
        %    continue;
        %end
        %%%%%%%%%%%%%%%
        for spikename = {'sig001a','sig001b'}
            if (~strcmp(spikename,stro.sum.rasterCells))
                continue;
            end
            filenames{counter} = dirstruct(i).name;
            spikenumbers(counter) = find(strcmp(spikename,{'sig001a','sig001b'}));
            counter = counter+1;
        end
    end
end
strippedsuffix = cellfun(@(x) x(1:end-4), filenames,'UniformOutput',0);
disp([char(strippedsuffix'),repmat('   ',counter-1,1),num2str(spikenumbers')])


%%
% Section 1.1
% Going down a directory tree looking for white noise data sets in which
% both gun and cone noise were interleaved. Also keeping track of numbers
% of trials and spike rates.

PARADIGMID = 100;
datapath = 'N:\NexFiles\Greg';
eval(['cd (''',datapath,''')']);
% First making a gigantic dirstruct
dirstruct = dir;
subdirs = {};
for i = 1:size(dirstruct,1) % Getting the initial list of directories
    if (dirstruct(i).isdir && ~strncmp(dirstruct(i).name, '.',1))
        subdirs{length(subdirs)+1} = [datapath, filesep,dirstruct(i).name];
    end
end
alldone = 0;
startsubdir = 1;
while ~alldone
    nsubdirs = length(subdirs);
    for i = startsubdir:length(subdirs)
        i
        eval(['cd (''',subdirs{i},''')']);
        dirstruct = dir;
        for j = 1:size(dirstruct,1)
            if (dirstruct(j).isdir && ~strncmp(dirstruct(j).name, '.',1));
                datapath = eval('pwd');
                subdirs{length(subdirs)+1} = [datapath, filesep,dirstruct(j).name];
            end
        end
    end
    if length(subdirs) == nsubdirs
        alldone = 1;
    else
        startsubdir = nsubdirs;
    end
end
% OK, now we have a complete list of subdirectories that we can go through,
% one by one, looking for whitenoise files.
data = [];
for i = 1:length(subdirs)
    i
    cd(subdirs{i});
    dirstruct = dir;
    for j = 1:length(dirstruct)
        if (length(dirstruct(j).name) < 4)
            continue;
        end
        if (dirstruct(j).name(end-3:end) ~= '.nex')
            continue;
        end
        try
            if (getparadigmID(dirstruct(j).name)~=PARADIGMID)
                continue;
            end
            stro = nex2stro(dirstruct(j).name);
        catch
            disp(['Problem with ',dirstruct(j).name]);
            continue
        end
        if (~isfield(stro,'sum'))  % If there's a error
            continue
        end
        noisetypeidx = strcmp(stro.sum.trialFields(1,:),'noise_type');
        noisetype = stro.trial(:,noisetypeidx);
        if ~any(noisetype == 1) | ~any(noisetype == 2) % both gun and cone noise in file?
            continue
        end
        
        % Stripping out any subunit stuff
        ntrials = size(stro.trial,1);
        subunitmaskidx = strcmp(stro.sum.rasterCells,'subunit_mask');
        synthimageidx = strcmp(stro.sum.rasterCells,'synth_image');
        subunitmask = stro.ras(:,subunitmaskidx);
        synthimage = stro.ras(:,synthimageidx);
        if ~isempty(subunitmask) & ~isnan(subunitmask{1})
            L = false(ntrials,1);
            for k = 1:ntrials
                if all(subunitmask{k} == 0) & isnan(synthimage{k})
                    L(k) = 1;
                end
            end
            stro.trial = stro.trial(L,:);
            stro.ras = stro.ras(L,:);
            ntrials = sum(L);
        end
        Lguntrials = stro.trial(:,noisetypeidx) == 1;
        Lconetrials = stro.trial(:,noisetypeidx) == 2;
        stimonidx = strcmp(stro.sum.trialFields(1,:),'stim_on');
        stimon_t = stro.trial(:,stimonidx);
        stimoffidx = strcmp(stro.sum.trialFields(1,:),'all_off');
        stimoff_t = stro.trial(:,stimoffidx);
        spikeidx = strcmp(stro.sum.rasterCells,'sig001a');
        firingrate = nan*ones(ntrials,1);
        for k = 1:ntrials
            sp = stro.ras{k,spikeidx};
            nsp = sum(sp>stimon_t(k) & sp<stimoff_t(k));
            firingrate(k) = nsp./(stimoff_t(k)-stimon_t(k));
        end
        nguntrials = sum(Lguntrials);
        nconetrials = sum(Lconetrials);
        gunfiringrate = mean(firingrate(Lguntrials));
        conefiringrate = mean(firingrate(Lconetrials));
        data(length(data)+1).name = dirstruct(j).name;
        data(length(data)).nguntrials = nguntrials;
        data(length(data)).nconetrials = nconetrials;
        data(length(data)).gunfiringrate = gunfiringrate;
        data(length(data)).conefiringrate = conefiringrate;
    end
end

% Getting filenames
for i = 1:length(data)
    disp([data(i).name,' ', num2str(data(i).nguntrials),' ', num2str(data(i).gunfiringrate),' ', num2str(data(i).nconetrials),' ', num2str(data(i).conefiringrate)]);
end
figure; axes; hold on;
plot([data.nguntrials],'m-');
plot([data.nconetrials],'k-');

%% 
% Section 2
% Finding DTspot data files that meet a variety of uniformity criteria
% (QUEST, RF location, color directions, etc.)

filenames = fnamesFromTxt;
data = zeros(size(filenames,1), 11);
epochdurs = zeros(size(filenames,1),3,2);
for i = 1:size(filenames,1)
    stro = nex2stro(findfile(filenames(i,:)));
    spatialPeriods = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'gabor_lambda'));
    
    if isempty(stro.sum.analog.sigid) || isempty(stro.trial)
        data(i,1) = 1;
    end
    if (stro.sum.exptParams.flash_length ~= 666)
        data(i,2) = 1;
    end
    if (stro.sum.exptParams.flash_size ~= 2)
        data(i,3) = 1;
    end
    if (stro.sum.exptParams.frame_offset ~= 5)
        data(i,4) = 1;
    end
    if (stro.sum.exptParams.rf_x ~= -50)
        data(i,5) = 1;
    end
    if (stro.sum.exptParams.rf_y ~= -35)
        data(i,6) = 1;
    end
    if (~all(stro.sum.exptParams.RF_colors([1:6]) == [7 7 7 3 -3 0]'))
        data(i,7) = 1;
    end
    if (~all(stro.trial(:,3) == 1) & ~all(isnan(stro.trial(:,3))))  % Frames shown
        data(i,8) = 1;
    end
    if (~all(stro.trial(:,20) == 8))  % sigma
        data(i,9) = 1;
    end
    if (stro.sum.exptParams.expt_meth ~= 3)  % QUEST
        data(i,10) = 1;
    end
    if (~all(ismember([206 58 16],spatialPeriods)))  % SF
        data(i,11) = 1;
    end
    frame_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'frame_on'));
    stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_on'));
    stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'flash_off'));
    targ_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targ_on'));
    epochdurs(i,1,[1 2]) = [min(stimon_t-frame_t) max(stimon_t-frame_t)];
    epochdurs(i,2,[1 2]) = [min(stimoff_t-stimon_t) max(stimoff_t-stimon_t)];
    epochdurs(i,3,[1 2]) = [min(targ_t-stimoff_t) max(targ_t-stimoff_t)];
    data
end

badidxs = find(sum(data,2))';
for i = badidxs
   stro = nex2stro(findfile(filenames(i,:)));
   keyboard
end

%%
% Section 3
% Finding cells with L-M opponent STAs from white noise.
% This is going to be problematic for double opponent neurons
% because we're averaging across space.
% Also a screen for single-opponent cells.

[filenames, spikenums] = fnamesFromTxt2;
maxT = 9;
data = [];
alpha = 0.01; % for cone weight analysis
for i = 1:size(filenames,1)
    stro = {};
    for f = 1:size(filenames{i},1)
        tmpstro = nex2stro(findfile(char(filenames{i}(f))));
        if (isempty(stro))
            stro = tmpstro;
        else
            stro = strocat(stro, tmpstro);
        end
    end
    noisetypeidx = find(strcmp(stro.sum.trialFields(1,:),'noise_type'));
    sigmaidxs = strmatch('sigma',stro.sum.trialFields(1,:));
    muidxs = strmatch('mu',stro.sum.trialFields(1,:));
    spikeidx = find(strcmp(stro.sum.rasterCells(1,:),'sig001a'));
    nstixperside = stro.sum.exptParams.nstixperside;
    
    % Reconstructing the M matrix and gamma table
    fundamentals = stro.sum.exptParams.fundamentals;
    fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
    mon_spd = stro.sum.exptParams.mon_spd;
    mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
    mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
    M = fundamentals'*mon_spd;
    gammaTable = stro.sum.exptParams.gamma_table;
    gammaTable = reshape(gammaTable, length(gammaTable)/3, 3);
    
    % Getting the background rgb/lms
    ridx = find(strcmp(stro.sum.trialFields(1,:),'bkgnd_r'));
    gidx = find(strcmp(stro.sum.trialFields(1,:),'bkgnd_g'));
    bidx = find(strcmp(stro.sum.trialFields(1,:),'bkgnd_b'));
    bkgndRGB = [mode(stro.trial(:,ridx)), mode(stro.trial(:,gidx)), mode(stro.trial(:,bidx))];
    bkgndrgb = [gammaTable(bkgndRGB(1)+1,1); gammaTable(bkgndRGB(2)+1,2); gammaTable(bkgndRGB(3)+1,3)];
    bkgndlms = M*bkgndrgb;
    
    L = stro.trial(:,noisetypeidx) == 1;  % Gun noise only.
    muvect = unique(stro.trial(L,muidxs),'rows')/1000;
    sigmavect = unique(stro.trial(L,sigmaidxs),'rows')/1000;
    sigmavect = sigmavect(~all(sigmavect == 0,2),:);
    gausslims = [stro.sum.exptParams.gauss_locut stro.sum.exptParams.gauss_hicut]/1000;
    mumat = repmat(muvect,[nstixperside^2,1]);
    sigmamat = repmat(sigmavect,[nstixperside^2,1]);
    
    stro.ras(~L,:) = [];
    stro.trial(~L,:) = [];
    out = getWhtnsStats(stro,maxT,'STCOVmex', {nstixperside^2, 3, maxT}, ['sig001',abs('a')+spikenums(i)-1]);
    
    STSs = out{1};
    STCs = out{2};
    nspikes = out{3};
    peakframe = find((sum(STSs.^2)== max(sum(STSs.^2))));
    if peakframe == 1 | peakframe == maxT
        error('peakframe is at an extreme. Cannot take a frame on either side.');
    end
    STA = STSs./nspikes;
    
     whichpix = getRFfromSTA(STA(:,peakframe+[-1 0 1]),2,.9,sigmavect(1).^2/nspikes);
%     if sum(whichpix(:)) == 1
%         single_opponent = nan;
%     else
    %threeframeSTA = reshape(sum(STA(:,peakframe+[-1 0 1]),2),[nstixperside nstixperside 3]);
    %threeframeSTA = reshape(threeframeSTA,[nstixperside^2 3]);
    %threeframeSTA = threeframeSTA(:,:);
    %[u,s,v] = svd(threeframeSTA);
    [u,s,v] = svd(reshape(STA(:,peakframe),[nstixperside^2,  3]));
    ft = fft2(reshape(u(:,1),10,10)); % 
    subplot(2,1,1); imagesc(reshape(u(:,1),10,10))
    subplot(2,1,2); imagesc(abs(ft));
    if max(abs(ft(:)))== abs(ft(1))
        single_opponent = 1;
    else
        single_opponent = 0;
    end
    
    STA = reshape(STA(:,peakframe),[nstixperside^2 3]);
    %STA = reshape(STSs.^2/nspikes,[100 3 maxT]);
    STC = reshape(STCs(:,peakframe),[300 300])./nspikes;  % Remember to divide by 'n' to get STC of STA
    %zscores = sqrt(nspikes)*(STA-mumat)./sigmamat;
    %whichpix = logical(sum(zscores.^2,2) > chi2inv(1-alpha,3));
        
%     % Checking for spatial opponency
%     % Q: how many STA stixels in the half of color space *opposite* the
%     % hottest stixel?
%     hottestpix = sum(STA.^2,2) == max(sum(STA.^2,2));
%     gunweights_hottestpix = STA(hottestpix,:);
%     dotprods = STA(find(whichpix),:)*gunweights_hottestpix';
%     if length(dotprods) > 1 & sum(dotprods < 0)./length(dotprods) <= 0.1 % Less than 10% of stixels on the "wrong" side
%         single_opponent = 1;
%     else
%         single_opponent = 0;
%     end
%     
    
    gunweights = sum(STA(find(whichpix),:),1);
    coneweights = inv(M')*gunweights';
    %   DKLmat = mkbasis([1 -1 0; 0 0 1; 1 1 0]')';
    %   DKLweights = DKLmat*coneweights;
    
    STCrgb = zeros(3);
    for j = 1:3
        for k = 1:3
            STCrgb(j,k) = sum(sum(STC(100*(j-1)+find(whichpix),100*(k-1)+find(whichpix))));
        end
    end
    STCrgb = STCrgb./nspikes;  % Dividing by 'n' because we're doing STA not STS
    STClms = inv(M')*STCrgb*inv(M);
    % STCdkl = DKLmat*STClms*DKLmat';
    crit = finv(1-alpha,3,nspikes-3);
    crit = crit*((nspikes-1)*3)/(nspikes-3);
    bounds = FindEllipsoidExtrema(STClms,crit);
    data = [data; coneweights' bounds abs(coneweights') > bounds sum(whichpix(:)) single_opponent] 
    % 1-3: cone weights
    % 4-6: critical values on cone weights
    % 7-9: significance for each cone weight
    % 10: number of stixels deemed to be in RF
    % 11: single opponent
end

accepted = find(sign(data(:,1)) ~= sign(data(:,2)) &  data(:,7) == 1 & data(:,8) == 1);
for i = 1:length(accepted)
    str = [];
    for j = 1:length(filenames{accepted(i)})
        str = [str, sprintf('%s ',char(filenames{accepted(i)}(j)))];
        
    end
    str = [str, sprintf('\t%d',spikenums(accepted(i)))];
    disp(str);
end
%%
% Section 4
% Looking through files for reasonable NeuroThresh data.
% Option to look through the entire directory or just a text file.
whichfiles = 'TEXTFILE'; % 'ALL or 'TEXTFILE';
fnamesthatpassed = {};
fnamesthatfailed = {};

if strcmp(whichfiles, 'ALL') 
    datapath = 'N:\NexFiles';
    dirstruct = subdir([datapath filesep '*.nex']);
    nfiles = size(dirstruct,1);
elseif strcmp(whichfiles, 'TEXTFILE')
   [fnames, spikeIdx] = fnamesFromTxt2();
   fnames = [fnames{:}]';
   nfiles = size(fnames,1);
end

for i = 1:nfiles
    i
    if strcmp(whichfiles, 'ALL')
        try
            filename = dirstruct(i).name;
            stro = nex2stro(filename);
        catch exception
            wh = findall(0, 'tag', 'TMWWaitbar');
            if ~isempty(wh), delete(wh); end
            continue
        end
    elseif strcmp(whichfiles, 'TEXTFILE')
        filename = char(fnames{i});
        stro = nex2stro(findfile(filename));
    end
    if (stro.sum.paradigmID ~= 103)
        continue;
    end
    
    out1 = NTpreprocess(stro,0,inf,inf); %prelim filter
    if ~isempty(out1)
        if length(unique(out1(:,1))) < 20 % at least n color dirs tested
            disp('Fewer than 20 color dirs');
            filename
            fnamesthatfailed = [fnamesthatfailed; filename];
            continue
        end
    else
        fnamesthatfailed = [fnamesthatfailed; filename];
        continue;
    end
    
    out2 = NTpreprocess(stro,0,inf,0); %stringent filter
    if ~isempty(out2)
        if length(out2(:,end)) - sum(out2(:,end) == 1) < 3 % at least n nOOG points
            fnamesthatfailed = [fnamesthatfailed; filename];
            continue
        end
    else
        fnamesthatfailed = [fnamesthatfailed; filename];
        continue;
    end
    
    spikeidx = strcmp(stro.sum.rasterCells(1,:), getSpikenum(stro, 'first'));
    coloridxs = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'coloridx'));
    levels = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'levelidx'));
    stepsize = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stepsize'));
    
    lmsidxs = [find(strcmp(stro.sum.trialFields(1,:),'lcont')) ...
        find(strcmp(stro.sum.trialFields(1,:),'mcont')) ...
        find(strcmp(stro.sum.trialFields(1,:),'scont'))];
    
    fpacq_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fp_acq'));
    stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_on'));
    stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_off'));
    
    lms = stro.trial(:,lmsidxs);

    try
        Lbaseline = ~levels & all(lms == zeros(size(lms)),2);
        Lprobe = ~levels & ~all(lms == zeros(size(lms)),2);
        Lreplay = (levels == max(levels)) & ~stepsize;
    catch exception
        keyboard;
        continue;
    end
    % In what proportion of trials does the baseline activity exceed threshold?
    L = ~isnan(fpacq_t) & ~Lreplay;
    prestimfr = [];
    for j = find(L)'
        spiketimes = stro.ras{j,spikeidx};
        nspikes = sum(spiketimes > fpacq_t(j) & spiketimes < stimon_t(j));
        prestimfr = [prestimfr; nspikes./(stimon_t(j)-fpacq_t(j))];
    end
    abovethresh = prestimfr > stro.sum.exptParams.threshold;
    if sum(abovethresh)/numel(prestimfr) > 0.05
        disp('5% threshold violations');
        filename
        fnamesthatfailed = [fnamesthatfailed; filename];
        continue
    end
    
    
    % t-tests on the probe trials.  This may have been too stringent a
    % criterion (and unfairly penalizes cells with few trials!)
  %  spikerates = zeros(size(stro.trial,1),1);
  %  for j = 1:size(stro.trial,1)
  %      spiketimes = stro.ras{j,spikeidx};
  %      nspikes = sum(spiketimes > stimon_t(j)+stro.sum.exptParams.latency/1000 & spiketimes < stimoff_t(j));
  %      spikerates(j) = nspikes./(stimoff_t(j)-stimon_t(j)-stro.sum.exptParams.latency/1000);
  %  end
    
  %  probefrs = spikerates(Lprobe);
  %  baselinefrs = spikerates(Lbaseline);
  %  if length(probefrs) > 1
  %      halflen = floor(length(probefrs)/2);
  %      [~,pval] = ttest2(probefrs(1:halflen), probefrs(halflen+1:end));
      %  if pval < 0.05
      %      continue
      %  end
  %  else
  %      continue
  %  end
    if (isnan(stro.sum.exptParams.monopolar))
        % Do nothing
    elseif (stro.sum.exptParams.monopolar == 0)
        % Do nothing
    else % If this is a monopolar data file
        fnamesthatfailed = [fnamesthatfailed; filename];
        continue;
    end
  
    fnamesthatpassed = [fnamesthatpassed; filename];
end

%%
% Section 5
% Finding files that differ mostly in threshold, or initial color
% directions, without too much variation in latency.  Also looking at the
% number of color directions test and the total number of data points.
[fnames, spikeIdx] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\NT\NTmultiples.txt');
%[fnames, spikeIdx] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\NT\altinitcoldirs.txt');

maxnfiles = 0;
for i = 1:size(fnames,1)
    maxnfiles = max(maxnfiles, length(fnames{i}));
end

threshs = nan*ones(length(fnames), maxnfiles);
lats = nan*ones(length(fnames), maxnfiles);
ncolordirstested = nan*ones(length(fnames), maxnfiles);
nnoogs = nan*ones(length(fnames), maxnfiles);
for cellcounter = 1:size(fnames,1)
    NT = {}; 
    for i = 1:size(fnames{cellcounter},2)
        filename = findfile(char(fnames{cellcounter}(i)));
        NT{i} = nex2stro(filename);
    end
    for i = 1:length(NT)  % basic checks
        out1 = NTpreprocess(NT{i},0,inf,inf); %prelim filter
        out2 = NTpreprocess(NT{i},0,inf,0); %stringent filter
        ncolordirstested(cellcounter,i) = length(unique(out1(:,1))); % # color dirs tested
        nnoogs(cellcounter,i) = length(out2(:,end)) - sum(out2(:,end) == 1); % # nOOG points
    end
    
    f = fieldnames(NT{1}.sum.exptParams);
    mismatchedfields = {};
    for i = 1:length(f)
        tmp = [];
        for j = 1:length(NT)
            tmp(j,:) = getfield(NT{j}.sum.exptParams,f{i});
        end
        if (strcmp(f(i),'threshold'))
            threshs(cellcounter, 1:length(NT)) = tmp';
        end
        if (strcmp(f(i),'latency'))
             lats(cellcounter, 1:length(NT)) = tmp';
        end

        if(size(unique(tmp,'rows'),1) ~= 1)
            if (~all(isnan(tmp)))
                % We've found a mismatch
                f(i)
                tmp
            end
        end
    end  
end

threshratio = max(threshs,[],2)./min(threshs,[],2);
latratio = max(lats,[],2)./min(lats,[],2);

figure; 
subplot(2,2,1);
hist(threshratio);
title('Threshold ratio');
subplot(2,2,2);
hist(latratio);
title('Latency ratio');
subplot(2,2,3); hold on;
for i = 1:length(fnames)
    h(i) = plot(latratio(i), threshratio(i),'ko','MarkerFaceColor','black');
    set(h(i),'ButtonDownFcn',['disp(''', char(fnames{i}(1)),''')']);
end
xlabel('latency ratio');
ylabel('threshold ratio');

% Seeing which files to cull based on # of color directions tested
L = nanmin(ncolordirstested,[],2) < 20;
cellfun (@(x) (x(2)), fnames(L));
L = nanmin(nnoogs,[],2) < 3;
cellfun (@(x) (x(2)), fnames(L))

% Which cells to cull based on differences in theshold?
figure; axes; hold on;
plot(max(threshs,[],2)./min(threshs,[],2),'k.')
plot([0 size(threshs,1)],[1.5 1.5])

% Looking for files that were abbreviated
fnames(any(ncolordirstested < 20,2))
fnames{find(nnoogs < 3)}
%%
% Section 6) Looking through the directory for Neurothresh files that 
% were tested with alternative color directions.  Also pulling out other
% parameters from the header

whichfiles = 'TEXTFILE'; % 'ALL or 'TEXTFILE';
if strcmp(whichfiles, 'ALL') 
    datapath = 'N:\NexFiles';
    dirstruct = subdir([datapath filesep '*.nex']);
    nfiles = size(dirstruct,1);
elseif strcmp(whichfiles, 'TEXTFILE')
   [fnames, spikeIdx] = fnamesFromTxt2();
   fnames = cellfun(@(x) x(:,2), fnames);
   nfiles = size(fnames,1);
end

data = {};
for i = 1:nfiles
    i
    if strcmp(whichfiles, 'ALL')
        try
            filename = dirstruct(i).name;
            stro = nex2stro(filename);
        catch exception
            wh = findall(0, 'tag', 'TMWWaitbar');
            if ~isempty(wh), delete(wh); end
            continue
        end
    elseif strcmp(whichfiles, 'TEXTFILE')
        filename = char(fnames{i});
        stro = nex2stro(findfile(filename));
    end
    if (stro.sum.paradigmID ~= 103)
        continue;
    end
    idx = length(data)+1;
    
    data(idx).filename = filename;
    data(idx).monkinitial = filename(1);
    data(idx).rf = [stro.sum.exptParams.rf_x stro.sum.exptParams.rf_y];
    data(idx).nreversals = stro.sum.exptParams.nreversals;
    data(idx).stepsize =  stro.sum.exptParams.stepsize;
    data(idx).scale = stro.sum.exptParams.scale;
    data(idx).linpredtol = stro.sum.exptParams.linpredtol;
    data(idx).oogscale = stro.sum.exptParams.oogscale;
    data(idx).initcoldirs = [stro.sum.exptParams.initLMS1 stro.sum.exptParams.initLMS2 stro.sum.exptParams.initLMS3];
end

initcoldirs = cat(3,data.initcoldirs);
altfnames = [];
stdfnames = [];
for i = 1:size(data,2)
    i
    if all(all(data(i).initcoldirs ~= [0.025 -0.01 0; 0.025 0.01 0; 0 0 0.1]))
        if (data(i).initcoldirs ~= [0.025 -0.01 0.01; 0.025 0.01 -0.01; 0.025 0.1 0.1])
            shortfn = data(i).filename(find(data(i).filename == '\',1,'last')+1:find(data(i).filename == '.',1,'last')-1);
            altfnames = [altfnames; shortfn];
        else
            data(i).initcoldirs
        end
    else
        shortfn = data(i).filename(find(data(i).filename == '\',1,'last')+1:find(data(i).filename == '.',1,'last')-1);
        stdfnames = [stdfnames; shortfn];
    end
end

% How many files from each monkey?
sum([data.monkinitial] == 'K')
sum([data.monkinitial] == 'S')

% nreversals
[data.nreversals]
% stepsize x(i+1)=x(i)*(1+/-stepsize)
[data.stepsize]
% How much the step size changes after each reversal
[data.scale]
% Linear prediction tolerance
[data.linpredtol]
% OOG scale
[data.oogscale]

%%
% Section 7 
% Finding NeuroThresh planar cells with L+M cone weights.  Also ellipsoidal
% cells.
[fnames, spikeIdx] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\NT\Gratings&NT2.txt');
data ={};
for cellcounter = 1:size(fnames,1)
    NT = {};
    for i = 1:size(fnames{cellcounter},2)
        filename = findfile(char(fnames{cellcounter}(i)));
        paradigmID = getparadigmID(filename);
        if paradigmID == 103
            NT = nex2stro(filename);
        end
    end
    
    % -------------------------------
    % Converting cone contrasts in nex file to 10 deg fundmentals.
    % Must go through excitations first - can't just transform contrasts.
    % -------------------------------
    lmsidxs = [find(strcmp(NT.sum.trialFields(1,:),'lcont'))...
        find(strcmp(NT.sum.trialFields(1,:),'mcont'))...
        find(strcmp(NT.sum.trialFields(1,:),'scont'))];
    fundamentals = NT.sum.exptParams.fundamentals;
    fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
    mon_spd = NT.sum.exptParams.mon_spd;
    mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
    mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
    M = fundamentals'*mon_spd;
    load ('T_cones_smj10');
    M10 = T_cones_smj10*mon_spd;
    NT.trial(:,lmsidxs) = ConvertConeContrastBasis(M, M10, NT.sum.exptParams.bkgndrgb, NT.trial(:,lmsidxs));
    % -------------------------------
    out = NTpreprocess(NT,0,Inf);
    
    scaled = out(:,[2:4]).*repmat(out(:,5), 1,3);
  
    Loog = logical(out(:,7));
    [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(scaled, Loog);
   
    coneweights = (xformmat*planeparams)';
    if (coneweights(2) < 0)
        coneweights = -coneweights;
    end
    data(cellcounter).coneweights = coneweights./sum(abs(coneweights));
    A = [quadparams(1) quadparams(4) quadparams(5);...
         quadparams(4) quadparams(2) quadparams(6);...
         quadparams(5) quadparams(6) quadparams(3)];

    % [evecs,evals] = eig(A); % Doing this analysis in the whitened space. Is this OK?
    % This space used *does* make a difference to the magnitudes (not
    % signs) of the eigenvalues. GDLH 3/14/12
    [evecs,evals] = eig(xformmat*A*xformmat'); % Back to the cone space
    [evals,i] = sort(diag(evals),1,'ascend');
    data(cellcounter).evecs = evecs(:,i);
    data(cellcounter).evals = evals;
    data(cellcounter).Loog = sum(Loog);
    data(cellcounter).filename = NT.sum.fileName(find(NT.sum.fileName == '\',1,'last')+1:find(NT.sum.fileName == '.',1,'last')-1);
    data(cellcounter).planeSSE = planeSSE;
    data(cellcounter).quadSSE = quadSSE;
end

EIGENTHRESH = 0.5;
SSERATIOTHRESH = 0;
CONEWEIGHTTHRESH = 0.9;
LOOGTHRESH = 5;
cw = reshape([data.coneweights],3,length(data))';
L = cw(:,1)+cw(:,2) > CONEWEIGHTTHRESH;
%L = L & cw(:,1) > .25 & cw(:,2) > .25;
evals = [data.evals];

normevals = sort(evals./repmat(max(abs(evals)),[3 1]),1);

%figure; axes; hold on;
%plot(normevals(2,:),'k.');
% A criterion based on eigenvalues

L = L & (normevals(3,:) == 1)'; 
L = L & (max(abs(normevals([1 2],:))) < EIGENTHRESH)';
L = L & [data.Loog]' >= LOOGTHRESH;  % make sure no overlap with ellipsoids

% a criterion based on SSE of model fits
SSEratio = [data.quadSSE]./[data.planeSSE];
figure; axes; hold on;
%plot(SSEratio,'k.');
plot(SSEratio,max(abs(normevals([1 2],:))),'k.');
xlabel('SSE ratio'); ylabel('2nd eigenvector');

L = L & SSEratio' > SSERATIOTHRESH;
plot(SSEratio(L),max(abs(normevals([1 2],L))),'r*');

figure; axes; hold on;
plot(cw(:,1),cw(:,2),'k.');
plot(cw(L,1),cw(L,2),'g.');
for i = find(L)'
    disp(data(i).filename)
end
% 
% % Finding ellipsoids
% L = [data.Loog]' <=4;
% L = L & [sum(sign([data.evals]))]' == 3;
% for i = find(L)'
%     disp(data(i).filename)
% end
% sum(L)
%%
% Section 8
% Finding DTspot files with aberrant thresholds

filenames = fnamesFromTxt2();
data = [];
load T_xyz1931;
Vlambda = T_xyz1931(2,:);
for fileidx = 1:size(filenames,1)
    whichfiles = filenames{fileidx};
    for i = 1:length(whichfiles)
        stro = nex2stro(findfile(char(whichfiles{i})));
        TFidx = find(strcmp(stro.sum.trialFields(1,:),'gabor_speed'));
        [thresh, colorDirs, sfs] = DTquestUnpack(stro, 'mode');
        pause
        close(gcf); % for when you just want the summary figures
        normcolordirs= colorDirs./repmat(sqrt(sum(colorDirs.^2,2)),1,3);
        guns = reshape(stro.sum.exptParams.mon_spect,81,3);
        bkgndrgb = [stro.sum.exptParams.bkgnd_r stro.sum.exptParams.bkgnd_g stro.sum.exptParams.bkgnd_b];
        M = reshape(stro.sum.exptParams.m_mtx,3,3);
        bkgndlms = M*bkgndrgb';
        bkgndlum = sum(bkgndlms([1 2]));
        % Quantifying thresholds in 2 deg luminance contrast units
        bkgndlum = bkgndrgb*guns'*Vlambda';
        threshrgb = inv(M)*(repmat(bkgndlms',2,1).*(1+repmat(thresh/100,[1 3]).*normcolordirs([1 2],:)))';
        deltagun = threshrgb-repmat(bkgndrgb',1,2);
        for j = 1:size(deltagun,2)
            thresholds(1,j) =  (guns*deltagun(:,j))'*Vlambda';
        end
        ecc = stro.sum.exptParams.rf_x;
        data = [data; ecc/10 thresholds'];
    end
end
% 
% % Iterative approach, chucking data points, one at a time, whose 
% % residuals are > 3 stdevs from the fit.
% whichfilesbad = false(length(filenames),1);
% stop = 0;
% tmpdata = data;
% while (stop == 0)
%     figure(1); clf; axes; hold on;
%     L = ~whichfilesbad;
%     plot(tmpdata(L,1),tmpdata(L,2)*100,'go','MarkerFaceColor','green','MarkerSize',3);
%     plot(tmpdata(L,1),tmpdata(L,3)*100,'bo','MarkerFaceColor','blue','MarkerSize',3);
%     pp_g  = csaps(tmpdata(L,1),log10(tmpdata(L,2)*100),.5);
%     pp_b  = csaps(tmpdata(L,1),log10(tmpdata(L,3)*100),.5);
%     x = linspace(0,8,100);
%     plot(x,10.^(fnval(x,pp_g)),'g-','Linewidth',2)
%     plot(x,10.^(fnval(x,pp_b)),'b-','Linewidth',2)
%     
%     % Looking at the residuals
%     resid_g = 10.^(fnval(tmpdata(:,1),pp_g))-tmpdata(:,2)*100;
%     resid_b = 10.^(fnval(tmpdata(:,1),pp_b))-tmpdata(:,3)*100;
%     resid_g(whichfilesbad) = 0;
%     resid_b(whichfilesbad) = 0;  
%     
%     all_resid = [resid_g;resid_b];
%     figure(2); clf; hist(all_resid);
%     CI = mean(all_resid)+[-1 1]*3*std(all_resid);
%     badfileidx = find(resid_g<CI(1) | resid_g > CI(2)|resid_b<CI(1) | resid_b > CI(2));
%     if (isempty(badfileidx))
%         stop = 1;
%     else
%         worstresid = max(max(abs([resid_g(badfileidx) resid_b(badfileidx)])));
%         worstidx = find(sum(abs([resid_g(badfileidx) resid_b(badfileidx)]) == worstresid,2)); 
%         whichfilesbad(badfileidx(worstidx)) = 1;
%     end
% end
%disp('Here are the bad files');
%disp(char([filenames{whichfilesbad}]'))

% Just chucking the data points with residuals in the top 95%
figure(1); clf; axes; hold on;
plot(data(:,1),data(:,2)*100,'go','MarkerFaceColor','green','MarkerSize',3);
plot(data(:,1),data(:,3)*100,'bo','MarkerFaceColor','blue','MarkerSize',3);
pp_g  = csaps(data(:,1),log10(data(:,2)*100),.5);
pp_b  = csaps(data(:,1),log10(data(:,3)*100),.5);
x = linspace(0,8,100);
plot(x,10.^(fnval(x,pp_g)),'g-','Linewidth',2)
plot(x,10.^(fnval(x,pp_b)),'b-','Linewidth',2)

% Looking at the residuals
resid_g = 10.^(fnval(data(:,1),pp_g))-data(:,2)*100;
resid_b = 10.^(fnval(data(:,1),pp_b))-data(:,3)*100;
all_resid = [resid_g;resid_b];
resid_threshold = prctile(abs(all_resid),95);
Lbadfiles = abs(all_resid) > resid_threshold;
whichfilesbad = any(reshape(Lbadfiles,size(resid_g,1),2),2);
disp('Here are the bad files');
disp(char([filenames{whichfilesbad}]'))

%%
% Section 8.1
% Culling DTscot files.
if(ispc)
    filelistpath = 'N:\NexFiles\nexfilelists\Greg\DTEM';
else
   filelistpath = '/Volumes/NO BACKUP/NexFiles/nexfilelists/Greg/DTEM';
end

listnames = {'ZackDTscot.txt','GregDTscot.txt','LeahDTscot.txt','FreyaDTscot.txt','SednaDTscot.txt' };
subjectnames = {'Human Z','Human G','Human L','Monkey F','Monkey S'};
%listnames = {'GregDTscot.txt'};
%subjectnames = {'Human G'};

data = []; allfnames = [];
for i = 1:length(listnames)
    [fnames, spikenums] = fnamesFromTxt2([filelistpath,filesep,listnames{i}]);
    allfnames = [allfnames; fnames];
    for j = 1:length(fnames)
        stro = nex2stro(findfile(char(fnames{j})));
        
        colordir_idx = strcmp(stro.sum.trialFields(1,:), 'color_dir');
        questmode_idx = strcmp(stro.sum.trialFields(1,:), 'quest_mode');
        thresholds = zeros(3,1);
        
        for gun = 1:3
            L = stro.trial(:, colordir_idx) == gun;
            all_modes = stro.trial(L, questmode_idx);
            subplot(3,1,gun);
            plot(all_modes,'k.');
            thresholds(gun) = all_modes(end);
        end
        data = [data; thresholds' i];
        title(subjectnames{i});
        %pause
    end
end

% Some helpful plots plus chucking files for which any one of the thresholds 
% is far from:  
% (a) median (5% of the residuals, relative to the median, log scaled)
% (b) mean (anything 2.5 standard deviations or more away from the mean)

for i = 1:length(listnames)
    figure; axes; hold on;
    Lobs = data(:,end) == i;
    logthresh = log10(data(Lobs,[1 2 3]));
    plot(logthresh(:,1),'r.');
    plot(logthresh(:,2),'g.');
    plot(logthresh(:,3),'b.');
    title(subjectnames{i}); xlabel('filenumber'); ylabel('log threshold');
    set(gca,'Ylim',[-2.5 0],'Xlim',[0 sum(Lobs)+1]);
    obsfnames = allfnames(Lobs);

   residuals = logthresh - repmat(median(logthresh), size(logthresh,1),1);
   L = abs(residuals(:)) > prctile(abs(residuals(:)),95);
   L = reshape(L,size(logthresh));
    
 %  L = logthresh<repmat(mean(logthresh)-2*std(logthresh),sum(Lobs),1) |...
 %       logthresh>repmat(mean(logthresh)+2*std(logthresh),sum(Lobs),1);
    
    whichfilesbad = any(L,2);
  %  plot(logthresh(whichfilesbad,:)','ro');
    disp(['Here are the bad files for ',char(listnames{i})]);   
    disp(char([obsfnames{whichfilesbad}]'));
end


%%
% Section 9
% Looking at the temporal frequency of DTNT files

if (ismac)
    listpath = '/Volumes/NO BACKUP/NexFiles/nexfilelists/Greg';
else
    listpath = 'N:\NexFiles\nexfilelists\Greg';
end
%[~,txt,raw]=xlsread([listpath,filesep,'DTNT', filesep, 'DTNT25Hz']);
%filenames = txt(:,3);
 
[filenames, spikenums] = fnamesFromTxt2([listpath,filesep,'DTEM', filesep,'KaliMacPig.txt']);

data = [];
for a = 1:size(filenames,1)
    if (length(char(filenames{a})) > 6)   % Only taking real files
        stro = nex2stro(findfile(char(filenames{a})));
        TF = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'gabor_speed'));
        data = [data; min(TF) max(TF) stro.sum.exptParams.rf_x]
    end
end

unique(data(:,3))


%%
% Section 10
% Finding SMurray files (size constancy experiment) in which the gray
% background or the corridor background were used.
PARADIGMID = 102;
if (ispc)
    datapath = 'N:\NexFiles\Amy\Apollo';
else
    datapath = '/Volumes/NO BACKUP/NexFiles/Amy/Apollo';    
end

eval(['cd (''',datapath,''')']);
startdate = '010113';  % month, day, year (last two digits)
startdatenum = datenum(str2double(['20',startdate([5 6])]), str2double(startdate([1 2])), str2double(startdate([3 4])));
% First making a gigantic dirstruct
dirstruct = dir;
for i = 1:size(dirstruct,1)
    if (dirstruct(i).isdir && ~strncmp(dirstruct(i).name, '.',1))
        cd(dirstruct(i).name);
        dirstruct = [dirstruct; dir];
        cd('../')
    end
end
counter = 1;
filenames = {};
bkgnds = [];
for i = 1:size(dirstruct,1)
    if (length(dirstruct(i).name) < 4)
        continue;
    end
    [num2str(i),' ',dirstruct(i).name]
    if (strcmp(dirstruct(i).name(end-3:end),'.nex'))
        if (length(dirstruct(i).name) >= 13)  % minimum # of characters for file named "<prefix><date><number>.nex"
            % Check the date
            filedate = dirstruct(i).name(end-12:end-7);
            filedatenum = datenum(str2double(['20',filedate([5 6])]), str2double(filedate([1 2])), str2double(filedate([3 4])));
            if (filedatenum <= startdatenum)
                continue
            end
        end
        try
            stro = nex2stro(findfile(dirstruct(i).name));
        catch
            disp(['skipping a file: ',dirstruct(i).name]);
            continue;
        end
        if ~isfield(stro,'sum')
             disp(['skipping a file: ',dirstruct(i).name]);
            continue;
        end
        if ~(stro.sum.paradigmID == PARADIGMID)
             disp(['skipping a file: ',dirstruct(i).name]);
            continue;
        end
        filenames{counter} = dirstruct(i).name;
        bkgnds(counter) = unique(stro.trial(:,strcmp(stro.sum.trialFields(1,:),'bkgnd')));
        disp([num2str(counter),': Got one: ',dirstruct(i).name,' bkgnd = ',num2str(bkgnds(counter))]);
        counter = counter+1;
    end
end

%%
% Section 10.1
% Getting some key statistics from Amy's V1 size illusion data files
% (SMurray). rfx always equaled rfy? Thats doesn't make sense.

if (ismac)
    nexfilepath = '/Volumes/NO BACKUP/NexFiles/Amy';
    filelist = '/Volumes/NO BACKUP/NexFiles/nexfilelists/Greg/AmyFreya.txt';
else
    nexfilepath = 'N:\NexFiles\Amy';
    filelist = 'N:\NexFiles\nexfilelists\Greg\AmyApollo.txt';
end

filenames = fnamesFromTxt2(filelist);

daylist = [];
currentdate = [];
for i = 1:length(filenames)
    tmp = char(filenames{i});
    if (i == 1)
        daylist = 1;
        currentdate = str2num(tmp(2:end-7));
    else
        if (str2num(tmp(2:end-7)) == currentdate)
            daylist(i) = daylist(i-1);
        else
            daylist(i) = daylist(i-1)+1;
            currentdate = str2num(tmp(2:end-7));
        end
    end
end


data = [];
for i = 1:length(filenames)
    stro = nex2stro(findfile(char(filenames{i})));
    ntrials = size(stro.trial,1);
    rfx = stro.sum.exptParams.rf_x;
    rfy = stro.sum.exptParams.rf_y;
    mondist = stro.sum.exptParams.mondist;
    bkgnd = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'bkgnd'));
    ubkgnd = unique(bkgnd);
    if (length(ubkgnd) > 1)
        keyboard
    end
    data = [data; daylist(i) ubkgnd rfx rfy ntrials mondist]
end

figure; axes; hold on;
plot(data(:,3)/10,data(:,4)/10,'ko')
plot(0,0,'k*')

% Getting ring sizes (mostly to compare psychophysics data sets between
% Apollo and Freya)

ir_near = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_ir_near'));
ir_far = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_ir_far'));

figure; 
subplot(2,1,1); hist(2*ir_near);
subplot(2,1,2); hist(2*ir_far);

%%
% Section 11
% How many trials per threshold measurement

%humanlist = {'GregDTNT15Hz.txt','ZackDTNT15Hz.txt','LeahDTNT15Hz.txt'};
%monkeylist = {'NutDTNT15Hz.txt','FreyaDTNT15Hz.txt','SednaDTNT15Hz.txt'};

humanlist = {'GregMacPig.txt','ZackMacPig.txt','LeahMacPig.txt'};
monkeylist = {'NutMacPig.txt','FreyaMacPig.txt','SednaMacPig.txt','KaliMacPig.txt','ApolloMacPig.txt'};

if (ismac)
    filelistpath = '/Volumes/NO BACKUP/NexFiles/nexfilelists/Greg/DTNT';
else
%    filelistpath = 'N:\NexFiles\nexfilelists\Greg\DTNT';
    filelistpath = 'N:\NexFiles\nexfilelists\Greg\DTEM';
end

data = [];
for HUMANS = 0:1
    if (HUMANS)
        listnames = humanlist;
    else
        listnames = monkeylist;
    end
    for i = 1:length(listnames)
        filenames = fnamesFromTxt2([filelistpath,filesep,char(listnames{i})]);
        for fileidx = 1:size(filenames,1)
            if (strncmp(filenames{fileidx},'sf',2))
                continue;
            end
            stro = nex2stro(findfile(char(filenames{fileidx})));
            [thresholds, colordirs, sfs, QuestTrajectories] = DTquestUnpackGH(stro, 'mode');
            data = [data; i size(QuestTrajectories{1},2)];
        end
    end
end

%%
% Section 12
% In LMTF data files, how long was the stimulus on the screen for? There
% are likely some problems here since the old version of LMTFslave.m
% accepted numbers of frames from REX, and these weren't changed when the
% refresh rate of the monitor changed.

textfilepath = [nexfilepath,filesep,'nexfilelists',filesep,'Greg',filesep,'LMTF'];
listnames = {'BigApolloTmpListLMTF'};
data = [];
for textfileidx = 1:length(listnames)
    filenames = fnamesFromTxt2([textfilepath,filesep,char(listnames{textfileidx})]);
    for nexfileidx = 1:length(filenames)
        try
            stro = nex2stro(findfile(char(filenames{nexfileidx})));
            stimonidx = strcmp(stro.sum.trialFields(1,:),'stimon_t');
            stimoffidx = strcmp(stro.sum.trialFields(1,:),'stimoff_t');
            stimdur = stro.trial(:,stimoffidx)-stro.trial(:,stimonidx);
            data = [data; textfileidx nexfileidx mode(stimdur)];
        catch
            disp(sprintf('problem with file %s',char(filenames{nexfileidx})));
            continue;
        end
    end
end

%%
% Section 13
% Getting RF positions for a bunch of LMTF data files
% Somehow link to Excel

textfilepath = [nexfilepath,filesep,'nexfilelists',filesep,'Greg',filesep,'LMTF'];
listnames = {'BigApolloTmpListLMTF'};
data = [];
for textfileidx = 1:length(listnames)
    filenames = fnamesFromTxt2([textfilepath,filesep,char(listnames{textfileidx})]);
    for nexfileidx = 1:length(filenames)
        try
            stro = nex2stro(findfile(char(filenames{nexfileidx})));
            rfx = abs(stro.sum.exptParams.stim_x); % For LMTF don't care about side of the screen
            rfy = stro.sum.exptParams.stim_y;         
            targon_t_idx = strcmp(stro.sum.trialFields(1,:),'targon_t');
            saccmade_t_idx = strcmp(stro.sum.trialFields(1,:),'saccmade_t');
            RT = stro.trial(:,saccmade_t_idx)-stro.trial(:,targon_t_idx);
            data = [data; textfileidx nexfileidx rfx rfy mean(RT) stro.sum.paradigmID];
        catch
            disp(sprintf('problem with file %s',char(filenames{nexfileidx})));
            continue;
        end
    end
end
% order of columns in data
% 1) Text file id
% 2) NEX file id
% 3) rfX
% 4) rfY
% 5) Mean saccade latency
% 6) Paradigm ID (should always be 157)

assert(all(data(:,6) == 157),'Non LMTF files found: %s', char([filenames{data(:,6) ~= 157}])');
    
uniqueRFxy = unique(data(:,[3 4]),'rows','stable');
ns = zeros(size(uniqueRFxy,1),1);
RTs = zeros(size(uniqueRFxy,1),1);
for i = 1:size(uniqueRFxy,1)
    L = data(:,3)==uniqueRFxy(i,1) & data(:,4)==uniqueRFxy(i,2);
    ns(i) = sum(L);
    RTs(i) = mean(data(L,5));
    if ~all(diff(find(L))==1) % Making sure the files with a common RF and adjacent in the file
        disp('Some files out of order');
        sprintf('Files with eccentricity %d, %d include', uniqueRFxy(i,:))
        [filenames{L}]'
    end
end

% Figure showing how much data we have at each location
figure; axes; hold on;
for i = 1:size(uniqueRFxy,1)
    h = plot(uniqueRFxy(i,1),uniqueRFxy(i,2),'ko','MarkerFaceColor','black');
    set(h,'MarkerSize',20*ns(i)/max(ns))
end

% Figure showing average latency at each location
figure; axes; hold on;
for i = 1:size(uniqueRFxy,1)
    h = plot(uniqueRFxy(i,1),uniqueRFxy(i,2),'ko','MarkerFaceColor','black');
    set(h,'MarkerSize',20*RTs(i)/max(RTs))
end


%%
% Section 14
% Isosamp. Culling neurons on the basis of their RF locations
% and number of trials per condition.
NEARTHRESHOLD = 2;
FARTHRESHOLD = 14;
YTHRESHOLD = 2;

conn = database('Nex_Paradigm_Sort','','','Vendor','MySql','Server','128.95.153.12');
flistplus = fetch(conn, 'SELECT fileID, rfX, rfY, neuron FROM IsoSamp_LGN WHERE quality = ''1''');
close(conn);
data = [];
for i = 1:size(flistplus,1)
   rfx = flistplus{i,2};
   rfy = flistplus{i,3};
   data = [data; rfx/10 rfy/10];
end
ecc = sqrt(sum(data.^2,2));
Ltooclosetofovea = ecc < NEARTHRESHOLD;
Ltoofarfromfovea = ecc > FARTHRESHOLD;
Ltooclosetovertmerdian = abs(data(:,1)) < YTHRESHOLD;

sum(Ltooclosetofovea)
sum(Ltoofarfromfovea)
sum(Ltooclosetovertmerdian)

% For now ignoring distance to vertical meridian
% Getting a list of the cell whose RFs are too close to the fixation point
if any(Ltooclosetofovea)
    disp('These cells have RFs that are too close to the fovea');
    for i = find(Ltooclosetofovea)'
        disp(flistplus{i,1})
    end
end

% Getting a list of the cell whose RFs are too far to the fixation point
if any(Ltoofarfromfovea)
    disp('These cells have RFs that are too far to the fovea');
    for i = find(Ltoofarfromfovea)'
        disp(flistplus{i,1})
    end
end

% Culling on the basis of number of trials per condition
uniqueids = unique([flistplus{:,4}]);
data=[];
for i = 1:length(uniqueids)
    L = [flistplus{:,4}] == uniqueids(i);
    stro =[];
    for j=find(L)
        stro = strocat(stro,nex2stro(findfile(flistplus{j,1},[nexfilepath,filesep,'Greg'])));
    end
    Lcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_l'));
    Mcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_m'));
    TF = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'tf'));
    uniquetrialtypes = unique([Lcc Mcc TF],'rows');
    tmpdata = zeros(size(uniquetrialtypes,1),1);
    for j = 1:length(tmpdata)
       L = Lcc == uniquetrialtypes(j,1) & Mcc == uniquetrialtypes(j,2) & TF == uniquetrialtypes(j,3);
       tmpdata(j) = sum(L);
    end
    data = [data; min(tmpdata) max(tmpdata) length(tmpdata)];
    if min(tmpdata) < 5
        keyboard
        disp(stro.sum.fileName)
    end
end

%%
% Section 15
% Getting fixation window and fixation point sizes
% and max TF.
conn = database('Nex_Paradigm_Sort','','','Vendor','MySql','Server','128.95.153.12');
flist = fetch(conn, 'SELECT fileID neuron FROM IsoSamp_LGN WHERE quality = ''1''');
close(conn);
data = [];
for i = 1:size(flist,1)
    stro = nex2stro(findfile(flist{i}));
    tf = max(stro.trial(:,strcmp(stro.sum.trialFields(1,:), 'tf')));
    data = [data; stro.sum.exptParams.eyewin_x stro.sum.exptParams.eyewin_y stro.sum.exptParams.fp_size tf];
end

%%
% Section 16
% Getting the number and sizes of stixels used in the white noise
% experiments.

[fnames_isosamp,spikecds,~,neuronid_isosamp] = fnamesFromTxt('IsoSamp_LGN','cellClass',{'M','P'});

data = [];
for i = neuronid_isosamp' 
    fname_WN = fnamesFromTxt('WhiteNoiseLGN_forIS','neuron',{num2str(i)});
    stro = nex2stro(findfile(fname_WN{1}{1}));
    data = [data; stro.sum.exptParams.nstixperside stro.sum.exptParams.npixperstix];
end

% Propixx screen is 51 cm wide (rigconts.h) and viewed from 60-61 cm away.
% Propixx screen is 1024 pixels wide (double check!)
unq_stim = unique(data,'rows'); % 
for i = 1:size(unq_stim,1)
    L = data(:,1) == unq_stim(i,1) & data(:,2) == unq_stim(i,2);
    ncells = sum(L);
    stixwidthincm = unq_stim(i,2)/1024*51;
    stixwidthindeg = atan2(stixwidthincm,61)*180/pi;
    sprintf('%d cells used %d x %d stimuli with %d° stixels',ncells,unq_stim(i,1),unq_stim(i,1),stixwidthindeg)
end


