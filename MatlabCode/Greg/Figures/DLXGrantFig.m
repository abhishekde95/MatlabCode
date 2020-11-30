% DLX in V1 grant figs
%
% For DToneloc interleaved spatial control, see Abhishek's "Popfigures2.m",
% Figure 4.
%
% Section 1: SC_stimcue showing saccades to a row of target locations using
% datafile M051518001.1 (also can be used to show four directions of
% saccade)
%
% Section 1.1 Old FixStim data for comparison with SC_stimcue data.
%
% Section 2: DToneloc data from M0518180XX showing some spread in the
% spatial distribution of optoeffect. Only picking three target locations.
%
% Section 2.1: Like section 2 but using a single file that had multiple
% randomly interleaved stimulus locations
%
% Section 3: A single DToneloc example data file.
%
% Section 4: Pseudo ROC curve showing an effect of the laser on the noise
%
% Section 5: Schematic of the pseudo ROC curve with multiple criteria.
%
% Section 6: Spatial uncertainty simulation.
%
% Section 7: Achromatic and S-cone detection supressed by laser
%
% Section 8: SC_stimcue fig for Sackler talk
%
% Section 9: Rasters for example excited and suppressed cells
%%
% Section 1: SC_stimcue data. A row of targets showing the spatial
% specificity of the opto effect.

stro = nex2stro(findfile('M042618002.nex',[nexfilepath,filesep,'Abhishek']));

% M042618002

targ_x = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targ_x'));
targ_y = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targ_y'));
uniquetargxy = unique([targ_x targ_y],'rows');
ntrials = size(stro.trial,1);
Lcatchtrials = targ_x == 0 & targ_y == 0;
laseron_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimon_t'));
laseroff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimoff_t'));
Llaser = ~isnan(laseron_t);
fpoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fpoff_t'));
fix_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fpacq_t'));
targon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targon_t'));
targoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targoff_t'));
sacinit_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'saccinit_t'));
sacmade_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'saccmade_t'));
H = stro.ras(:,strcmp(stro.sum.rasterCells,'AD11'));
V = stro.ras(:,strcmp(stro.sum.rasterCells,'AD12'));
ep_t = stro.ras(:,strcmp(stro.sum.rasterCells,'anlgStartTime'));
ADfreq = stro.sum.analog.storeRates{1};

%uniquetargxy(uniquetargxy(:,2) >= 0,:) = []; % Getting rid of targets in the upper hemifield

%middlerow = median(uniquetargxy(:,2));
%middlerow = median(uniquetargxy(:,2));
%middlecolumn = median(uniquetargxy(:,1));
%uniquetargxy(uniquetargxy(:,2) ~= middlerow,:) = [];% Getting rid of targets not in the middle row
%uniquetargxy(uniquetargxy(:,2) > middlerow,:) = [];% Getting rid of targets in the upper row

L = logical([0 1 1 0 0 0 1 1 0]);
uniquetargxy = uniquetargxy(L,:);

data = [];
for i = 1:size(uniquetargxy,1)
   Ltarg = targ_x == uniquetargxy(i,1) & targ_y == uniquetargxy(i,2);
   lat_laser = sacinit_t(Ltarg&Llaser)-fpoff_t(Ltarg&Llaser);
   lat_nolaser = sacinit_t(Ltarg&~Llaser)-fpoff_t(Ltarg&~Llaser);
   data(i,:) = [nanmean(lat_laser) sum(~isnan(lat_laser)) nanmean(lat_nolaser) sum(~isnan(lat_nolaser)) ];
end

barwidth = 2.1;
maxbarheight = 6;
barxoffset = 2.4;
figure; axes; hold on;
for i = 1:size(uniquetargxy,1)
   if data(i,2) > 0 % Laser trials
       h = rectangle('Position',[uniquetargxy(i,1),uniquetargxy(i,2) barwidth maxbarheight/max(max(data(:,[2 4])))*data(i,2)]);
       set(h,'FaceColor','blue','EdgeColor','black')
   end
   if data(i,4) > 0 % No laser trials
       h = rectangle('Position',[uniquetargxy(i,1)+barxoffset,uniquetargxy(i,2) barwidth maxbarheight/max(max(data(:,[2 4])))*data(i,4)]);
       set(h,'FaceColor',[.5 .5 .5])
   end
   %set(h,'EdgeColor','white');
end
% Adding a fixation point
plot(0,0,'k+','MarkerSize',15,'LineWidth',3);
axis equal;

% Plots in X,Y
figure; axes; hold on;
for i = 1:size(uniquetargxy,1)
    Ltarg = targ_x == uniquetargxy(i,1) & targ_y == uniquetargxy(i,2);
    trialsid = find(Ltarg);
    for counter = 1:sum(Ltarg)
        trialid = trialsid(counter);
        h = H{trialid}*4096/400;
        v = V{trialid}*4096/400;
        t = [0:length(h)-1]./ADfreq+ep_t{trialid};
        t = t-targon_t(trialid);
        Lt = t>0 & t<.30;
        if ~Llaser(trialid)
            plot(h(Lt),v(Lt),'-','Color',[0 0 0],'Linewidth',2);
        else
            plot(h(Lt),v(Lt),'-','Color',[0 1 1 ],'Linewidth',2);
        end
    end
end
axis equal;

% Now a few plots of (X/Y vs t)
figure;
for i = 1:size(uniquetargxy,1)
    Ltarg = targ_x == uniquetargxy(i,1) & targ_y == uniquetargxy(i,2);
    trialsid = find(Ltarg);
    subplot(length(unique(uniquetargxy(:,2))),length(unique(uniquetargxy(:,1))),i); hold on;
    
    for counter = 1:sum(Ltarg)
        trialid = trialsid(counter);
        h = H{trialid}*4096/400;
        v = V{trialid}*4096/400;
        t = [0:length(h)-1]./ADfreq+ep_t{trialid};
        t = t-targon_t(trialid);
        Lt = t>0.0 & t<.3;
        if ~Llaser(trialid)
            plot(t(Lt),h(Lt),'-','Color',[0 0 0],'Linewidth',2);
            plot(t(Lt),v(Lt),'-','Color',[.7 .7 .7],'Linewidth',2);
        else
            plot(t(Lt),h(Lt),'--','Color',[0 0 1],'Linewidth',2);
            plot(t(Lt),v(Lt),'--','Color',[0 .9 .9],'Linewidth',2);
        end
    end
    set(gca,'Ylim',[-8 8],'Xlim',[0 .33]);   
end

%%
% Section 1.1 
% Some old Fixstim data for comparison with the SC_stimcue data
stro = nex2stro(findfile('S081111023.nex'));
% S081111018
stim_type = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_type'));
targ_shown = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targ_shown'));
fpoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fpoff_t'));
H = stro.ras(:,strcmp(stro.sum.rasterCells,'AD11'));
V = stro.ras(:,strcmp(stro.sum.rasterCells,'AD12'));
ep_t = stro.ras(:,strcmp(stro.sum.rasterCells,'anlgStartTime'));
ADfreq = stro.sum.analog.storeRates{1};
Llaser = logical(stim_type);
Ltarget = logical(targ_shown);

% Plots in X,Y
figure; axes; hold on;
for i = 0:1
    L = Llaser == i & ~Ltarget;
    trialsid = find(L);
    for counter = 1:sum(L)
        trialid = trialsid(counter);
        h = H{trialid}*4096/400;
        v = V{trialid}*4096/400;
        t = [0:length(h)-1]./ADfreq+ep_t{trialid};
        t = t-fpoff_t(trialid);
        Lt = t>0 & t<.30;
        if ~Llaser(trialid)
            plot(h(Lt),v(Lt),'-','Color',[0 0 0],'Linewidth',2);
        else
            plot(h(Lt),v(Lt),'-','Color',[0 1 1 ],'Linewidth',2);
        end
    end
end
axis square;
set(gca,'Xlim',[-4 4],'Ylim',[-4 4]);


% Now a few plots of (X/Y vs t)
figure;
for i = 0:1
    L = Llaser == i & ~Ltarget;
    trialsid = find(L);
    subplot(2,2,i+1); hold on;
    for counter = 1:sum(L)
        trialid = trialsid(counter);
        h = H{trialid}*4096/400;
        v = V{trialid}*4096/400;
        t = [0:length(h)-1]./ADfreq+ep_t{trialid};
        t = t-fpoff_t(trialid);
        Lt = t>0.0 & t<.3;
        if ~Llaser(trialid)
            plot(t(Lt),h(Lt),'-','Color',[0 0 0],'Linewidth',2);
            plot(t(Lt),v(Lt),'-','Color',[.7 .7 .7],'Linewidth',2);
        else
            plot(t(Lt),h(Lt),'-','Color',[0 0 1],'Linewidth',2);
            plot(t(Lt),v(Lt),'-','Color',[0 .9 .9],'Linewidth',2);
        end
    end
    set(gca,'Ylim',[-3 3],'Xlim',[0.1 .33]);   
end

%%
% Section 2
% DToneloc showing spatial spread
% These file names should get pulled from DB
filenamestem = 'M051818';
stros ={};
stim_locations_xy = [];
for i = 2:15 % Hack
    suffix = num2str(i);
    suffix = [repmat('0',1,3-length(suffix)) suffix];
    stro = nex2stro(findfile([filenamestem,suffix,'.nex']));
    stim_locations_xy = [stim_locations_xy; stro.sum.exptParams.rf_x stro.sum.exptParams.rf_y];
    stros{length(stros)+1} = stro;
end
% Now culling the data I don't want for the figure
desired_locations = [8 -51; 27 -57; 47 -57];

% Testing 
%desired_locations = unique(stim_locations_xy,'rows');

L = zeros(size(stim_locations_xy,1),1);
for i = 1:size(stim_locations_xy)
    if ismember(stim_locations_xy(i,:),desired_locations,'rows')
        L(i) = 1;
    end
end
stim_locations_xy(~L,:) = [];
stros(~L) = [];
uniquexys = unique(stim_locations_xy,'rows');

% Doing the work
for i = 1:size(uniquexys,1)
    stro = [];
    % Pooling data from experiments with the same stimulus location
    for j = 1:size(stim_locations_xy,1)
        if stim_locations_xy(j,1) == uniquexys(i,1) & stim_locations_xy(j,2) == uniquexys(i,2)
            stro = strocat(stro,stros{j});
        end
    end
    % Setting up a bunch of crap that I might need
    Lstimpresent = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimpresent'));
    Lcorrect = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'correct'));
    Llaser = logical(stro.trial(:,strcmp(stro.sum.trialFields(1,:),'optstim')));
    Lseenresp = (Lcorrect&Lstimpresent) |(~Lcorrect&~Lstimpresent);
    lcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'lcc'));
    mcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'mcc'));
    scc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'scc'));
    Loog = logical(stro.trial(:,strcmp(stro.sum.trialFields(1,:),'oog')));
    
    approxgamutedgeLcc = max(lcc(~Loog));
   % nhit_laser = sum(Lcorrect(Lstimpresent&Llaser));
   % nhit_nolaser = sum(Lcorrect(Lstimpresent&~Llaser));
    nFA_laser = sum(~Lcorrect(~Lstimpresent&Llaser));
    nFA_nolaser = sum(~Lcorrect(~Lstimpresent&~Llaser));
   % nmiss_laser = sum(~Lcorrect(Lstimpresent&Llaser));
   % nmiss_nolaser = sum(~Lcorrect(Lstimpresent&~Llaser));
    nCR_laser = sum(Lcorrect(~Lstimpresent&Llaser));
    nCR_nolaser = sum(Lcorrect(~Lstimpresent&~Llaser));
    nbins = 10;
    lccbinedges = logspace(log10(min(lcc(lcc > 0))),log10(max(lcc(~Loog))),nbins+1);
    lccbincenters = 10.^(log10(lccbinedges(1:end-1))+diff(log10(lccbinedges))/2);
    figure; axes; hold on;
    for k = 0:1 % No laser/laser
        if k == 0
            L = ~Llaser&Lstimpresent;
        else
            L = Llaser&Lstimpresent;
        end
        x = zeros(nbins,2); % hits, total #trials
        if any(Loog)
            gamutedge = max(lcc(~Loog));
            lcc(Loog) = gamutedge;
        end
        rawresps = [discretize(lcc(L),lccbinedges) Lseenresp(L)];
        for j = 1:nbins
            Lbin = rawresps(:,1) == j;
            x(j,1) = sum(rawresps(Lbin,2));
            x(j,2) = sum(Lbin);
        end
        % z scores for d-prime calculation: d' = Z(p(hit)) - Z(p(FA))
        if k
            FAterm = norminv(nFA_laser./(nFA_laser+nCR_laser));
        else
            FAterm = norminv(nFA_nolaser./(nFA_nolaser+nCR_nolaser));
        end
        for j = 1:nbins
            h = plot(lccbincenters(j),min(norminv(x(j,1)./x(j,2))-FAterm,10),'ko','MarkerFaceColor','black','MarkerSize',max(x(j,2),.01));
            if k
                set(h,'Color','blue','MarkerFaceColor','blue');
            end
        end
    end
    set(gca,'Xscale','log','Xlim',[lccbinedges(1) lccbinedges(end)]);
    title(num2str(uniquexys(i,:)));
    ylabel('d''');
    xlabel('lcc');
end

% -----------------
% Finding the first contrast at which d' exceeds a threshold
data = [];
DPRIMETHRESH = 1;
for i = 1:size(uniquexys,1)
    stro = [];
    % Pooling data from experiments with the same stimulus location
    for j = 1:size(stim_locations_xy,1)
        if stim_locations_xy(j,1) == uniquexys(i,1) & stim_locations_xy(j,2) == uniquexys(i,2)
            stro = strocat(stro,stros{j});
        end
    end
    % Setting up a bunch of crap that I might need
    Lstimpresent = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimpresent'));
    Lcorrect = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'correct'));
    Llaser = logical(stro.trial(:,strcmp(stro.sum.trialFields(1,:),'optstim')));
    Lseenresp = (Lcorrect&Lstimpresent) |(~Lcorrect&~Lstimpresent);
    lcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'lcc'));
    mcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'mcc'));
    scc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'scc'));
    Loog = logical(stro.trial(:,strcmp(stro.sum.trialFields(1,:),'oog')));
    
    thresholds = []
    for k = 0:1 % No laser/laser
        if k == 0
            Llasercond = ~Llaser;
        else
            Llasercond = Llaser;
        end
        if any(Loog)
            gamutedge = max(lcc(~Loog));
            lcc(Loog) = gamutedge;
        end
        uniquecontrasts = unique(lcc(Llasercond));
        tmp = [];
        for m = 1:length(uniquecontrasts)
            % Below: contrast, nhits, ntotal
            tmp = [tmp; uniquecontrasts(m) sum(lcc == uniquecontrasts(m) & Llasercond & Lcorrect) sum(lcc == uniquecontrasts(m) & Llasercond)];
        end
        if tmp(1,1) ~= 0 % making sure lowest contrast is zero
            error('lowest contrast is not zero');
        end
        FAterm = norminv((tmp(1,3)-tmp(1,2))./tmp(1,3));
        HITterm = norminv(tmp(:,2)./tmp(:,3));
        HITterm(1) = FAterm; % hits don't make sense for zero contrast
        dprime = HITterm-FAterm;
        whichcontrast = find(dprime > DPRIMETHRESH,1,'first');
        if ~isempty(whichcontrast)
            thresholds(k+1) = uniquecontrasts(whichcontrast);
        else
            thresholds(k+1) = nan;
        end
        nhit(k+1) = sum(Lstimpresent & Llasercond & Lcorrect);
        nmiss(k+1) = sum(Lstimpresent & Llasercond & ~Lcorrect);
        nFA(k+1) = sum(~Lstimpresent & Llasercond & ~Lcorrect);
        nCR(k+1) = sum(~Lstimpresent & Llasercond & Lcorrect);
        
    end
    data = [data; uniquexys(i,:) thresholds nhit nmiss nFA nCR];
end

% Doing the plotting
figure; axes; hold on;
for i = 1:size(data,1)
    plot(data(i,1),data(i,2),'ko','MarkerSize',10)
end
plot(0,0,'k+','MarkerSize',15,'LineWidth',3);
set(gca,'Ylim',[-61 0],'Xlim',[-7.9208 69.4208]); % To match SC_stimcue fig

sens = [1./data(:,3)'; 1./data(:,4)']; % contrast sensitivity = 1/threshold
sens(isnan(sens)) = 0;
figure;
for i = 1:size(data,1)
    lasertrialstats = data(i,6:2:12);
    nolasertrialstats = data(i,5:2:11);
    n_lasertrials = sum(lasertrialstats);
    n_nolasertrials = sum(nolasertrialstats);
    subplot(2,3,i); hold on;
    props = [lasertrialstats./n_lasertrials; nolasertrialstats./n_nolasertrials];
   % h = bar(props');
   h = bar([lasertrialstats;nolasertrialstats]');
    set(h(1),'FaceColor','blue');
    set(h(2),'FaceColor',[.5 .5 .5])
   % sep = sqrt([(lasertrialstats./n_lasertrials).*(1-lasertrialstats./n_lasertrials)./n_lasertrials;...
   %             (nolasertrialstats./n_nolasertrials).*(1-nolasertrialstats./n_nolasertrials)./n_nolasertrials]);
   % plot([[1:4]-.15;[1:4]-.15],[props(1,:)-sep(1,:);props(1,:)+sep(1,:)],'b-','LineWidth',2,'Color',[0 0 .5])
   % plot([[1:4]+.15;[1:4]+.15],[props(2,:)-sep(2,:);props(2,:)+sep(2,:)],'k-','LineWidth',2)
   % [~,p] = equalproptest([lasertrialstats(1),nolasertrialstats(1)],[n_lasertrials/2 n_nolasertrials/2],.01)
    set(gca,'Xlim',[0 5],'Xtick',[],'Ylim',[0 .5]);
    subplot(2,3,i+3); hold on;
    h = bar(1,sens(1,i),'k'); set(h,'FaceColor',[.5 .5 .5]);
    h = bar(2,sens(2,i));    set(h,'FaceColor',[0 0 1]);
    set(gca,'Xlim',[0 16],'Xtick',[],'YScale','log','Ylim',[1 10]);
end
%%
% Section 2A
% Repeating above analysis but for a file with stimuli whose positions are
% randomly interleaved.

stro = nex2stro(findfile('M052518009.nex'));
stim_x = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_x'));
stim_y = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_y'));
stim_locations_xy = [stim_x stim_y];
uniquexys = unique(stim_locations_xy,'rows');
Lstimpresent = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimpresent'));
Lcorrect = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'correct'));
Llaser = logical(stro.trial(:,strcmp(stro.sum.trialFields(1,:),'optstim')));
Lseenresp = (Lcorrect&Lstimpresent) |(~Lcorrect&~Lstimpresent);
lcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'lcc'));
mcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'mcc'));
scc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'scc'));
Loog = logical(stro.trial(:,strcmp(stro.sum.trialFields(1,:),'oog')));
approxgamutedgeLcc = max(lcc(~Loog));

DPRIMETHRESH = 1;
for i = 1:size(uniquexys,1)
    Lxy = stim_x == uniquexys(i,1) & stim_y == uniquexys(i,2);
    figure; axes; hold on;
    for k = 0:1 % No laser/laser
        if k == 0
            Llasercond = ~Llaser;
        else
            Llasercond = Llaser;
        end
        if any(Loog)
            gamutedge = max(lcc(~Loog));
            lcc(Loog) = gamutedge;
        end
        uniquecontrasts = unique(lcc(Llasercond & Lxy));
        data = [];
        for m = 1:length(uniquecontrasts)
            ntrials = sum(lcc == uniquecontrasts(m) & Llasercond & Lxy);
            nhits = sum(lcc == uniquecontrasts(m) & Llasercond & Lxy & Lcorrect);
            data = [data; uniquecontrasts(m) nhits ntrials];
        end
        if data(1,1) ~= 0 % making sure lowest contrast is zero
            error('lowest contrast is not zero');
        end
        FAterm = norminv((data(1,3)-data(1,2))./data(1,3));
        HITterm = norminv(data(:,2)./data(:,3));
        HITterm(1) = FAterm; % hits don't make sense for zero contrast
        dprime = HITterm-FAterm;
        whichcontrast = find(dprime > DPRIMETHRESH,1,'first');
        if ~isempty(whichcontrast)
            threshold = uniquecontrasts(whichcontrast);
        else
            threshold = nan;
        end
        threshold
        pause
    end
end


nbins = 10;
for i = 1:size(uniquexys,1)
    Lxy = stim_x == uniquexys(i,1) & stim_y == uniquexys(i,2);

    figure; axes; hold on;
    for k = 0:1 % No laser/laser
        if k == 0
            L = ~Llaser&Lstimpresent&Lxy;
        else
            L = Llaser&Lstimpresent&Lxy;
        end
        lccbinedges = logspace(log10(min(lcc(lcc > 0)))-eps,log10(max(lcc(~Loog & Lxy))),nbins+1);
        lccbincenters = 10.^(log10(lccbinedges(1:end-1))+diff(log10(lccbinedges))/2);

        x = zeros(nbins,2);
        rawresps = [discretize(lcc(L),lccbinedges) Lseenresp(L)];
        for j = 1:nbins+1
            if j <= nbins
                Lbin = rawresps(:,1) == j;
            else
                Lbin = isnan(rawresps(:,1));
            end
            x(j,1) = sum(rawresps(Lbin,2));
            x(j,2) = sum(Lbin);
        end
        % Putting all the out of gamut responses in the highest contrast
        % bin. Ony do this once!
        %x(nbins,:) = x(nbins,:)+x(nbins+1,:);
        x(nbins+1,:) = []; % Chucking the OOG points
        for j = 1:nbins
            if x(j,2) > 0 % Don't plot anything if the bin is empty
                h = plot(lccbincenters(j),x(j,1)./x(j,2),'ko','MarkerSize',x(j,2));
                if k
                    set(h,'Color','blue');
                end
            end
        end
    end
    title(num2str(uniquexys(i,:)));
end

%%
% Section 3
% Example DToneloc data file
filename = 'M051818012.nex';
stro = nex2stro(findfile(filename,[nexfilepath,filesep,'Abhishek']));
Lstimpresent = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimpresent'));
stim_idx = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_idx'));
correct = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'correct'));
Llaser = logical(stro.trial(:,strcmp(stro.sum.trialFields(1,:),'optstim')));
Lseenresp = (correct&stimpresent) |(~correct&~stimpresent);
lcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'lcc'));
Loog = logical(stro.trial(:,strcmp(stro.sum.trialFields(1,:),'oog')));

% numbers of different types of corrects/incorrects
nhit_laser = sum(correct(stimpresent&Llaser));
nhit_nolaser = sum(correct(stimpresent&~Llaser));
nFA_laser = sum(~correct(~stimpresent&Llaser));
nFA_nolaser = sum(~correct(~stimpresent&~Llaser));
nmiss_laser = sum(~correct(stimpresent&Llaser));
nmiss_nolaser = sum(~correct(stimpresent&~Llaser));
nCR_laser = sum(correct(~stimpresent&Llaser));
nCR_nolaser = sum(correct(~stimpresent&~Llaser));

figure;
h = bar([nhit_laser nhit_nolaser; nmiss_laser nmiss_nolaser; nFA_laser nFA_nolaser; nCR_laser nCR_nolaser]);
set(h(1),'FaceColor','blue');
set(h(2),'FaceColor',[.5 .5 .5])
set(gca,'Xticklabel',{'Hits','Misses','FAs','CRs'});
ylabel('trial count');

% Staircase plotting
figure; axes; hold on;
approxgamutedgeLcc = max(lcc(~Loog));
lcc(Loog) = approxgamutedgeLcc;
plot(lcc(~Llaser& Lstimpresent),'ko-','MarkerFaceColor','black','Linewidth',2);
plot(lcc(Llaser& Lstimpresent),'bo-','MarkerFaceColor','blue','Linewidth',2);
set(gca,'Yscale','log','Ylim',[.07 .6],'TickDir','out');
xlabel('Trial index');
ylabel('Contrast');

%%
% Section 3.1
% Eye movement analysis for the data file in section 3
filename = 'M051818012.nex';

stro = nex2stro(findfile(filename,[nexfilepath,filesep,'Abhishek']));
Lstimpresent = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimpresent'));
Lcorrect = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'correct'));
Llaser = logical(stro.trial(:,strcmp(stro.sum.trialFields(1,:),'optstim')));
Lseenresp = (Lcorrect&Lstimpresent) |(~Lcorrect&~Lstimpresent);
lcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'lcc'));
Loog = logical(stro.trial(:,strcmp(stro.sum.trialFields(1,:),'oog')));
stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimon_t'));
stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimoff_t'));
laseron_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'laseron'));
samplerate = stro.sum.analog.storeRates{strcmp(stro.sum.analog.sigid,'AD11')};
ntrials = size(stro.trial,1);
h_ep = cellfun(@(x)x*4096/400,stro.ras(:,strcmp(stro.sum.rasterCells,'AD11')),'UniformOutput',0);
v_ep = cellfun(@(x)x*4096/400,stro.ras(:,strcmp(stro.sum.rasterCells,'AD12')),'UniformOutput',0);
t_analogstart = cell2mat(stro.ras(:,strcmp(stro.sum.rasterCells,'anlgStartTime')));

h_data = cell(1,2);
v_data = cell(1,2);
t_offsets = [-.5, .5];
for j = 1:ntrials
    t = (t_analogstart(j)-stimon_t(j))+(1./samplerate)*[0:length(h_ep{j})-1];
    L = t >= t_offsets(1) & t < t_offsets(2);
    if sum(L) < (t_offsets(2)-t_offsets(1))*samplerate % avoiding round off errors
        L(find(L,1,'last')+1) = 1;
    end
    h_data{Llaser(j)+1} = [h_data{Llaser(j)+1}; h_ep{j}(L)'];
    v_data{Llaser(j)+1} = [v_data{Llaser(j)+1}; v_ep{j}(L)'];
end

titlestrs = ['H';'V'];
colors = [.7 .7 .7; 0 1 1];
figure;
t_plot = linspace(t_offsets(1),t_offsets(2),size(h_data{1},2));
for i = 1:2 % H, V
    subplot(2,2,i); hold on;
    for laser_idx = 1:2 % one means no laser
        if i == 1
            x = h_data{laser_idx}';
        else
            x = v_data{laser_idx}';
        end
        plot(t_plot, x,'k-','Color',colors(laser_idx,:),'LineWidth',.5');
        handle(laser_idx) = plot(t_plot, mean(x'),'r-','LineWidth',2,'Color',colors(laser_idx,:)/2);
        title(titlestrs(i));
        ylabel('deg');
    end
    uistack(handle,'top');
end

% Now microsaccades
sacstats = getSacData(stro); close;
msactimebins = linspace(t_offsets(1),t_offsets(2),50);
binwidth = msactimebins(2)-msactimebins(1); % s
msactimehist = zeros(2,length(msactimebins));
for j = 1:ntrials
    tmptimes = sacstats.starttimes{j}-stimon_t(j);
    tmptimes(tmptimes<t_offsets(1) | tmptimes>t_offsets(2)) = [];
    [n,x] = hist(tmptimes,msactimebins);
    msactimehist(Llaser(j)+1,:) = msactimehist(Llaser(j)+1,:) + n;
end

subplot(2,2,3); hold on;
for laser_idx = 1:2 % one means no laser
    h = bar(msactimebins,msactimehist(laser_idx,:)/sum(Llaser == laser_idx-1)./binwidth);
    set(h,'FaceColor',colors(laser_idx,:),'EdgeColor',colors(laser_idx,:)/2);
end
set(gca,'Xlim',[t_offsets(1)-binwidth/2,t_offsets(2)+binwidth/2],'TickDir','out');
xlabel('Time (s)');
ylabel('Microsaccades/s');
set(gcf,'Renderer','painters');
%%
% Section 4
% Pseudo ROC

conn = database('Abhishek','','','Vendor','MySql','Server','128.95.153.12');
location_query = 'SELECT filename FROM DToneloc WHERE training = ''no''';
filenames = fetch(conn, location_query);
close(conn);
nfp = [nexfilepath,filesep,'Abhishek'];

data = [];

for fileidx = 1:length(filenames)
    disp([fileidx length(filenames)]);
    stro = nex2stro(findfile(filenames{fileidx}, nfp));
    Lstimpresent = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimpresent'));
    Lcorrect = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'correct'));
    Llaser = logical(stro.trial(:,strcmp(stro.sum.trialFields(1,:),'optstim')));
    
    lcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'lcc'));
    mcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'mcc'));
    scc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'scc'));
    oog = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'oog'));
    
    % numbers of different types of corrects/incorrects
    nhit_laser = sum(Lcorrect&Lstimpresent&Llaser);
    nhit_nolaser = sum(Lcorrect&Lstimpresent&~Llaser);
    nFA_laser = sum(~Lcorrect&~Lstimpresent&Llaser);
    nFA_nolaser = sum(~Lcorrect&~Lstimpresent&~Llaser);
    nmiss_laser = sum(~Lcorrect&Lstimpresent&Llaser);
    nmiss_nolaser = sum(~Lcorrect&Lstimpresent&~Llaser);
    nCR_laser = sum(Lcorrect&~Lstimpresent&Llaser);
    nCR_nolaser = sum(Lcorrect&~Lstimpresent&~Llaser);
    
    data = [data; nhit_laser nhit_nolaser nFA_laser nFA_nolaser nmiss_laser nmiss_nolaser nCR_laser nCR_nolaser];
end

% --------------------------------------------
% Binning CRs to make an "ROC" type function.
figure; axes; hold on;
nblocks_stimabsent_laser = data(:,7)+data(:,3); % number of blocks
nblocks_stimabsent_nolaser = data(:,8)+data(:,4);
p_fa_laser = data(:,3)./nblocks_stimabsent_laser;
p_fa_nolaser = data(:,4)./nblocks_stimabsent_nolaser;
binedges = [0:.12:1]; % binning on y (nolaser)
for i = 1:length(binedges)-1
    L = p_fa_nolaser > binedges(i) & p_fa_nolaser < binedges(i+1);
    if any(L)
        mn_fa_nolaser = mean(p_fa_nolaser(L));
        mn_fa_laser = mean(p_fa_laser(L));
        % Variance of mean across blocks is the sum of variances across blocks
        % divided by number of blocks (squared). Variance within each block is
        % computed under binomial assumptions: (p*(1-p))/n.
        % So I'm computing the SE of the mean assuming that the things that
        % went into the mean are binomial proportions.
        var_fa_nolaser = sum((p_fa_nolaser(L).*(1-p_fa_nolaser(L)))./nblocks_stimabsent_nolaser(L))/sum(L)^2;
        var_fa_laser = sum((p_fa_laser(L).*(1-p_fa_laser(L)))./nblocks_stimabsent_laser(L))/sum(L)^2;
       % plot(mn_fa_nolaser,mn_fa_laser,'ko','MarkerFaceColor','black');
       % plot(mn_fa_nolaser+(sqrt(var_fa_nolaser).*[-1 1]),mn_fa_laser.*[1 1],'k-');
       % plot(mn_fa_nolaser.*[1 1],mn_fa_laser+(sqrt(var_fa_laser).*[-1 1]),'k-');
         plot(mn_fa_laser, mn_fa_nolaser,'ko','MarkerFaceColor','black');
         plot(mn_fa_laser.*[1 1],mn_fa_nolaser+(sqrt(var_fa_nolaser).*[-1 1]),'k-');
         plot(mn_fa_laser+(sqrt(var_fa_laser).*[-1 1]),mn_fa_nolaser.*[1 1],'k-');
    end
end

plot(p_fa_laser,p_fa_nolaser,'k.','MarkerEdgeColor',[.75 .75 .75])
axis square;
xlabel('p(FA) laser');
ylabel('p(FA) no laser');
plot([0 1],[0 1],'k-');

% Fitting ROCs. 
% ROC data = [nFAs(laser) nFAs(no laser) nstimabsent(laser) nstimabsent(nolaser)]
% substitute FA(no laser) for hits on the y-axis
% substitute FA(laser) for FA on the x-axis
ROC_data = [data(:,3) data(:,4) data(:,3)+data(:,7) data(:,4)+data(:,8)];
[muhat, sigmahat, initialguesses, stats] = FitROC(ROC_data); % here's the call to the fitting routine
tmp = linspace(-4,4,100);
plot(normcdf(tmp,muhat,sigmahat),normcdf(tmp,0,1),'b-');

200*(sum(normcdf(x,-muhat,1).*(y./sum(y)))-.5)
200*(sum(normcdf(x,-stats.muCI(1),1).*(y./sum(y)))-.5)
200*(sum(normcdf(x,-stats.muCI(2),1).*(y./sum(y)))-.5)


%%
% Section 5
% Schematic of the pseudo ROC curve

x = linspace(-4,6,100);
mu_noise1 = 0; % no laser
sigma_noise1 = 1;
mu_noise2 = -.4; % laser
sigma_noise2 = .8;
mu_signal = 2; % signal
sigma_signal = 1;

signalpdf = normpdf(x,mu_signal, sigma_signal)./sum(normpdf(x,mu_signal, sigma_signal));
noise1pdf = normpdf(x,mu_noise1, sigma_noise1)./sum(normpdf(x,mu_noise1, sigma_noise1)); % no laser
noise2pdf = normpdf(x,mu_noise2, sigma_noise2)./sum(normpdf(x,mu_noise2, sigma_noise2)); % laser

cyan = [0 174 239]/255;
gray = [.4 .4 .4];
figure; axes; hold on;
plot(x,signalpdf,'k--','Color',[.5 .5 .5],'LineWidth',2);
plot(x,noise1pdf,'k','LineWidth',3);
plot(x,noise2pdf,'c','LineWidth',3,'Color',cyan);
plot([-2 -2],[0 0.05],'k:','LineWidth',1,'Color',gray*.9);
plot([-1 -1],[0 0.05],'k:','LineWidth',1,'Color',gray*.7);
plot([0 0],[0 0.05],'k:','LineWidth',1,'Color',gray*.5);
plot([1 1],[0 0.05],'k:','LineWidth',1,'Color',gray*.3);
%plot([2 2],[0 0.05],'k:','LineWidth',1,'Color',gray*.3);
set(gca,'TickDir','out');

figure; axes; hold on;
plot([0 1],[0 1],'k--');
plot(1-cumsum(noise2pdf),1-cumsum(noise1pdf),'c-','LineWidth',3,'Color',[0 174 239]/255);
plot(1-sum(noise2pdf(x<-2)),1-sum(noise1pdf(x<-2)),'*','MarkerSize',11,'Color',gray*.9);
plot(1-sum(noise2pdf(x<-1)),1-sum(noise1pdf(x<-1)),'*','MarkerSize',11,'Color',gray*.7);
plot(1-sum(noise2pdf(x<0)),1-sum(noise1pdf(x<0)),'*','MarkerSize',11,'Color',gray*.5);
plot(1-sum(noise2pdf(x<1)),1-sum(noise1pdf(x<1)),'*','MarkerSize',11,'Color',gray*.3);
%plot(1-sum(noise2pdf(x<2)),1-sum(noise1pdf(x<2)),'*','MarkerSize',11,'Color',gray*.3);

axis square;
set(gca,'TickDir','out','Xlim',[0 1],'Ylim',[0 1],'Xtick',[0 .5 1],'Ytick',[0 .5 1]);
xlabel('p(FA) laser');
ylabel('p(FA) no laser');

%%
% Section 6
% Spatial uncertainty simulation.


npix = 5;
signal_template = zeros(npix);
signal_template(ceil(npix/2),ceil(npix/2)) = 1;
noise_template = ones(npix);
center_surround_filter = signal_template-(noise_template./sum(noise_template(:)));
signalamplitudes = [0 logspace(-1.5,1.5,10)];
niter = 10000;
optomask = ones(npix);
optomask(ceil(npix/2)+1:end,:) = 0;
optomask = signal_template;
optomask(1,1) = 1;

data = [];
for i = 1:length(signalamplitudes)
    for j = 1:niter
        noise = 10+normrnd(0,2.5)*noise_template+normrnd(0,.7,size(noise_template));
        sig_plus_noise = signalamplitudes(i)*signal_template+noise;
        
        dv1 = center_surround_filter(:)'*sig_plus_noise(:);
        dv2 = max(sig_plus_noise(:));
        dv3 = sum(sum(center_surround_filter.*sig_plus_noise.*optomask));
        dv4 = max(max(sig_plus_noise.*optomask));
        data = [data; signalamplitudes(i) dv1 dv2 dv3 dv4];
    end
end

% Order of columns
% 1) signal amplitude
% 2) matched filter without optomask
% 3) maximum without optomask
% 4) matched filter with optomask
% 5) maximum with optomask

Lblank = data(:,1) == 0;
uniquestim = unique(data(:,1));
uniquestim(uniquestim == 0) = [];
titles = {'matched filter, no mask','maximum, no mask','matched filter, with mask','maximum with mask'};
metadata = [];
for i = 1:length(uniquestim)
  %  figure;
    for j = 1:4
        subplot(2,2,j); hold on;
        L = data(:,1) == uniquestim(i);
        hist([data(Lblank,j+1) data(L,j+1)]);
        r = roc(data(Lblank,j+1), data(L,j+1));
        title([titles{j},' roc = ',num2str(r)]);
        metadata(i,j) = r;
    end
end

figure; 
subplot(2,2,1); hold on;
plot(uniquestim,metadata(:,1),'ko-','LineWidth',2,'MarkerFaceColor','black');
plot(uniquestim,metadata(:,3),'co-','LineWidth',2,'MarkerFaceColor','cyan');
set(gca,'Xscale','log','Xtick',[],'Ytick',[.5 .75 1],'Ylim',[.5 1],'Tickdir','out')
subplot(2,2,2); hold on;
plot(uniquestim,metadata(:,2),'ko-','LineWidth',2,'MarkerFaceColor','black');
plot(uniquestim,metadata(:,4),'co-','LineWidth',2,'MarkerFaceColor','cyan');
set(gca,'Xscale','log','Xtick',[],'Ytick',[.5 .75 1],'Ylim',[.5 1],'Tickdir','out')


%%
% Section 7
% Achromatic vs S-cone detection

stro = strocat(nex2stro(findfile('M032518019.nex')), nex2stro(findfile('M032518020.nex')));

stim_idx = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_idx'));
stro.trial(stim_idx == 1,:) = []; % Getting rid of L-M trials which the monkey couldn't see
stim_idx(stim_idx == 1,:) = []; % Getting rid of L-M trials which the monkey couldn't see
Lstimpresent = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimpresent'));
Lcorrect = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'correct'));
Llaser = logical(stro.trial(:,strcmp(stro.sum.trialFields(1,:),'optstim')));
lcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'lcc'));
mcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'mcc'));
scc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'scc'));
Loog = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'oog'));
for i = [0 2] % skipping L-M
    % Didn't end up using this stuff since S-cone went out of gamut on
    % no-laser trials
     L = stim_idx == i;
%     stimuli_laser = [lcc(L&Lstimpresent&Llaser) mcc(L&Lstimpresent&Llaser) scc(L&Lstimpresent&Llaser)];
%     oog = Loog(L&Lstimpresent&Llaser);
%     laser_contrasts = sqrt(sum(stimuli_laser.^2,2));
%     approxgamutedge = max([max(laser_contrasts(~oog)), min(laser_contrasts)]);
%     laser_contrasts(laser_contrasts > approxgamutedge) = approxgamutedge;
     stimuli_nolaser = [lcc(L&Lstimpresent&~Llaser) mcc(L&Lstimpresent&~Llaser) scc(L&Lstimpresent&~Llaser)];
%     nolaser_contrasts = sqrt(sum(stimuli_nolaser.^2,2));
%     nolaser_contrasts(nolaser_contrasts > approxgamutedge) = approxgamutedge;
    colordir = stimuli_nolaser(1,:)./norm(stimuli_nolaser(1,:));

    % numbers of different types of corrects/incorrects
    % Pooling FAs and CRs across color directions since they are just two
    % identical sets of interleaved blank trials
    nhit_laser = sum(L&Lcorrect&Lstimpresent&Llaser);
    nhit_nolaser = sum(L&Lcorrect&Lstimpresent&~Llaser);
    nFA_laser = sum(~Lcorrect&~Lstimpresent&Llaser);
    nFA_nolaser = sum(~Lcorrect&~Lstimpresent&~Llaser);
    nmiss_laser = sum(L&~Lcorrect&Lstimpresent&Llaser);
    nmiss_nolaser = sum(L&~Lcorrect&Lstimpresent&~Llaser);
    nCR_laser = sum(Lcorrect&~Lstimpresent&Llaser);
    nCR_nolaser = sum(Lcorrect&~Lstimpresent&~Llaser);
    
    n_stim_laser = nhit_laser+nmiss_laser;
    n_stim_nolaser = nhit_nolaser+nmiss_nolaser;
    n_nostim_laser = sum(~Lstimpresent&Llaser);
    n_nostim_nolaser = sum(~Lstimpresent&~Llaser);
    
    lasertrialstats = [nhit_laser nmiss_laser nFA_laser nCR_laser];
    n_lasertrials = [n_stim_laser n_stim_laser n_nostim_laser n_nostim_laser];
    n_lasertrials = repmat(n_stim_laser+n_nostim_laser,1,4);
    nolasertrialstats = [nhit_nolaser nmiss_nolaser nFA_nolaser nCR_nolaser];
    n_nolasertrials = [n_stim_nolaser n_stim_nolaser n_nostim_nolaser n_nostim_nolaser];
    n_nolasertrials = repmat(n_stim_nolaser+n_nostim_nolaser,1,4);
    props = [lasertrialstats./n_lasertrials; nolasertrialstats./n_nolasertrials];

    figure; axes; hold on;
    h = bar([lasertrialstats; nolasertrialstats]');  
    set(h(1),'FaceColor',[40 170 226]/255);
    set(h(2),'FaceColor',[.5 .5 .5]);
    
    %sep = sqrt([(lasertrialstats./n_lasertrials).*(1-lasertrialstats./n_lasertrials)./n_lasertrials;...
    %            (nolasertrialstats./n_nolasertrials).*(1-nolasertrialstats./n_nolasertrials)./n_nolasertrials]);
    %plot([[1:4]-.15;[1:4]-.15],[props(1,:)-sep(1,:);props(1,:)+sep(1,:)],'b-','LineWidth',2,'Color',[0 0 .5])
    %plot([[1:4]+.15;[1:4]+.15],[props(2,:)-sep(2,:);props(2,:)+sep(2,:)],'k-','LineWidth',2)
        
    set(gca,'Xticklabel',{'Hits','Misses','FAs','CRs'},'Xlim',[0 5],'Xtick',1:4,'Ylim',[0 25],'TickDir','out');
    ylabel('Trial proportion');
end

% Plotting some Gabors to go with it
isostim = [.3 .3 .3; 0 0 .6; 0 0 0];

M = [ 0.0608    0.1219    0.0175
    0.0220    0.1266    0.0257
    0.0019    0.0095    0.0976];
bkgndrgb = [.7 .7 .7];
bkgndlms = M*bkgndrgb';

edgergb = [0 0 0];
theta = 0;
lambda = 2;
sigma = .5;
gamma = 1;
phi = 0;
xoff = 0;
yoff = 0;
etheta = 0;
edisp = 0;
gausslim = .999;
pixperdeg = 40;

figure;
for i = 1:size(isostim,1)
    gaborrgb = inv(M)*(bkgndlms.*(1+isostim(i,:))') - bkgndrgb;
    im = DrawGaborEdge(bkgndrgb, gaborrgb, edgergb, theta, lambda, sigma, gamma, phi, xoff, yoff, etheta, edisp, gausslim, pixperdeg);
    subplot(ceil(sqrt(size(isostim,1))),ceil(sqrt(size(isostim,1))),i)
    image(im);
    set(gca,'Visible','off')
    axis image
end
set(gcf,'Renderer','painters');

%%
% Section 8
% SC_stimcue slide for Sackler meeting
% Potentials M042518001
% M042618002;
% M042718005 <-- good one


WINDOWDUR = 0.05;  % How long wait after SACINIT to plot trajectories (non catch trials)
RTDUR = 0.4; % Max RT

stro = nex2stro(findfile('M042718005.nex',[nexfilepath,filesep,'Abhishek']));
targ_x = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targ_x'));
targ_y = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targ_y'));
uniquetargxy = unique([targ_x targ_y],'rows');
ntrials = size(stro.trial,1);
Lcatchtrials = targ_x == 0 & targ_y == 0;
laseron_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimon_t'));
laseroff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimoff_t'));
fpoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fpoff_t'));
fix_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fpacq_t'));
targon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targon_t'));
targoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targoff_t'));
sacinit_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'saccinit_t'));
sacmade_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'saccmade_t'));

Llaser = ~isnan(laseron_t);
samplerate = stro.sum.analog.storeRates{1};
figure; axes; hold on;
sactimes = []; % start, stop relative to fpoff
for i = 1:size(stro.trial,1)
    x = stro.ras{i,2}*4096/400;
    y = stro.ras{i,3}*4096/400;
    t = stro.ras{i,5}+[0:1:length(x)-1]/samplerate;
    if ~Lcatchtrials(i)
        if sacinit_t(i) < fpoff_t(i) + RTDUR
            Lt = t>fpoff_t(i) & t < sacinit_t(i) + WINDOWDUR;
            sactimes = [sactimes; sacinit_t(i)-fpoff_t(i)];
            h = plot(x(Lt),y(Lt),'k','LineWidth',2);
            if Llaser(i)
                set(h,'color',[102 204 255]/255);
            else
                set(h,'color',[.7 .7 .7]);
            end
        end
    end
end
axis square
set(gca,'Xlim',[-10 10],'Ylim',[-10 10]);

% Laser catch trials
if any(Lcatchtrials&Llaser)
    figure; axes; hold on;
    set(gca,'GridLineStyle','--');
    for i = 1:size(stro.trial,1)
        x = stro.ras{i,2}*4096/400;
        y = stro.ras{i,3}*4096/400;
        t = stro.ras{i,5}+[0:1:length(x)-1]/samplerate;
        if Lcatchtrials(i) && Llaser(i)
            Lt1 = t>fpoff_t(i) & t<fpoff_t(i)+min(sactimes);
            Lt2 = t>fpoff_t(i)+max(sactimes) & t<fpoff_t(i)+max(sactimes)+WINDOWDUR;
            mn_x = mean(x(Lt2))-mean(x(Lt1));
            mn_y = mean(y(Lt2))-mean(y(Lt1));
            if abs(mn_x) < stro.sum.exptParams.eyewin_x/20 && abs(mn_x) < stro.sum.exptParams.eyewin_y/20
                mn_x = 0;
                mn_y = 0;
            end
            plot([0 mn_x],[0 mn_y],'k-','Color',[.7 .7 .7]);
            h = plot(mn_x,mn_y,'ko','MarkerSize',10,'Linewidth',1);
            set(h,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[.7 .7 .7]);
        end
    end
    plot(cos(linspace(0,2*pi,100))/2,sin(linspace(0,2*pi,100))/2,'r-','Linewidth',2)
    axis square
    set(gca,'Xlim',[-10 10],'Ylim',[-10 10],'Xtick',[-5 0 5],'Ytick',[-5 0 5]);
    grid on;
end

%%
% Section 9 Rasters
% (Taken from Abhishek's "Popfiguresforeditor.m".)

filename = {'A012319005.nex';'A020119004.nex'};
binwidth = .005;
bins = -0.2:binwidth:0.5;
figprefs(1);
for jj = 1:length(filename)
    stro = nex2stro(findfile(char(filename(jj,:))));
    stimpresentidx  = strcmp(stro.sum.trialFields(1,:),'stimpresent');
    Lcc = strcmp(stro.sum.trialFields(1,:),'lcc');
    Mcc = strcmp(stro.sum.trialFields(1,:),'mcc');
    Scc = strcmp(stro.sum.trialFields(1,:),'scc');
    tf = strcmp(stro.sum.trialFields(1,:),'tf');
    oog = strcmp(stro.sum.trialFields(1,:),'oog');
    optstim = strcmp(stro.sum.trialFields(1,:),'optstim');
    correct = strcmp(stro.sum.trialFields(1,:),'correct');
    laseron = strcmp(stro.sum.trialFields(1,:),'laseron_t');
    laseroff = strcmp(stro.sum.trialFields(1,:),'laseroff_t');
    analogstrtimeind = strcmp(stro.sum.rasterCells(1,:),'anlgStartTime');
    spikeind = strcmp(stro.sum.rasterCells(1,:),'sig001a');
    lasertraceind = strcmp(stro.sum.rasterCells(1,:),'AD13');
    stimpresent = logical(stro.trial(:,stimpresentidx));
    lasertrials = logical(stro.trial(:,optstim));
    samplerate = stro.sum.analog.storeRates{3}; % sample rate at which laser analog pulses are stored in file
    
    % Stimulus absent trials
    idxs = find(~stimpresent);
    counter = 0;
    PSTHlaser = zeros(1,length(bins));
    PSTHnolaser = zeros(1,length(bins));
    for ii = 1:numel(idxs) % looping over stimulus absent trial
        ind = idxs(ii);
        analogstartime = stro.ras{ind,analogstrtimeind};
        spiketimes = stro.ras{ind,spikeind};
        if lasertrials(ind)
            % Pulling out the laseron and off timings from analog traces
            t = stro.ras{ind,analogstrtimeind}+[0:1:length(stro.ras{ind,lasertraceind})-1]/samplerate;
            laserontime = t(find(stro.ras{ind,lasertraceind}>0.1,1));
            laserofftime = t(find(stro.ras{ind,lasertraceind}>0.1,1,'last'));
            subplot(2,4,2*jj-1); hold on;
            Ltime = spiketimes>(laserontime+bins(1)-binwidth) & spiketimes<(laserontime+bins(end)+binwidth);
            
            tmpspikes = spiketimes(Ltime)-laserontime;
            nspikestot = length(tmpspikes);
            plot([tmpspikes tmpspikes]',[zeros(nspikestot,1) .7*ones(nspikestot,1)]'+counter,'y-','linewidth',1);
            
            PSTHlaser = PSTHlaser + hist(tmpspikes, bins);

            counter = counter + 1;
        end
    end
    PSTHlaser = PSTHlaser/counter/binwidth;
    set(gca,'Xlim',[-0.15 0.45],'Tickdir','out','XTick',[-0.15:0.15:0.45]); xlabel('time(s)'); ylabel('Number of trials'); axis square; 
    
    subplot(2,4,2*jj); plot(bins,PSTHlaser,'color',[0 200/255 1],'Linewidth',3); 
    set(gca,'Xlim',[-0.15 0.45],'Tickdir','out','XTick',[-0.15:0.15:0.45],'YTick',[0:100:1000]);
    if max(PSTHlaser)<20
        set(gca,'Ylim',[0 20],'YTick',[0:5:20]);
    end
    xlabel('time(s)'); ylabel('Response (ips)'); axis square;     
end

