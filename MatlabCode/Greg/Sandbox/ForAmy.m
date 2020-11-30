% A potential background

x = linspace(0,10,500);
f = logspace(-1,.6,500);
a = linspace(1,.3,500);
a = linspace(1,1,500);
y = a.*cos(x.*f);
%plot(x,y);

imagesc(y'*y);
colormap(gray);
axis image; axis xy;
set(gca,'XTick',[],'YTick',[]);

%%
% Another potential background (a radial grating)

freq = 18;
DISKRAD = 8000;
DISKPOS1 = [ 781  500]
DISKPOS2 = [ 400  200]
DISKCOL = .6;
BKAMP = 1;


figure; axes;
set(gcf,'Position',[270 44 1171 740]);

set(gca,'Units','pixels');
pos = get(gca,'OuterPosition')

[x,y] = meshgrid([1:pos(3)],[1:pos(4)]);
centerpoint = [pos(3) pos(4)]*2/3;

theta = atan2(y-centerpoint(2),x-centerpoint(1));
colormap(gray(255))
im = BKAMP*cos(freq*theta)/2+.5;

disks = ((x-DISKPOS1(1)).^2+(y-DISKPOS1(2)).^2 < DISKRAD |...
        (x-DISKPOS2(1)).^2+(y-DISKPOS2(2)).^2 < DISKRAD);
im = im.*(1-disks);
    
image(255*max(cat(3,DISKCOL*disks,im),[],3))
axis xy;

%%
% Analysis of gratings data
PRINT = 0;
GT = nex2stro;
orients = GT.trial(:,strcmp(GT.sum.trialFields(1,:),'orient'));
sfs = GT.trial(:,strcmp(GT.sum.trialFields(1,:),'sf'));
tfs = GT.trial(:,strcmp(GT.sum.trialFields(1,:),'tf'));
unq_sfs = unique(sfs');
unq_tfs = unique(tfs');
if (length(unq_sfs) > length(unq_tfs))
    fs = sfs;
    unq_fs = unq_sfs;
else
    fs = tfs;
    unq_fs = unq_tfs;
end

unq_orients = unique(orients');
stimon_t = GT.trial(:,strcmp(GT.sum.trialFields(1,:),'stim_on'));
stimoff_t= GT.trial(:,strcmp(GT.sum.trialFields(1,:),'stim_off'));
spikeidxs = strncmp(GT.sum.rasterCells,'sig',3);
for spikeidx = find(spikeidxs)
    spikerates = [];
    baselines = [];  baseline_t = 0.5;
    for i = 1:size(GT.trial,1)
        spiketimes = GT.ras{i,spikeidx};
        nspikes = sum(spiketimes > stimon_t(i) & spiketimes < stimoff_t(i));
        spikerates = [spikerates; nspikes./(stimoff_t(i)-stimon_t(i))];
        nspikes = sum(spiketimes > stimon_t(i)-baseline_t & spiketimes < stimon_t(i));
        baselines = [baselines; nspikes./baseline_t];
    end
    
    frmap = zeros(length(unq_orients),length(unq_fs));
    semap = zeros(length(unq_orients),length(unq_fs));
    for j = 1:length(unq_fs)
        for k = 1:length(unq_orients)
            L = fs == unq_fs(j) & orients == unq_orients(k);
            frmap(k,j) = mean(spikerates(L))
            semap(k,j) = std(spikerates(L))./sqrt(sum(L));
        end
    end
    figure; 
    subplot(3,1,1); hold on;
    plot(unq_orients*180/pi,frmap,'Linewidth',2);
    plot(unq_orients*180/pi,frmap+semap,':');
    plot(unq_orients*180/pi,frmap-semap,':');
    legend(num2str(unq_fs'));
    xlabel('orientation (deg)');
    ylabel('firing rate (sp/s)');
    title(GT.sum.rasterCells{spikeidx});
    set(gca,'Xlim',[0 315]);
    subplot(3,1,2);
    polar([unq_orients unq_orients(1)],[mean(frmap') mean(frmap(1,:))]);
    subplot(3,1,3);
    semilogx(unq_fs, mean(frmap),'k.-');
    set(gca,'Xlim',[min(fs) max(fs)])
    title(GT.sum.rasterCells{spikeidx});
    if (PRINT)
        eval(['print -dpsc junk',num2str(spikeidx)]);
    end
end


%%
% Population orientation tuning

filenames = {'A080313009','A080413010','A080513008','A080613008','A082213009','A082313006','A082413007','A082713007','A090913007','A091113009','A091713005','A092013005','A092313005'};

metadata = [];
for filecounter = 1:length(filenames)
    GT = nex2stro(fullfile(nexfilepath,'Amy','Apollo',[filenames{filecounter},'.nex']));
    orients = GT.trial(:,strcmp(GT.sum.trialFields(1,:),'orient'));
    sfs = GT.trial(:,strcmp(GT.sum.trialFields(1,:),'sf'));
    unq_sfs = unique(sfs');
    unq_orients = unique(orients');
    stimon_t = GT.trial(:,strcmp(GT.sum.trialFields(1,:),'stim_on'));
    stimoff_t= GT.trial(:,strcmp(GT.sum.trialFields(1,:),'stim_off'));
    spikeidxs = strncmp(GT.sum.rasterCells,'sig',3);
    hs = [];
    for spikeidx = find(spikeidxs)
        spikerates = [];
        baselines = [];  baseline_t = 0.5;
        for i = 1:size(GT.trial,1)
            spiketimes = GT.ras{i,spikeidx};
            nspikes = sum(spiketimes > stimon_t(i) & spiketimes < stimoff_t(i));
            spikerates = [spikerates; nspikes./(stimoff_t(i)-stimon_t(i))];
            nspikes = sum(spiketimes > stimon_t(i)-baseline_t & spiketimes < stimon_t(i));
            baselines = [baselines; nspikes./baseline_t];
        end
        
        frmap = zeros(length(unq_orients),length(unq_sfs));
        semap = zeros(length(unq_orients),length(unq_sfs));
        for j = 1:length(unq_sfs)
            for k = 1:length(unq_orients)
                L = sfs == unq_sfs(j) & orients == unq_orients(k);
                frmap(k,j) = mean(spikerates(L));
                semap(k,j) = std(spikerates(L))./sqrt(sum(L));
            end
        end
        %subplot(1,sum(spikeidxs),spikeidx); hold on;
        figure; axes; hold on;
        plot(unq_orients*180/pi,frmap,'Linewidth',2);
        plot(unq_orients*180/pi,frmap+semap,':');
        plot(unq_orients*180/pi,frmap-semap,':');
        legend(num2str(unq_sfs'));
        xlabel('orientation (deg)');
        ylabel('firing rate (sp/s)');
        title([filenames{filecounter},': ',GT.sum.rasterCells{spikeidx}]);
        set(gca,'Xlim',[0 315]);
    end
    
    % Fitting a *real* period spline
    fr = frmap(:,logical(sum(frmap == max(frmap(:)))))';
    [x,y] = pol2cart([unq_orients unq_orients(1)],[fr fr(1)]);
    pp = cscvn([x;y]);
    tmph = figure; % clunky way to evaluating the spline at a bunch of points
    fnplt(pp);
    yy = get(get(gca,'Children'),'YData');
    xx = get(get(gca,'Children'),'XData');
    close(tmph);
    [th,r] = cart2pol(xx,yy);
    idx = find(r == max(r),1);
    
    %metadata = [metadata; atan2(yy(idx),xx(idx))];
    metadata = [metadata; unq_orients(logical(sum(frmap == max(frmap(:)),2)))];
end
figure;
hist(metadata*180/pi,10)
%hist(mod(metadata*180/pi,180),10)

%%
% What magnitude of artifact do we expect from assuming a fixed
% proportion between cm on the computer screen and degrees of visual angle?
% Can thins explain Why Amy needed to adjust the monitor for each cell to
% get the size-tuning curves (on the gray background) to be aligned?

% Bottom line - this effect is very small. I don't think that this could be
% the entire explanation for the discrepansies between the size-tuning
% curves measured on the greay background.

% Shortest distance from eye to screen
VIEWINGDIST = 50; % cm 
% How far the eye has to rotate to obtain an eccentric fixation point
ECCPOS = 20; % deg

SIZEATFP = tan(1*pi/180)*VIEWINGDIST; % distance on screen (in cm) subtended by 1 deg straight ahead
SIZEATECCPOS = tan((ECCPOS+1)*pi/180)*VIEWINGDIST - tan((ECCPOS)*pi/180)*VIEWINGDIST; % distance on screen (in cm) subtended by 1 deg off to the side

% Difference in screen distance travelled by a 1 deg rotation
SIZEATECCPOS-SIZEATFP  % in cm

% Now for a fixed distance on the screen, what is the difference in angle?
DVAATFP = atan2(1,VIEWINGDIST)*180/pi
% Where (in cm) on the monitor is the stimulus that is ECCPOS dva from the
% fp?
ECCPOSCM = tan(ECCPOS*pi/180)*VIEWINGDIST;
DVAATECCPOS = (atan2(ECCPOSCM+1, VIEWINGDIST)-atan2(ECCPOSCM, VIEWINGDIST))*180/pi;
DVAATFP-DVAATECCPOS
% This is the difference (in DVA) betwen the angular subtense of a 1 cm
% stimulus at the FP and at the eccentric location

%%
% Eye movement analysis testing ground
stro = nex2stro(findfile('A081713006.nex'));
offset = [0 .1]; % in sec
hepidx = strcmp(stro.sum.rasterCells(1,:),'AD11');
vepidx = strcmp(stro.sum.rasterCells(1,:),'AD12');
starttimeidx = strcmp(stro.sum.rasterCells(1,:),'anlgStartTime');
stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_on'));
stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_off'));
fpx = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fp_x'));
fpy = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fp_y'));
uniquefppos = unique([fpx fpy],'rows');
bkgnd = unique(stro.trial(:,strcmp(stro.sum.trialFields(1,:),'bkgnd')));
samplerate = stro.sum.analog.storeRates{1};

ntrials = size(stro.ras,1);
% figure; axes; hold on;
data = nan*ones(ntrials,ceil(max(stimoff_t+offset(2)-stimon_t-offset(1))*samplerate),2);  % Eye position snippet
for i = 1:ntrials
    h = stro.ras{i,hepidx}.*4096/400;
    v = stro.ras{i,vepidx}.*4096/400;
    
    % 4096 A/D levels = 10 V
    % 1 degree = 40 A/D levels (according to REX)
    % (4096 levels/10 V) * (1 degree/40 levels) = 4096/400.
    % The key is that both REX and PLEXON use 12 bit A/D boards configured
    % for +/- 5 V.
    
   	e1t = stro.ras{i,starttimeidx};
    neyesamp = size(h,1);
    x = linspace(e1t,neyesamp/samplerate+e1t,neyesamp);
    L = x >=stimon_t(i)+offset(1) & x <=stimoff_t(i)+offset(2);
 %   plot(x(L)-stimon_t(i),h(L),'r-');
 %   plot(x(L)-stimon_t(i),v(L),'g-');
    data(i,1:sum(L),1) = h(L);
    data(i,1:sum(L),2) = v(L);
end

figure;
x = linspace(offset(1), offset(1)+2*(size(data,2)-1), size(data,2));
for i = 1:2
    for j = 1:2
        subplot(2,2,j+2*(i-1)); hold on;
        L = fpx == uniquefppos(i,1) & fpy == uniquefppos(i,2);
        plot(x,squeeze(data(L,:,j))','r-');
        plot(x,mean(squeeze(data(L,:,j))),'k-','LineWidth',2);
        if (j == 1)
            title(['Horz ',num2str(bkgnd)]);
        else
            title(['Vert ',num2str(bkgnd)]);
        end
        set(gca,'Xlim',[x(1) x(end)]);
        % Lnotnan = ~any(isnan(squeeze(data(L,:,1))));
        % [coeff, score, latent] = princomp(squeeze(data(L,Lnotnan,1)));
        % plot(coeff(:,1));
    end
end

%%
% Collecting a few files and seeing if I can eke out any difference in eye
% position between the corridor and the gray backgrounds.

filelist = {'F092412001.nex','F092412002.nex','F092412003.nex','F092412004.nex','F092412005.nex'};  % Size tuning shift, standard geometry
% filelist = {'F111412001.nex','F111412002.nex','F111412003.nex','F111412004.nex','F111412005.nex','F111412006.nex'}; % Size tuning shift, "eye movement control"
%filelist = {'A080313003.nex','A080313004.nex','A080313006.nex','A080313007.nex','A080313008.nex'};  % Size tuning shift, standard geometry
%filelist = {'A081113003.nex','A081113004.nex','A081113005.nex','A081113006.nex','A081113007.nex','A081113008.nex','A081113009','A081113010','A081113011'};  % Size tuning shift, "eye movement control"
%filelist = {'A082213001.nex','A082213002.nex','A082213003.nex','A082213004.nex','A082213005.nex','A082213006.nex','A082213007.nex'}; % standard geometry - shift in EP but not in tuning curve?
%filelist = {'A082113001.nex','A082113002.nex','A082113003.nex','A082113004.nex'};
%filelist = {'F092512001.nex','F092512002.nex','F092512003.nex','F092512004.nex'};
filelist = {'F112012001.nex','F112012002.nex','F112012003.nex','F112012004.nex'};

offset = [.05 .05]; % in sec
granddata = [];
grandidx = [];
grandspikes = [];
for fileidx = 1:length(filelist)
    stro = nex2stro(findfile(char(filelist{fileidx})));
    if (stro.sum.paradigmID ~= 102)
        continue
    end
    spikeidxs = find(strncmp(stro.sum.rasterCells,'sig',3));
    hepidx = strcmp(stro.sum.rasterCells(1,:),'AD11');
    vepidx = strcmp(stro.sum.rasterCells(1,:),'AD12');
    starttimeidx = strcmp(stro.sum.rasterCells(1,:),'anlgStartTime');
    stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_on'));
    stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_off'));
    e1t = [stro.ras{:,starttimeidx}]';
    stimir = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_ir'));
    fpx = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fp_x'));
    fpy = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fp_y'));
    uniquefppos = unique([fpx fpy],'rows');  % These are guaranteed to be sorted
    bkgnd = unique(stro.trial(:,strcmp(stro.sum.trialFields(1,:),'bkgnd')));
    samplerate = stro.sum.analog.storeRates{1};
    ntrials = size(stro.ras,1);
    % figure; axes; hold on;
    data = nan*ones(ntrials,ceil(max(stimoff_t+offset(2)-stimon_t-offset(1))*samplerate),2);  % Eye position snippet
    spikecounts = nan*ones(ntrials,length(spikeidxs));
    for i = 1:ntrials
        h = stro.ras{i,hepidx}.*4096/400;
        v = stro.ras{i,vepidx}.*4096/400;
        if (any(isnan(h)))
            keyboard
        end
        i
        % 4096 A/D levels = 10 V
        % 1 degree = 40 A/D levels (according to REX)
        % (4096 levels/10 V) * (1 degree/40 levels) = 4096/400.
        % The key is that both REX and PLEXON use 12 bit A/D boards configured
        % for +/- 5 V.
        
        neyesamp = size(h,1);
        x = linspace(e1t(i),neyesamp/samplerate+e1t(i),neyesamp);
        L = x >=stimon_t(i)+offset(1) & x <=stimoff_t(i)+offset(2);
        data(i,1:sum(L),1) = h(L);
        data(i,1:sum(L),2) = v(L);
        for j = 1:length(spikeidxs)
            spikecounts(i,j) = sum (stro.ras{i,spikeidxs(j)} > stimon_t(i)+offset(1) & stro.ras{i,spikeidxs(j)} < stimoff_t(i) + offset(2));
        end
    end
    grandidx = [grandidx; [fpx fpy (fpx == uniquefppos(1,1) & fpy == uniquefppos(1,2)) stro.trial(:,strcmp(stro.sum.trialFields(1,:),'bkgnd')) stimir repmat(fileidx,length(fpx),1)]];
    granddata = cat(1,granddata, data);
    grandspikes = cat(1,grandspikes, spikecounts);
end
% Columns of grandidx:
% 1) FPX   2) FPY   3) FPlocidx  4) bkgndidx   5) stimir   6) fileidx
% bkgndidx = 0 means corridor
%%
% Continued from above
% Just analysis here on down - no more disk acccess
if (size(grandspikes,2) > 1)
    whichspike = input(['Which spike? (n=',num2str(size(grandspikes,2)),')']);
else
    whichspike = 1;
end
% Now doing some analyses
Lfp = grandidx(:,3) == 1; % 1 = near
Lcorr = grandidx(:,4) == 0; % 0 means corridor background, 2 means gray.

% Generating a size tuning curve
figure; axes; hold on;
ringrad = unique(grandidx(:, 5));
for i = 0:1 % FP position
    for j = 0:1 % On the corridor?
        tmp = [];
        for k = 1:length(ringrad)
            Lbase = Lfp == i & Lcorr == j;
            L = grandidx(:,5) == ringrad(k);
            tmp(k,1) = mean(grandspikes(L&Lbase,whichspike));
            tmp(k,2) = std(grandspikes(L&Lbase,whichspike));
        end
       % plot([ringrad ringrad]'/10,[tmp(:,1)+tmp(:,2) tmp(:,1)-tmp(:,2)]','k-');
        h = plot(ringrad/10,tmp(:,1),'k-o','LineWidth',2);
        if (i == 1) % If in the near part of the corridor
           set(h,'Color','magenta') 
        end
        if (j == 0)  % if on gray bkgnd
            set(h,'LineStyle','--'); 
        end
    end
end
xlabel('ring radius (°)');
ylabel('spike count');
title('Dashed = gray background, pink = near');


% Let's just start with average fixation position in X and Y
% And pos vs time
avgpos = squeeze(nanmean(granddata,2));
x = linspace(offset(1),size(granddata,2)/samplerate,size(granddata,2));
figure; 
subplot(2,2,1); hold on;
plot(avgpos(Lfp&~Lcorr,1),avgpos(Lfp&~Lcorr,2),'o','Color',[.75 .75 .75]) % Gray means "on the corridor"
plot(avgpos(Lfp&Lcorr,1),avgpos(Lfp&Lcorr,2),'ko');
    plot([mean(avgpos(Lfp&~Lcorr,1)), mean(avgpos(Lfp&~Lcorr,1))+.01*stro.sum.exptParams.rf_x],...
         [mean(avgpos(Lfp&~Lcorr,2)), mean(avgpos(Lfp&~Lcorr,2))+.01*stro.sum.exptParams.rf_y],'b-','LineWidth',2);
plot([mean(avgpos(Lfp&~Lcorr,1)), mean(avgpos(Lfp&Lcorr,1))], [mean(avgpos(Lfp&~Lcorr,2)), mean(avgpos(Lfp&Lcorr,2))],'r-','LineWidth',2)
title('black = corridor, gray = gray bkgnd','Color','magenta');
axis square;
subplot(2,2,2); hold on;
plot(avgpos(~Lfp&~Lcorr,1),avgpos(~Lfp&~Lcorr,2),'o','Color',[.75 .75 .75])
plot(avgpos(~Lfp&Lcorr,1),avgpos(~Lfp&Lcorr,2),'ko')
plot([mean(avgpos(~Lfp&~Lcorr,1)), mean(avgpos(~Lfp&~Lcorr,1))+.01*stro.sum.exptParams.rf_x],...
    [mean(avgpos(~Lfp&~Lcorr,2)), mean(avgpos(~Lfp&~Lcorr,2))+.01*stro.sum.exptParams.rf_y],'b-','LineWidth',2);
plot([mean(avgpos(~Lfp&~Lcorr,1)), mean(avgpos(~Lfp&Lcorr,1))], [mean(avgpos(~Lfp&~Lcorr,2)), mean(avgpos(~Lfp&Lcorr,2))],'r-','LineWidth',2)
title('black = corridor, gray = gray bkgnd','Color','black');
axis square;

% Scaling the size of the mean fixation position to show the spike count

subplot(2,2,3); hold on;
for i = find(Lfp&Lcorr)'
    plot(avgpos(i,1),avgpos(i,2),'o','MarkerFaceColor','black','Color','black','MarkerSize',grandspikes(i,whichspike)*10./max(grandspikes(:,whichspike))+1);
end
for i = find(Lfp&~Lcorr)'
    plot(avgpos(i,1),avgpos(i,2),'o','MarkerFaceColor',[.65 .65 .65],'Color',[.65 .65 .65],'MarkerSize',grandspikes(i,whichspike)*10./max(grandspikes(:,whichspike))+1); % Gray means "on the gray bkgnd"
end
axis square;


subplot(2,2,4); hold on;
for i = find(~Lfp&Lcorr)'
    plot(avgpos(i,1),avgpos(i,2),'o','MarkerFaceColor','black','Color','black','MarkerSize',grandspikes(i,whichspike)*10./max(grandspikes(:,whichspike))+1);
end
for i = find(~Lfp&~Lcorr)'
    plot(avgpos(i,1),avgpos(i,2),'o','MarkerFaceColor',[.65 .65 .65],'Color',[.65 .65 .65],'MarkerSize',grandspikes(i,whichspike)*10./max(grandspikes(:,whichspike))+1); % Gray means "on the gray bkgnd"
end
axis square;

% This time course analysis didn't turn out to be very useful
% % commenting it out for now.
% subplot(2,2,3); hold on;
% plot(x,nanmean(granddata(Lfp&~Lcorr,:,1)),'r-');
% plot(x,nanmean(granddata(Lfp&Lcorr,:,1)),'m-');
% plot(x,nanmean(granddata(Lfp&~Lcorr,:,2)),'g-');
% plot(x,nanmean(granddata(Lfp&Lcorr,:,2)),'c-');
% subplot(2,2,4); hold on;
% plot(x,nanmean(granddata(~Lfp&~Lcorr,:,1)),'r-');
% plot(x,nanmean(granddata(~Lfp&Lcorr,:,1)),'m-');
% plot(x,nanmean(granddata(~Lfp&~Lcorr,:,2)),'g-');
% plot(x,nanmean(granddata(~Lfp&Lcorr,:,2)),'c-');

% Displaying some text
fprintf([filelist{1},' : ',filelist{end},'   Spike #',num2str(whichspike),'\n'])
disp([]);
% Permutation test
niter = 2000;
tmpdata = zeros(niter+1,2);
for fploc = 1:-1:0
    for i = 1:niter+1
        tmp = avgpos(Lfp == fploc,:);
        L = Lcorr(Lfp == fploc);
        if (i > 1)
            L = L(randperm(length(L)));
        end
        v = [mean(tmp(~L,1))-mean(tmp(L,1)), mean(tmp(~L,2))-mean(tmp(L,2))];
        tmpdata(i,:) = v;
    end
    dist2 = sum(tmpdata.^2,2);
    p = sum(dist2(2:end)>dist2(1))./niter;
    if (fploc == 1)
        fprintf(['Near: ']);
    else
        fprintf(['Far: ']);
    end
    disp(['amp: ',num2str(sqrt(dist2(1))),' angle: ',num2str(atan2(tmpdata(1,2),tmpdata(1,1))*180/pi), ' p: ',num2str(p)])
end

% Regression models
% Doing this on a fploc by fploc basis. Doesn't makes sense to model 
% both fixation point locations together
for fploc = 1:-1:0 % first is near, second is far.
    [p,t,stats] = anovan(grandspikes(Lfp == fploc,whichspike),[grandidx(Lfp == fploc,5) avgpos(Lfp == fploc,:)],'varnames',{'stimir','EPx','EPy'},'model','interaction','continuous',[2 3],'display','off');
    p = anovan(stats.resid, grandidx(Lfp == fploc,[4 5]),'varnames',{'bkgnidx','stimir'},'model','interaction','display','off');
    if (fploc == 1)
        disp(['Near: ',num2str(p')]);
    else
        disp(['Far: ',num2str(p')]);    
    end
    [p,t,stats]  = anovan(grandspikes(Lfp == fploc,whichspike),[grandidx(Lfp == fploc,[4 5]) avgpos(Lfp == fploc,:)],'varnames',{'bkgndidx','stimir','EPx','EPy'},'model','full','continuous',[3 4],'display','off');
    if (fploc == 1)
        disp('Near:');
        disp(t(:,[1 3 5 7]));
    else
        disp('Far:');
        disp(t(:,[1 3 5 7]))
    end
    disp(' ');
end

%p = anovan(grandspikes(:,whichspike),grandidx(:,[3 4 5]),'varnames',{'FPloc','bkgnidx','stimir'},'model','interaction');
% Spike count as a fucntion of Predictors: 3) FPlocidx  4) bkgndidx  5) stimir 


%%
% Permutation tests on every neighboring pair of fp locations
niter = 2000;
metadata =[];
for k = 0:1 % 1 = corridor background
    for i = 1:length(uniquefpx)
        Lx = grandidx(:,1) == uniquefpx(i);
        uniquefpy = unique(grandidx(Lx,2));
        for j = 1:length(uniquefpy)-1
            Ly = grandidx(:,2) == uniquefpy(j);
            sample1 = [zeros(sum(Lx&Ly&Lcorr==k),1) avgpos(Lx&Ly&Lcorr==k,1),avgpos(Lx&Ly&Lcorr==k,2)];
            Ly = grandidx(:,2) == uniquefpy(j+1);
            sample2 = [ones(sum(Lx&Ly&Lcorr==k),1) avgpos(Lx&Ly&Lcorr==k,1),avgpos(Lx&Ly&Lcorr==k,2)];
            sample = [sample1; sample2];
            data=nan*ones(niter,2);
            for iter = 1:niter
                if (iter == 1)
                    L = logical(sample(:,1));
                else
                    L = logical(sample(randperm(size(sample,1)),1));
                end
                data(iter,:) = [mean(sample(L,2)-sample(~L,2)) mean(sample(L,3)-sample(~L,3))];
            end
            stat = sqrt(sum(data.^2,2))
            p = sum(stat>=stat(1))./niter;
            metadata = [metadata; p data(1,:) (uniquefpy(j+1)-uniquefpy(j))/10]
        end
    end
end


%%
% Identical to the above cell but now working on a bigger filelist

% Note regarding eye position units
% 4096 A/D levels = 10 V
% 1 degree = 40 A/D levels (according to REX)
% (4096 levels/10 V) * (1 degree/40 levels) = 4096/400.
% The key is that both REX and PLEXON use 12 bit A/D boards configured
% for +/- 5 V.

[fnames, spikenums] = fnamesFromTxt2('/Volumes/NO BACKUP/NexFiles/nexfilelists/Greg/AmyFreya.txt');
daylist = [];
currentdate = [];
for i = 1:length(fnames)
    tmp = char(fnames{i});
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

offset = [.05 .05]; % in sec
metadata = nan*ones(daylist(end),7);
for daycounter = 1:daylist(end)
    granddata = [];
    grandidx = [];
    for fileidx = min(find(daylist==daycounter)):max(find(daylist==daycounter))
        stro = nex2stro(findfile(char(fnames{fileidx})));
        if (stro.sum.paradigmID ~= 102)
            continue
        end
        hepidx = strcmp(stro.sum.rasterCells(1,:),'AD11');
        vepidx = strcmp(stro.sum.rasterCells(1,:),'AD12');
        starttimeidx = strcmp(stro.sum.rasterCells(1,:),'anlgStartTime');
        stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_on'));
        stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_off'));
        stimx_near = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimx_near'));
        fpx = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fp_x'));
        fpy = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fp_y'));
        uniquefppos = unique([fpx fpy],'rows');  % These are guaranteed to be sorted
        bkgnd = unique(stro.trial(:,strcmp(stro.sum.trialFields(1,:),'bkgnd')));
        samplerate = stro.sum.analog.storeRates{1};
        ntrials = size(stro.ras,1);
        % figure; axes; hold on;
        data = nan*ones(ntrials,ceil(max(stimoff_t+offset(2)-stimon_t-offset(1))*samplerate),2);  % Eye position snippet
        for i = 1:ntrials
            h = stro.ras{i,hepidx}.*4096/400;
            v = stro.ras{i,vepidx}.*4096/400; 
            e1t = stro.ras{i,starttimeidx};
            neyesamp = size(h,1);
            x = linspace(e1t,neyesamp/samplerate+e1t,neyesamp);
            L = x >=stimon_t(i)+offset(1) & x <=stimoff_t(i)+offset(2);
            data(i,1:sum(L),1) = h(L);
            data(i,1:sum(L),2) = v(L);
        end
        grandidx = [grandidx; [fpx fpy (fpx == uniquefppos(1,1) & fpy == uniquefppos(1,2)) stro.trial(:,strcmp(stro.sum.trialFields(1,:),'bkgnd')) repmat([nneurons fileidx],length(fpx),1)]];
        granddata = cat(1,granddata, data);
    end
    
    % Now doing some analyses
    Lfp = grandidx(:,3) == 1;  % 1 = lower position (near)
    Lcorr = grandidx(:,4) == 0;  % 0 means corridor. 2 means gray. I double checked.
    
    % Let's just start with average fixation position in X and Y
    % And pos vs time
    avgpos = squeeze(nanmean(granddata,2));
    x = linspace(offset(1),size(granddata,2)/samplerate,size(granddata,2));
    figure;
    subplot(2,2,1); hold on;
    plot(avgpos(Lfp&~Lcorr,1),avgpos(Lfp&~Lcorr,2),'ro');
    plot(avgpos(Lfp&Lcorr,1),avgpos(Lfp&Lcorr,2),'go'); % Green means "on the corridor"
    plot([mean(avgpos(Lfp&~Lcorr,1)), mean(avgpos(Lfp&~Lcorr,1))+.01*stro.sum.exptParams.rf_x],...
         [mean(avgpos(Lfp&~Lcorr,2)), mean(avgpos(Lfp&~Lcorr,2))+.01*stro.sum.exptParams.rf_y],'b-','LineWidth',2);
    plot([mean(avgpos(Lfp&~Lcorr,1)), mean(avgpos(Lfp&Lcorr,1))], [mean(avgpos(Lfp&~Lcorr,2)), mean(avgpos(Lfp&Lcorr,2))],'k-','LineWidth',2);
    axis square;
    
    subplot(2,2,2); hold on;
    plot(avgpos(~Lfp&~Lcorr,1),avgpos(~Lfp&~Lcorr,2),'ro')
    plot(avgpos(~Lfp&Lcorr,1),avgpos(~Lfp&Lcorr,2),'go')
    plot([mean(avgpos(~Lfp&~Lcorr,1)), mean(avgpos(~Lfp&~Lcorr,1))+.01*stro.sum.exptParams.rf_x],...
         [mean(avgpos(~Lfp&~Lcorr,2)), mean(avgpos(~Lfp&~Lcorr,2))+.01*stro.sum.exptParams.rf_y],'b-','LineWidth',2);
    plot([mean(avgpos(~Lfp&~Lcorr,1)), mean(avgpos(~Lfp&Lcorr,1))], [mean(avgpos(~Lfp&~Lcorr,2)), mean(avgpos(~Lfp&Lcorr,2))],'k-','LineWidth',2)
    axis square;
    
    subplot(2,2,3); hold on;
    tmpmn = nanmean(nanmean(granddata(Lfp,:,1)));
    plot(x,nanmean(granddata(Lfp&~Lcorr,:,1))-tmpmn,'r-');
    plot(x,nanmean(granddata(Lfp&Lcorr,:,1))-tmpmn,'m-');
   
    tmpmn = nanmean(nanmean(granddata(Lfp,:,2)));
    plot(x,nanmean(granddata(Lfp&~Lcorr,:,2))-tmpmn,'g-');
    plot(x,nanmean(granddata(Lfp&Lcorr,:,2))-tmpmn,'c-');
    
    subplot(2,2,4); hold on;
    tmpmn = nanmean(nanmean(granddata(~Lfp,:,1)));
    plot(x,nanmean(granddata(~Lfp&~Lcorr,:,1))-tmpmn,'r-');
    plot(x,nanmean(granddata(~Lfp&Lcorr,:,1))-tmpmn,'m-');
    tmpmn = nanmean(nanmean(granddata(~Lfp,:,2)));
    plot(x,nanmean(granddata(~Lfp&~Lcorr,:,2))-tmpmn,'g-');
    plot(x,nanmean(granddata(~Lfp&Lcorr,:,2))-tmpmn,'c-');
    set(gcf,'Name', char(fnames{fileidx}));

    % Permutation test
    niter = 2000;
    tmpdata = zeros(niter+1,2);
    for fploc = 1:-1:0
        for i = 1:niter+1
            tmp = avgpos(Lfp == fploc,:);
            L = Lcorr(Lfp == fploc);
            if (i > 1)
                L = L(randperm(length(L)));
            end
            v = [mean(tmp(L,1))-mean(tmp(~L,1)), mean(tmp(L,2))-mean(tmp(~L,2))];
            tmpdata(i,:) = v;
        end
        dist2 = sum(tmpdata.^2,2);
        p = sum(dist2(2:end)>dist2(1))./niter;
        disp(char(fnames{fileidx}));
        disp(['amp: ',num2str(sqrt(dist2(1))),' angle: ',num2str(atan2(tmpdata(1,2),tmpdata(1,1))*180/pi), ' p: ',num2str(p)])
        disp('------------');
        metadata(daycounter, (1:3)+3*fploc) = [sqrt(dist2(1)) atan2(tmpdata(1,2),tmpdata(1,1))*180/pi p];  % far 1:3, near 4:6
    end
    stro.sum.exptParams.eye_pos_control
    if (isnan(stro.sum.exptParams.eye_pos_control) | stro.sum.exptParams.eye_pos_control == 0)
        metadata(daycounter, 7) = 0;
    else
        metadata(daycounter, 7) = 1;
    end
end
% Angles are from average saccade endpoint on gray to average saccade
% endpoint on corridor. Angular histograms.  
% Angle is from gray EP centroid to corridor EP centroid
figure;
for control = 0:1
    L = metadata(:,7 ) == control;
    subplot(2,2,2+2*control);
    rose(metadata(L,2)*pi/180);
    subplot(2,2,1+2*control);
    rose(metadata(L,5)*pi/180);
end

% Compass plots
figure;
for control = 0:1
    L = metadata(:,7 ) == control;
    subplot(2,2,2+2*control); polar(0,.2,'w'); hold on; 
    x = metadata(L,1).*cos(metadata(L,2)*pi/180);
    y = metadata(L,1).*sin(metadata(L,2)*pi/180);
    compass(x,y);
    subplot(2,2,1+2*control); polar(0,.2,'w'); hold on; 
    x = metadata(L,4).*cos(metadata(L,5)*pi/180);
    y = metadata(L,4).*sin(metadata(L,5)*pi/180);
    compass(x,y);
end

%%
% Putting together some data for Scott. For each neuron pulling out the 
% mean eye positions (on the gray and corridor backgrounds) and the
% direction of the RF. Also grabbing the ring sizes so we can see how much
% the local curvature changes.

[fnames, spikenums] = fnamesFromTxt2('/Volumes/NO BACKUP/NexFiles/nexfilelists/Greg/AmyFreya.txt');
daylist = [];
currentdate = [];
for i = 1:length(fnames)
    tmp = char(fnames{i});
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

offset = [.05 .05]; % in sec
metadata = [];
for daycounter = 1:daylist(end)
    granddata = [];
    grandidx = [];
    for fileidx = min(find(daylist==daycounter)):max(find(daylist==daycounter))
        stro = nex2stro(findfile(char(fnames{fileidx})));
        if (stro.sum.paradigmID ~= 102)
            continue
        end
        filedatestr = stro.sum.fileName(end-12:end-4);
        month = str2num(stro.sum.fileName(end-12:end-11))
        day = str2num(stro.sum.fileName(end-10:end-9))
        year = str2num(stro.sum.fileName(end-8:end-7))
        if (month >= 8 & day >=21 & year >=13)
            EPoffset = [.058 .055];  % something strange happened to rig 2 on or before 8/21
            disp('correcting eye position');
        else
            EPoffset = [0 0];
        end
        
        nneurons = sum(strncmp(stro.sum.rasterCells,'sig',3));
        hepidx = strcmp(stro.sum.rasterCells(1,:),'AD11');
        vepidx = strcmp(stro.sum.rasterCells(1,:),'AD12');
        starttimeidx = strcmp(stro.sum.rasterCells(1,:),'anlgStartTime');
        stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_on'));
        stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_off'));
       % stimx_near = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimx_near'));
        stimir = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_ir'));
        stimor = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_or'));

        fpx = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fp_x'));
        fpy = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fp_y'));
        uniquefppos = unique([fpx fpy],'rows');  % These are guaranteed to be sorted
        bkgnd = unique(stro.trial(:,strcmp(stro.sum.trialFields(1,:),'bkgnd')));
        samplerate = stro.sum.analog.storeRates{1};
        ntrials = size(stro.ras,1);
        % figure; axes; hold on;
        data = nan*ones(ntrials,ceil(max(stimoff_t+offset(2)-stimon_t-offset(1))*samplerate),2);  % Eye position snippet
        

        for i = 1:ntrials
            h = (stro.ras{i,hepidx}+EPoffset(1)).*4096/400;
            v = (stro.ras{i,vepidx}+EPoffset(2)).*4096/400; 
            e1t = stro.ras{i,starttimeidx};
            neyesamp = size(h,1);
            x = linspace(e1t,neyesamp/samplerate+e1t,neyesamp);
            L = x >=stimon_t(i)+offset(1) & x <=stimoff_t(i)+offset(2);
            data(i,1:sum(L),1) = h(L);
            data(i,1:sum(L),2) = v(L);
        end
        if ~isempty(granddata)
            data(:,size(granddata,2)+1:end,:) = []; % Cutting off any trials that went too long
        end
        rfx = stro.sum.exptParams.rf_x;
        rfy = stro.sum.exptParams.rf_y;
        Lfp = (fpx == uniquefppos(1,1) & fpy == uniquefppos(1,2));
        bkgnd = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'bkgnd'));
        grandidx = [grandidx; [fpx fpy Lfp bkgnd repmat([rfx rfy],length(Lfp),1) stimir stimor repmat(fileidx,length(fpx),1)]];
        granddata = cat(1,granddata, data);
    end
    
    % Now doing some analyses
    Lfp = grandidx(:,3) == 1;  % 1 = lower position (near)
    Lcorr = grandidx(:,4) == 0;  % 0 means corridor. 2 means gray. I double checked.
    
    % Let's just start with average fixation position in X and Y
    % And pos vs time
    avgpos = squeeze(nanmean(granddata,2));
    Lfp = grandidx(:,3) == 1;  % 1 = lower position (near)
    Lcorr = grandidx(:,4) == 0;  % 0 means corridor. 2 means gray. I double checked.
    fpnear = [unique(grandidx(Lfp,1)) unique(grandidx(Lfp,2))]/10;
    fpfar = [unique(grandidx(~Lfp,1)) unique(grandidx(~Lfp,2))]/10;
    rfpos = [unique(grandidx(:,5)) unique(grandidx(:,6))]/10;
    minringir = min(grandidx(:,7))/10;
    maxringir = min(grandidx(:,8))/10;
    if (isnan(stro.sum.exptParams.eye_pos_control) | stro.sum.exptParams.eye_pos_control == 0)
        whichgeometry = 0;
    else
        whichgeometry = 1;
    end
    
    metadata(daycounter,:) = [fpnear mean(avgpos(Lfp&~Lcorr,:)) mean(avgpos(Lfp&Lcorr,:)) fpfar mean(avgpos(~Lfp&~Lcorr,:)) mean(avgpos(~Lfp&Lcorr,:)) rfpos whichgeometry nneurons]
end
% columns
% 1) near FPx
% 2) near FPy
% 3) near EPx (gray)
% 4) near EPy (gray)
% 5) near EPx (corridor)
% 6) near EPy (corridor)
% 7) far FPx
% 8) far FPy
% 9) far EPx (gray)
% 10) far EPy (gray)
% 11) far EPx (corridor)
% 12) far EPy (corridor)
% 13) RFx
% 14) RFy
% 15) which geometry 0 = standard, 1 = eye movement control
% 16) number of neurons recorded in this session

nearEPdispX = metadata(:,5)-metadata(:,3);
nearEPdispY = metadata(:,6)-metadata(:,4);
farEPdispX = metadata(:,11)-metadata(:,9);
farEPdispY = metadata(:,12)-metadata(:,10);
nearEPdisplacement = sqrt(nearEPdispX.^2 + nearEPdispY.^2);
farEPdisplacement = sqrt(farEPdispX.^2 + farEPdispY.^2);
% computing a few statistics to quantify the magnitude of the eye position
% displacement vector
unitvectortowardsRF = mkbasis(metadata(:,[13 14])')';
% projection onto vector pointing toward RF
nearproj = sum([nearEPdispX nearEPdispY].*unitvectortowardsRF,2);
farproj = sum([farEPdispX farEPdispY].*unitvectortowardsRF,2);

% Like above, but representing mean eye position relative to the fixation
% point position. First on corridor, then on gray.
L = metadata(:,15) == 0; % 0 = Standard geometry only
nneurons = metadata(L,16);
%unitvectortowardsRF = mkbasis(metadata(L,[13 14])')';
rfxy = [metadata(L,13) metadata(L,14)];
cornearEPdispX = metadata(L,5)-metadata(L,1);
cornearEPdispY = metadata(L,6)-metadata(L,2);
corfarEPdispX = metadata(L,11)-metadata(L,7);
corfarEPdispY = metadata(L,12)-metadata(L,8);
graynearEPdispX = metadata(L,3)-metadata(L,1);
graynearEPdispY = metadata(L,4)-metadata(L,2);
grayfarEPdispX = metadata(L,9)-metadata(L,7);
grayfarEPdispY = metadata(L,10)-metadata(L,8);
tmp = [];
for j = 1:size(nneurons,1)  % Ugly nested loop
    for k = 1:nneurons(j)
        tmp = [tmp; rfxy(j,:) cornearEPdispX(j) cornearEPdispY(j) corfarEPdispX(j) corfarEPdispY(j) graynearEPdispX(j) graynearEPdispY(j) grayfarEPdispX(j) grayfarEPdispY(j)];
    end
end


% How much does the orientation of the tangent to the ring change
% due to changes in eye position?
% These are RF positions relative to the fixation point location
rfneargray = [metadata(:,13)+metadata(:,3)-metadata(:,1) metadata(:,14)+metadata(:,4)-metadata(:,2)]
rfnearcorridor = [metadata(:,13)+metadata(:,5)-metadata(:,1) metadata(:,14)+metadata(:,6)-metadata(:,2)]
thetaneargray = (pi-(pi/2+abs(atan(rfneargray(:,2)./rfneargray(:,1)))))*180/pi; % tangent to ring (deg)
thetanearcorridor = (pi-(pi/2+abs(atan(rfnearcorridor(:,2)./rfnearcorridor(:,1)))))*180/pi;
% Let see hw much the ring tangent changes with the mean shift in EP
thetaneargray-thetanearcorridor

% Again but now for the far fixation point
% These are RF positions relative to the fixation point location
rffargray = [metadata(:,13)+metadata(:,9)-metadata(:,7) metadata(:,14)+metadata(:,10)-metadata(:,8)]
rffarcorridor = [metadata(:,13)+metadata(:,11)-metadata(:,7) metadata(:,14)+metadata(:,12)-metadata(:,8)]
thetafargray = (pi-(pi/2+abs(atan(rffargray(:,2)./rffargray(:,1)))))*180/pi; % tangent to ring (deg)
thetafarcorridor = (pi-(pi/2+abs(atan(rffarcorridor(:,2)./rffarcorridor(:,1)))))*180/pi;
% Let see hw much the ring tangent changes with the mean shift in EP
thetafargray-thetafarcorridor

%EPstats = [nearEPdisplacement farEPdisplacement nearproj farproj abs(thetaneargray-thetanearcorridor) abs(thetafargray-thetafarcorridor)];

nearXgray = metadata(:,3);
nearYgray = metadata(:,4);
nearXcorridor = metadata(:,5);
nearYcorridor = metadata(:,6);
farXgray = metadata(:,9);
farYgray = metadata(:,10);
farXcorridor = metadata(:,11);
farYcorridor = metadata(:,12);

% Below, corridor relative to gray
v1 = [nearXcorridor-nearXgray nearYcorridor-nearYgray];
v2 = [farXcorridor-farXgray farYcorridor-farYgray];



EPstats = sum([v1-v2].*unitvectortowardsRF,2)
% Negative numbers mean that the monkey is looking down at the back and up at
% the front of the corridor. Positive numbers are consistent with the shift in the size
% tuning curves we've been seeing. Most of the numbers are negative because
% the monkey's behavior tends to work against the predicted shift.
% nearproj-farproj 
% This is the same as the thing above.

% Making a plot of eye position displacements between gray and corridor
% (separately for near and far and the two stimulus geometries - concentric
% and eccentric)
Lepcontrol = metadata(:,15) == 1;
figure;
AXLIMS = .2;
for i = 1:4
    subplot(2,2,i);
    if (i == 1)
        tmp = v1(~Lepcontrol,:);
    elseif (i == 2)
        tmp = v2(~Lepcontrol,:);
    elseif (i == 3)
        tmp = v1(Lepcontrol,:);
    elseif (i == 4)
        tmp = v2(Lepcontrol,:);
    end
    plot([zeros(size(tmp,1)) tmp(:,1)]',[zeros(size(tmp,1)), tmp(:,2)]','k-')
    set(gca,'Xlim',[-1 1]*AXLIMS,'Ylim',[-1 1]*AXLIMS); axis square;
    % A Rayleigh test
    unitvects = tmp./repmat(sqrt(sum(tmp.^2,2)),1,2)
    xy = mean(unitvects);
    r = norm(xy);
    n = size(unitvects,1)
    R = n*r;
    p = exp(sqrt(1+4*n+4*(n^2-R^2))-(1+2*n));
    title(num2str(p));
end

% A two-sample test to determine whether the angular distribution of mean eye
% position displacements depends on the stimulus geometry (concentric vs
% eccentric).

[th1,r1] = cart2pol(v1(~Lepcontrol,1), v1(~Lepcontrol,2));
[th2,r2] = cart2pol(v1(Lepcontrol,1), v1(Lepcontrol,2));
[~,p] = WatsonWheelerTest(th1,th2)  % near position comparison
[~,p] = ttest2(r1,r2)

[th1,r1] = cart2pol(v2(~Lepcontrol,1), v2(~Lepcontrol,2));
[th2,r2] = cart2pol(v2(Lepcontrol,1), v2(Lepcontrol,2));
[~,p] = WatsonWheelerTest(th1,th2) % far position comparison
[~,p] = ttest2(r1,r2)
[mean(r1) mean(r2)]  % For Freya eye displacements are little bigger in the
% standard version of the task than the eye position control version.


%save('ForScott','EPstats');
% % Projecting onto a few different unit vectors
% rs = []; ps = [];
% for i = 0:pi/10:2*pi
%     EPstats = [v1-v2]*[cos(i); sin(i)];
%     [r,p] = corrcoef([shift EPstats])
%     rs(end+1) = r(1,2);
%     ps(end+1) = p(1,2);
% end
% 

%%% Plot for Scott
Lepcontrol = metadata(:,15) == 1;
figure;
AXLIMS = 3;
for i = 1:4
    for j = 1:2
        if (j == 1)
            %  Below, all near (first gray then corridor)
            v1 = [metadata(:,3)-metadata(:,1) metadata(:,4)-metadata(:,2)];
            v2 = [metadata(:,5)-metadata(:,1) metadata(:,6)-metadata(:,2)];
        else
            % Below, all far (first gray then corridor) 
            v1 = [metadata(:,9)-metadata(:,7) metadata(:,10)-metadata(:,8)];
            v2 = [metadata(:,11)-metadata(:,7) metadata(:,12)-metadata(:,8)];
        end
        subplot(2,2,i); hold on;
        if (i == 1)
            tmp = v1(~Lepcontrol,:);
        elseif (i == 2)
            tmp = v2(~Lepcontrol,:);
        elseif (i == 3)
            tmp = v1(Lepcontrol,:);
        elseif (i == 4)
            tmp = v2(Lepcontrol,:);
        end
        h = plot(tmp(:,1),tmp(:,2),'ko');
        if (j == 1)
            set(h,'MarkerFaceColor','blue');
        else
            set(h,'MarkerFaceColor','red');    
        end
        set(gca,'Xlim',[-1 1]*AXLIMS,'Ylim',[-1 1]*AXLIMS); axis square;
        plot(0,0,'k+','LineWidth',2)
    end
end

% Number of units
sum(metadata(metadata(:,15) == 0,16))
sum(metadata(metadata(:,15) == 1,16))

%% Looking at fix move data

filelist = {'F090912001.nex','F090912002.nex','F090912003.nex','F090912004.nex','F090912005.nex','F090912006.nex'}% Fix move

offset = [.05 .05]; % in sec
granddata = [];
grandidx = [];
grandspikes = [];
for fileidx = 1:length(filelist)
    stro = nex2stro(findfile(char(filelist{fileidx})));
    if (stro.sum.paradigmID ~= 102)
        continue
    end
    spikeidxs = find(strncmp(stro.sum.rasterCells,'sig',3));
    hepidx = strcmp(stro.sum.rasterCells(1,:),'AD11');
    vepidx = strcmp(stro.sum.rasterCells(1,:),'AD12');
    starttimeidx = strcmp(stro.sum.rasterCells(1,:),'anlgStartTime');
    stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_on'));
    stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_off'));
    e1t = [stro.ras{:,starttimeidx}]';
    stimir = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_ir'));
    fpx = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fp_x'));
    fpy = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fp_y'));
    uniquefppos = unique([fpx fpy],'rows');  % These are guaranteed to be sorted
    bkgnd = unique(stro.trial(:,strcmp(stro.sum.trialFields(1,:),'bkgnd')));
    samplerate = stro.sum.analog.storeRates{1};
    ntrials = size(stro.ras,1);
    % figure; axes; hold on;
    data = nan*ones(ntrials,ceil(max(stimoff_t+offset(2)-stimon_t-offset(1))*samplerate),2);  % Eye position snippet
    spikecounts = nan*ones(ntrials,length(spikeidxs));
    for i = 1:ntrials
        h = stro.ras{i,hepidx}.*4096/400;
        v = stro.ras{i,vepidx}.*4096/400;
        if (any(isnan(h)))
            keyboard
        end
        i
        % 4096 A/D levels = 10 V
        % 1 degree = 40 A/D levels (according to REX)
        % (4096 levels/10 V) * (1 degree/40 levels) = 4096/400.
        % The key is that both REX and PLEXON use 12 bit A/D boards configured
        % for +/- 5 V.
        
        neyesamp = size(h,1);
        x = linspace(e1t(i),neyesamp/samplerate+e1t(i),neyesamp);
        L = x >=stimon_t(i)+offset(1) & x <=stimoff_t(i)+offset(2);
        data(i,1:sum(L),1) = h(L);
        data(i,1:sum(L),2) = v(L);
        for j = 1:length(spikeidxs)
            spikecounts(i,j) = sum (stro.ras{i,spikeidxs(j)} > stimon_t(i)+offset(1) & stro.ras{i,spikeidxs(j)} < stimoff_t(i) + offset(2));
        end
    end
    grandidx = [grandidx; [fpx fpy (fpx == uniquefppos(1,1) & fpy == uniquefppos(1,2)) stro.trial(:,strcmp(stro.sum.trialFields(1,:),'bkgnd')) stimir repmat(fileidx,length(fpx),1)]];
    granddata = cat(1,granddata, data);
    grandspikes = cat(1,grandspikes, spikecounts);
end

Lcorr = grandidx(:,4) == 0; % 0 means corridor background, 2 means gray.
avgpos = squeeze(nanmean(granddata,2));
fp_pos = unique(grandidx(:,[1 2]),'rows');
uniquefpx = unique(fp_pos(:,1));
cmap = colormap(jet(size(fp_pos(:,2),1)/2));
data = nan*ones(2,2,2,size(cmap,1));

% Plotting individual fixation positions in different colors
figure;
for k = 0:1 % 1 = corridor background
    for i = 1:length(uniquefpx)
        subplot(2,2,i+k*2); hold on;
        Lx = grandidx(:,1) == uniquefpx(i);
        uniquefpy = unique(grandidx(Lx,2));        
        for j = 1:length(uniquefpy)
            Ly = grandidx(:,2) == uniquefpy(j);
            h = plot(avgpos(Lx&Ly&Lcorr==k,1),avgpos(Lx&Ly&Lcorr==k,2),'o');
            set(h,'Color',cmap(j,:));
            centroid = [mean(avgpos(Lx&Ly&Lcorr==k,1)), mean(avgpos(Lx&Ly&Lcorr==k,2))];
            h = plot(uniquefpx(i)/10,uniquefpy(j)/10,'s','MarkerSize',5); 
            set(h,'MarkerFaceColor',cmap(j,:),'MarkerEdgeColor','black');
            h = plot(centroid(1),centroid(2),'p','MarkerSize',12); 
            set(h,'MarkerFaceColor',cmap(j,:)/2,'MarkerEdgeColor','black');
            
            data(:,k+1,i,j) = centroid;
        end
        axis equal
        if (k)
            title('corridor');
        else
            title('gray');
        end
    end
end

% The the plot above, but drawing 1 SD ellipses around each mean
STD = 1;                     %# 2 standard deviations
conf = 2*normcdf(STD)-1;     %# covers around 95% of population
scale = chi2inv(conf,2);     %# inverse chi-squared with dof=#dimensions
figure;
for k = 0:1 % 1 = corridor background
    for i = 1:length(uniquefpx)
        subplot(2,2,i+k*2); hold on;
        Lx = grandidx(:,1) == uniquefpx(i);
        uniquefpy = unique(grandidx(Lx,2));        
        for j = 1:length(uniquefpy)
            Ly = grandidx(:,2) == uniquefpy(j);
            
            % Plotting centroid
            centroid = [mean(avgpos(Lx&Ly&Lcorr==k,1)), mean(avgpos(Lx&Ly&Lcorr==k,2))];
            h = plot(centroid(1),centroid(2),'p','MarkerSize',12);
            set(h,'MarkerFaceColor',cmap(j,:)/2,'MarkerEdgeColor','black');
            
            % Plotting a 2 SD ellipse
            S = cov(avgpos(Lx&Ly&Lcorr==k,[1 2]));
            [v,d] = eig(S);
            [d, order] = sort(diag(d), 'descend');
            d = diag(d)*scale;
            v = v(:, order);
            t = linspace(0,2*pi,100);
            e = [cos(t) ; sin(t)];        %# unit circle
            vv = v*sqrt(d);               %# scale eigenvectors
            e = bsxfun(@plus, vv*e, centroid'); %#' project circle back to orig space
            h = plot(e(1,:), e(2,:), 'LineWidth',2);
            set(h,'Color',cmap(j,:)/1.5);
   
        end
        axis equal
        if (k)
            title('corridor');
        else
            title('gray');
        end
    end
end

% Need to do this last to make sure the nominal FP locations end up on top
for k = 0:1 % 1 = corridor background
    for i = 1:length(uniquefpx)
        subplot(2,2,i+k*2); hold on;
        Lx = grandidx(:,1) == uniquefpx(i);
        uniquefpy = unique(grandidx(Lx,2));        
        for j = 1:length(uniquefpy)
            % Plotting nominal fp position
            h = plot(uniquefpx(i)/10,uniquefpy(j)/10,'s','MarkerSize',5); 
            set(h,'MarkerFaceColor',cmap(j,:),'MarkerEdgeColor','black');         
        end
    end
end



% Y position as a function of y FP position
figure;
for k = 0:1 % 1 = corridor background
    for i = 1:length(uniquefpx)
        subplot(2,2,i+k*2); hold on;
        Lx = grandidx(:,1) == uniquefpx(i);
        uniquefpy = unique(grandidx(Lx,2));
        for j = 1:length(uniquefpy)
            Ly = grandidx(:,2) == uniquefpy(j);
            plot(uniquefpy(j)./10,avgpos(Lx&Ly&Lcorr==k,2),'ko');
        end
        plot([uniquefpy(1) uniquefpy(end)]/10,[uniquefpy(1) uniquefpy(end)]/10,'g-','LineWidth',2);
        axis square;
        set(gca,'Xlim',[uniquefpy(1)-5 uniquefpy(end)+5]/10);
        set(gca,'Ylim',[uniquefpy(1)-5 uniquefpy(end)+5]/10);
        if (k)
            title('corridor');
        else
            title('gray');
        end
    end
end

%%
% Population analysis of fix move data
% Also a supplementary figure for the paper


[fnames, spikenums] = fnamesFromTxt2('/Volumes/NO BACKUP/NexFiles/nexfilelists/Greg/AmyApolloFixMove.txt');
daylist = [];
currentdate = [];
for i = 1:length(fnames)
    tmp = char(fnames{i});
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

offset = [.05 .05]; % in sec
metadata =[]; % Currently unused
metadata2 = []; % eye positions for Scott
for daycounter = 1:daylist(end)
    granddata = [];
    grandidx = [];
    for fileidx = min(find(daylist==daycounter)):max(find(daylist==daycounter))
        stro = nex2stro(findfile(char(fnames{fileidx})));
        if (stro.sum.paradigmID ~= 102)
            continue
        end
        nneurons = sum(strncmp(stro.sum.rasterCells,'sig',3));
        filedatestr = stro.sum.fileName(end-12:end-4);
        month = str2num(stro.sum.fileName(end-12:end-11));
        day = str2num(stro.sum.fileName(end-10:end-9));
        year = str2num(stro.sum.fileName(end-8:end-7));
        if (month >= 8 & day >=21 & year >=13 | month >= 9 & day >= 17 & year >=13) % Ugly
            EPoffset = [.058 .055];  % something strange happened to rig 2 on or before 8/21, resolved from 8/28-9/8, then bad again 9/9-???
            disp('correcting eye position');
        else
            EPoffset = [0 0];
        end
        spikeidxs = find(strncmp(stro.sum.rasterCells,'sig',3));
        hepidx = strcmp(stro.sum.rasterCells(1,:),'AD11');
        vepidx = strcmp(stro.sum.rasterCells(1,:),'AD12');
        starttimeidx = strcmp(stro.sum.rasterCells(1,:),'anlgStartTime');
        stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_on'));
        stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_off'));
        e1t = [stro.ras{:,starttimeidx}]';
        stimir = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_ir'));
        fpx = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fp_x'));
        fpy = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fp_y'));
        uniquefppos = unique([fpx fpy],'rows');  % These are guaranteed to be sorted
        bkgnd = unique(stro.trial(:,strcmp(stro.sum.trialFields(1,:),'bkgnd')));
        samplerate = stro.sum.analog.storeRates{1};
        ntrials = size(stro.ras,1);
        % figure; axes; hold on;
        data = nan*ones(ntrials,ceil(max(stimoff_t+offset(2)-stimon_t-offset(1))*samplerate),2);  % Eye position snippet
        spikecounts = nan*ones(ntrials,length(spikeidxs));
        for i = 1:ntrials
            h = (stro.ras{i,hepidx}+EPoffset(1)).*4096/400;
            v = (stro.ras{i,vepidx}+EPoffset(2)).*4096/400;
            if (any(isnan(h)))
                keyboard
            end
            
            % 4096 A/D levels = 10 V
            % 1 degree = 40 A/D levels (according to REX)
            % (4096 levels/10 V) * (1 degree/40 levels) = 4096/400.
            % The key is that both REX and PLEXON use 12 bit A/D boards configured
            % for +/- 5 V.
            
            neyesamp = size(h,1);
            x = linspace(e1t(i),neyesamp/samplerate+e1t(i),neyesamp);
            L = x >=stimon_t(i)+offset(1) & x <=stimoff_t(i)+offset(2);
            data(i,1:sum(L),1) = h(L);
            data(i,1:sum(L),2) = v(L);
            for j = 1:length(spikeidxs)
                spikecounts(i,j) = sum (stro.ras{i,spikeidxs(j)} > stimon_t(i)+offset(1) & stro.ras{i,spikeidxs(j)} < stimoff_t(i) + offset(2));
            end
        end
        grandidx = [grandidx; [fpx fpy (fpx == uniquefppos(1,1) & fpy == uniquefppos(1,2)) stro.trial(:,strcmp(stro.sum.trialFields(1,:),'bkgnd')) stimir repmat(fileidx,length(fpx),1)]];
        if (~isempty(granddata) & size(granddata,2) > size(data,2))
            data = [data, nan*ones(size(data,1),size(granddata,2)-size(data,2),size(data,3))];
        end
        if (~isempty(granddata) & size(granddata,2) < size(data,2))
            granddata = [granddata, nan*ones(size(granddata,1),size(data,2)-size(granddata,2),size(granddata,3))];
        end
        granddata = cat(1,granddata, data);
    end % End of filelist loop

    Lcorr = grandidx(:,4) == 0; % 0 means corridor background, 2 means gray.
    avgpos = squeeze(nanmean(granddata,2));
    fp_pos = unique(grandidx(:,[1 2]),'rows');
    uniquefpx = unique(fp_pos(:,1));

    % Y position as a function of y FP position
    figure;
    for k = 0:1 % 1 = corridor background
        for i = 1:length(uniquefpx)
            subplot(2,2,i+k*2); hold on;
            Lx = grandidx(:,1) == uniquefpx(i);
            uniquefpy = unique(grandidx(Lx,2));
            for j = 1:length(uniquefpy)
                Ly = grandidx(:,2) == uniquefpy(j);
                plot(uniquefpy(j)./10,avgpos(Lx&Ly&Lcorr==k,2),'ko');
            end
            plot([uniquefpy(1) uniquefpy(end)]/10,[uniquefpy(1) uniquefpy(end)]/10,'g-','LineWidth',2);
            axis square;
            set(gca,'Xlim',[uniquefpy(1)-5 uniquefpy(end)+5]/10);
            set(gca,'Ylim',[uniquefpy(1)-5 uniquefpy(end)+5]/10);
            if (k)
                title('corridor');
            else
                title('gray');
            end
        end
    end
    
    % Getting some numbers together for Scott
    tmp1 = [];
    for i = 1:-1:0 % 1 = corridor background
        for j = 1:length(uniquefpx)
            Lx = grandidx(:,1) == uniquefpx(j);
            tmp2 = [];
            uniquefpy = unique(grandidx(Lx,2));
            for k = 1:length(uniquefpy)           
                Ly = grandidx(:,2) == uniquefpy(k);
                tmp2 =[tmp2; avgpos(Lx&Ly&Lcorr==i,:)-grandidx(Lx&Ly&Lcorr==i,[1 2])/10];
            end
            if isempty(tmp2)
                keyboard
            end
            tmp1 = [tmp1; mean(tmp2)]; % cor-near, cor-far, gray-near, gray-far
        end
    end
    for i = 1:nneurons
        metadata2 = [metadata2; stro.sum.exptParams.rf_x/10 stro.sum.exptParams.rf_y/10 reshape(tmp1,1,8)];
    end

    % Permutation tests on every neighboring pair of fp locations
    niter = 2000;
    for k = 0:1 % 1 = corridor background
        for i = 1:length(uniquefpx)
            Lx = grandidx(:,1) == uniquefpx(i);
            uniquefpy = unique(grandidx(Lx,2));
            for j = 1:length(uniquefpy)-1
                Ly = grandidx(:,2) == uniquefpy(j);
                sample1 = [zeros(sum(Lx&Ly&Lcorr==k),1) avgpos(Lx&Ly&Lcorr==k,1),avgpos(Lx&Ly&Lcorr==k,2)];
                Ly = grandidx(:,2) == uniquefpy(j+1);
                sample2 = [ones(sum(Lx&Ly&Lcorr==k),1) avgpos(Lx&Ly&Lcorr==k,1),avgpos(Lx&Ly&Lcorr==k,2)];
                sample = [sample1; sample2];
                data=nan*ones(niter,2);
                for iter = 1:niter
                    if (iter == 1)
                        L = logical(sample(:,1));
                    else
                        L = logical(sample(randperm(size(sample,1)),1));
                    end
                    data(iter,:) = [mean(sample(L,2))-mean(sample(~L,2)) mean(sample(L,3))-mean(sample(~L,3))];
                end
                stat = sqrt(sum(data.^2,2))
                p = sum(stat>=stat(1))./niter;
                metadata = [metadata; p data(1,:) (uniquefpy(j+1)-uniquefpy(j))/10 i k daycounter]
            end
        end
    end
end

% Getting some data for Scott
% [RFx, RFy, corr_near_x, corr_near_y, corr_far_x, corr_far_y, gray_near_x, gray_near_y, gray_far_x, gray_far_y]
figure; axes; hold on;
for i = 3:2:size(metadata2,2)
    if i == 1
        symb = 'ks';
    elseif i == 3
        symb = 'k.';
    elseif i == 5
        symb = 'gs';
    else
        symb = 'g.';
    end
    plot(metadata2(:,i),metadata2(:,i+1),symb)
end

% How many fixation point locations used on each day in the "front" and the
% "back" positions. Seems like almost always 7 verically offset fixation 
% point locations (except one file). This gives 6 comparisons per file.
nfpcomparisons = [];
Lcor = logical(metadata(:,6));
for i = 1:metadata(end,end)
    Lday = metadata(:,end) == i;
    nfpcomparisons(i,1) = sum(Lday&~Lcor&metadata(:,5)==1); % near position, gray
    nfpcomparisons(i,2) = sum(Lday&~Lcor&metadata(:,5)==2); % far position, gray
    nfpcomparisons(i,3) = sum(Lday&Lcor&metadata(:,5)==1); % near position, corridor
    nfpcomparisons(i,4) = sum(Lday&Lcor&metadata(:,5)==2); % far position, corridor
end

% How many times do we reject the null hypothesis of same eye position distribution
% as a function of separation between fixation point y offsets?

uniquefpy = unique(metadata(:,4));
for i = 1:length(uniquefpy)
    L = metadata(:,4) == uniquefpy(i);
    uniquefpy(i)
    [sum(metadata(L&Lcor,1) < 0.05) sum(L&Lcor)]
    [sum(metadata(L&~Lcor,1) < 0.05) sum(L&~Lcor)]
end


figure; subplot(2,2,1);
hist(metadata(:,1),30);
subplot(2,2,2);
hist(metadata(:,3)-metadata(:,4),30);
subplot(2,2,3);
hist(metadata(:,3)./metadata(:,4),30);
geomean(metadata(:,3)./metadata(:,4))
geomean(metadata(Lcor,3)./metadata(Lcor,4))
geomean(metadata(~Lcor,3)./metadata(~Lcor,4))
[h,p] = ttest2(metadata(Lcor,3)./metadata(Lcor,4),metadata(~Lcor,3)./metadata(~Lcor,4))

% Making a plot where each mean eye position displacement is represented by
% a vector. Only looking at the 0.2 degree comparisons. Only looking at
% corridor data.
L = metadata(:,6) == 1; % Corridor background only
L = L & metadata(:,4) == .2; % Magnitude of fp displacements to consider
Lfar = metadata(:,5) == 2;
AXLIMS = .3;
figure;
subplot(2,2,1);
plot([zeros(sum(L&~Lfar),1) metadata(L&~Lfar,2)]',[zeros(sum(L&~Lfar),1)  metadata(L&~Lfar,3)]','k-');
set(gca,'Xlim',[-1 1]*AXLIMS,'Ylim',[-1 1]*AXLIMS); axis square;
subplot(2,2,2);
plot([zeros(sum(L&Lfar),1) metadata(L&Lfar,2)]',[zeros(sum(L&Lfar),1)  metadata(L&Lfar,3)]','k-');
set(gca,'Xlim',[-1 1]*AXLIMS,'Ylim',[-1 1]*AXLIMS); axis square;

% Figure for the paper
figure;
whichmag = 0.2;
Lmag = metadata(:,4) == whichmag; % Magnitude of fp displacements to consider
Lfar = metadata(:,5) == 2;
for backgnd = 0:1 %1 means corridor
    subplot(2,2,1+backgnd); hold on;
    Lcor = logical(metadata(:,6)) == backgnd;
    for loc = 1:2 % 2 = far
        if (loc == 1)
            plot(metadata(Lcor&Lmag&~Lfar,2),metadata(Lcor&Lmag&~Lfar,3),'bo','Markerfacecolor','blue');
        else
            plot(metadata(Lcor&Lmag&Lfar,2),metadata(Lcor&Lmag&Lfar,3),'ro','Markerfacecolor','red');            
        end
    end
    if (backgnd == 0)
        title('near');
    else
        title('corridor');
    end
    axis square;
    
    if strcmp(stro.sum.fileName(end-13),'F')
        extent = 0.3;
    else
        extent = 0.5;
    end
    set(gca,'YLim',[0 extent],'Xlim',[-1 1]*extent/2);
    plot([-1 1]*extent/2,[1 1]* whichmag,'k--');
    plot([0 0],[0, extent],'k--');
end
% axes are delta mean eye positions (one fixation point position versus the
% next).

% How many sessions?
length(unique(metadata(Lmag, end)))


%%
% Microsaccade analysis


[fnames, spikenums] = fnamesFromTxt2('/Volumes/NO BACKUP/NexFiles/nexfilelists/Greg/AmyApollo.txt');
daylist = [];
currentdate = [];
for i = 1:length(fnames)
    tmp = char(fnames{i});
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

offset = [.05 .05]; % in sec
granddata = [];
grandidx = [];
for daycounter = 1:daylist(end)
    for fileidx = min(find(daylist==daycounter)):max(find(daylist==daycounter))
        stro = nex2stro(findfile(char(fnames{fileidx})));
        if (stro.sum.paradigmID ~= 102)
            continue
        end
        sacstats = getSacData(stro); close;
        stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_on'));
        stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_off'));
        stimx_near = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimx_near'));
        fpx = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fp_x'));
        fpy = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fp_y'));
        uniquefppos = unique([fpx fpy],'rows');  % These are guaranteed to be sorted
        bkgnd = unique(stro.trial(:,strcmp(stro.sum.trialFields(1,:),'bkgnd')));
        if (isnan(stro.sum.exptParams.eye_pos_control) | stro.sum.exptParams.eye_pos_control == 0)
            whichgeometry = 0;
        else
            whichgeometry = 1;
        end
        
        ntrials = size(stro.ras,1);
        for i = 1:ntrials
            [stimon_t(i)+offset(1) sacstats.starttimes{i}' stimoff_t(i)+offset(2)]
           
            L = sacstats.starttimes{i} > stimon_t(i)+offset(1) & sacstats.starttimes{i} < stimoff_t(i)+offset(2);
            if any(L)
                disp('found one!');
                for j = find(L)
                    j
                    grandidx = [grandidx; fpx(i) fpy(i) (fpx(i) == uniquefppos(1,1) & fpy(i) == uniquefppos(1,2)) stro.trial(i,strcmp(stro.sum.trialFields(1,:),'bkgnd')) whichgeometry daycounter];
                    granddata = [granddata; stimon_t(i) stimoff_t(i) sacstats.starttimes{i}(j) sacstats.directions{i}(j) sacstats.amplitudes{i}(j)];
                end
            end
        end
    end
end
% Indices into grandidx
% 1) fpx
% 2) fpy
% 3) near or far fixation point locations (1 = near)
% 4) What kind of background (0 = corridor, 2 = gray)
% 5) which geometry (0 = concentric, 1 = eccentric stimulus)
% 6) index for day on which files were recorded

% Indices into granddata
% 1) stim on time
% 2) stim off time
% 3) saccade start time
% 4) saccade direction
% 5) saccade amplitude

% First, amplitudes
figure;
subplotcounter = 1;
metadata = [];
EPCONTROL = 0;  % 0 = regular geometry, 1 = eye positon control.
for i = 1:-1:0 % 1 = near, 0 = far
    for j = [2 0] % 2 = gray, 0 = corridor
        L = grandidx(:,3) == i & grandidx(:,4) == j;
        L = L & grandidx(:,5) == EPCONTROL;  
        sum(L)
        subplot(2,2,subplotcounter);
        hist(granddata(L,5));
        if (subplotcounter == 1)
            title('near, gray');
        elseif (subplotcounter == 2)
            title('near, corridor');
        elseif (subplotcounter == 3)
            title('far, gray');
        elseif (subplotcounter == 4)
            title('far, corridor');
        end
        ylabel('count'); xlabel('µsac amp.');
        metadata = [metadata; granddata(L,5) repmat([i j],sum(L),1)];
        subplotcounter = subplotcounter+1;
    end
end
p =anovan(metadata(:,1),metadata(:,[2 3]),'model','interaction','varnames',{'near/far','bkgnd'});


% Second, directions
figure;
subplotcounter = 1;
metadata = [];
for i = 1:-1:0 % 1 = near, 0 = far
    for j = [2 0] % 2 = gray, 0 = corridor
        L = grandidx(:,3) == i & grandidx(:,4) == j;
        L = L & grandidx(:,5) == EPCONTROL;
        sum(L)
        subplot(2,2,subplotcounter);
        [t,r] = rose(granddata(L,4));
        h = polar(t,r,'k-');
        if (subplotcounter == 1)
            title('near, gray');
        elseif (subplotcounter == 2)
            title('near, corridor');
        elseif (subplotcounter == 3)
            title('far, gray');
        elseif (subplotcounter == 4)
            title('far, corridor');
        end
        metadata = [metadata; granddata(L,4) repmat([i j],sum(L),1)];
        subplotcounter = subplotcounter+1;
    end
end

% Third, directions and amplitudes jointly
figure;
subplotcounter = 1;
for i = 1:-1:0 % 1 = near, 0 = far
    for j = [2 0] % 2 = gray, 0 = corridor
        L = grandidx(:,3) == i & grandidx(:,4) == j;
        L = L & grandidx(:,5) == EPCONTROL;
        subplot(2,2,subplotcounter);
        [x,y] = pol2cart(granddata(L,4), granddata(L,5))
        compass(x,y);

        if (subplotcounter == 1)
            title('near, gray');
        elseif (subplotcounter == 2)
            title('near, corridor');
        elseif (subplotcounter == 3)
            title('far, gray');
        elseif (subplotcounter == 4)
            title('far, corridor');
        end
        subplotcounter = subplotcounter+1;
    end
end

% histogram of when microsaccades occur
subplotcounter = 1;
metadata =[];
for i = 1:-1:0 % 1 = near, 0 = far
    for j = [2 0] % 2 = gray, 0 = corridor
        L = grandidx(:,3) == i & grandidx(:,4) == j;
        L = L & grandidx(:,5) == EPCONTROL;
        subplot(2,2,subplotcounter);
        hist(granddata(L,3)-granddata(L,1),[0:.01:.3]);
        set(gca,'Xlim',[0 .3]);
        mean(granddata(L,3)-granddata(L,1))
        
        if (subplotcounter == 1)
            title('near, gray');
        elseif (subplotcounter == 2)
            title('near, corridor');
        elseif (subplotcounter == 3)
            title('far, gray');
        elseif (subplotcounter == 4)
            title('far, corridor');
        end
        xlabel('time from stimon (s)');
        ylabel('µ saccade count');
        
        subplotcounter = subplotcounter+1;
        metadata = [metadata; granddata(L,3)-granddata(L,1) repmat([i j],sum(L),1)];
    end
end
p =anovan(metadata(:,1),metadata(:,[2 3]),'model','interaction','varnames',{'near/far','bkgnd'});
mean(metadata(metadata(:,3) == 2),1) % gray
mean(metadata(metadata(:,3) == 0),1) % corridor
%%
% Getting RF positions, ring diameters, monitor distances
% stro.,trial columns 12 and 13 are inner radius and outer radius.
[fnames, spikenums] = fnamesFromTxt2('/Volumes/NO BACKUP/NexFiles/nexfilelists/Greg/AmyFreyaFixMove.txt');
data = [];
for fileidx = 1:length(fnames)
    stro = nex2stro(findfile(char(fnames{fileidx})));
    if (stro.sum.paradigmID ~= 102)
        continue
    end
    rfx = stro.sum.exptParams.rf_x/10;
    rfy = stro.sum.exptParams.rf_y/10;
    mondist = stro.sum.exptParams.mondist;
    eyewinx = 2*stro.sum.exptParams.eyewin_x/10;
    eyewiny = 2*stro.sum.exptParams.eyewin_y/10;
    
    ir = stro.trial(:,12)/10;
    or = stro.trial(:,13)/10;
    ringwidths = unique(or-ir);
    data = [data; rfx rfy eyewinx eyewiny mondist mean(ringwidths)]
end

mean(data)
std(data)


%%
% Getting data for time course analysis

EYEPOSCONTROL = 1;

% Good parameters
aw = .2; ws = 0.05; sw = 11;

aw = .2; % large analysis window time, starting at stimulus onset
ws = 0.02; % small window length
sw = 31; % number of time bins
binleftedges = linspace(-ws,aw,sw)


for listidx = 1:2
    if listidx == 1
        listnames = {'AmyFreya.txt','AmyFreyaFixMove.txt'};
    else
        
        listnames = {'AmyApollo.txt','AmyApolloFixMove.txt'};
    end
    %listnames = {'AmyFreya.txt'};
    
    fnames = [];
    for i = 1:length(listnames)
        fnames = [fnames; fnamesFromTxt2(['/Volumes/NO BACKUP/NexFiles/nexfilelists/Greg/',listnames{i}])];
    end
    
    daylist = [];
    currentdate = [];
    for i = 1:length(fnames)
        tmp = char(fnames{i});
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
    
    radii = []; name = [];
    Gnear = []; Gfar = []; GnearSEM = []; GfarSEM = []; rowGall = [];
    Cnear = []; Cfar = []; CnearSEM = []; CfarSEM = []; rowCall = [];
    rowGall = 1; rowCall = 1;
    for daycounter = 1:daylist(end)
        nearG=[]; farG=[]; nearC=[]; farC=[];
        rowG = 1; % Gray
        rowC = 1; % Corridor
        for fileidx = find(daylist==daycounter, 1 ):find(daylist==daycounter, 1, 'last' )
            stro = nex2stro(findfile(char(fnames{fileidx})));
            if (stro.sum.paradigmID ~= 102)
                continue
            end
            if (EYEPOSCONTROL == 1 & stro.sum.exptParams.eye_pos_control == 1) | (EYEPOSCONTROL == 0 & stro.sum.exptParams.eye_pos_control ~= 1)
                [bkgnd, near, far, uniqueors, spikename] = SMurray_unit_spikes_all(stro, binleftedges); %spike rates per small analysis windows
                if bkgnd == 2 %gray background
                    nearG(rowG:rowG+size(near,1)-1, :, :, :) = near;
                    farG(rowG:rowG+size(far,1)-1, :, :, :) = far;
                    rowG = rowG + size(near,1);
                elseif bkgnd == 0 %corridor background
                    nearC(rowC:rowC+size(near,1)-1, :, :, :) = near;
                    farC(rowC:rowC+size(far,1)-1, :, :, :) = far;
                    rowC = rowC + size(near,1);
                end
            end
        end
        % order of dimensions in "near" 1:Trial, 2:TimeBin ,3:OuterRadius ,4:Spike
        
        % come here once per day (not once per neuron)
        if (~isempty(nearG) & ~isempty(nearC))
            for k = 1:size(nearG,4)
                %store radii of rings used for each unit
                radii(rowGall, :) = uniqueors/10;
                name(rowGall, :) = spikename(k,:);
                
                %store mean data and z-scored data for each unit (GRAY)
                Gnear(rowGall, :, :) = mean(nearG(:,:,:,k));
                Gfar(rowGall, :, :) = mean(farG(:,:,:,k));
                GnearSEM(rowGall, :, :) = std(nearG(:,:,:,k))./sqrt(size(nearG,1));
                GfarSEM(rowGall, :, :)  = std(farG(:,:,:,k))./sqrt(size(farG,1));
                rowGall = rowGall + 1;
                
                %store mean data and z-scored data for each unit (CORRIDOR)
                Cnear(rowCall, :, :) = mean(nearC(:,:,:,k));
                Cfar(rowCall, :, :) = mean(farC(:,:,:,k));
                CnearSEM(rowCall, :, :) = std(nearC(:,:,:,k))./sqrt(size(nearC,1));
                CfarSEM(rowCall, :, :)  = std(farC(:,:,:,k))./sqrt(size(farC,1));
                rowCall = rowCall + 1;
            end
        end
    end
    % Saving some data for Scott
    savefilename = 'MXEX';
    if findstr(listnames{1},'Freya')
        savefilename(2) = '1';
    else
        savefilename(2) = '2';
    end
    
    if (EYEPOSCONTROL == 1)
        savefilename(4) = '2';
    else
        savefilename(4) = '1';
    end
    
    
    save (savefilename, 'Cfar', 'CfarSEM', 'Cnear', 'CnearSEM', 'Gfar', 'GfarSEM', 'Gnear', 'GnearSEM', 'binleftedges', 'radii')
    
end  % end of monkey (list) loop


% dimensions of Gnear: 1) Neuron, 2) TimeBin, 3) RingDiameter

% Just some simple analysis to make sure I've done this correctly
% UNNORMALIZED PSTHS
% figure; % PSTHs as a function of ring diam and FP position (CORRIDOR)
% for i = 1:size(Cnear,3)
%     subplot(size(Cnear,3),1,i); hold on;
%     plot(binleftedges, mean(Cnear(:,:,i))','b-') % PSTH at "i"th ring size
%     plot(binleftedges, mean(Cfar(:,:,i))','r-') % PSTH at "i"th ring size
%     set(gca,'Xlim',[min(binleftedges) max(binleftedges)],'Ylim',[10 70]);
% end
% 
% figure; % PSTHs as a function of ring diam and FP position (GRAY)
% for i = 1:size(Gnear,3)
%     subplot(size(Gnear,3),1,i); hold on;
%     plot(binleftedges, mean(Gnear(:,:,i))','b-') % PSTH at "i"th ring size
%     plot(binleftedges, mean(Gfar(:,:,i))','r-') % PSTH at "i"th ring size
%     set(gca,'Xlim',[min(binleftedges) max(binleftedges)],'Ylim',[10 70]);
% end

% NORMALIZED PSTHS
figure; 
set(gcf,'Position',[440 70 240 730],'PaperPositionMode','auto'); % PSTHs as a function of ring diam and FP position (CORRIDOR)
for i = 1:size(Cnear,3)
    subplot(size(Cnear,3),1,i); hold on;
    plot(binleftedges, mean(Cnear(:,:,i))./max([mean(Cnear(:,:,i)) mean(Cfar(:,:,i))]),'b.-') % PSTH at "i"th ring size
    plot(binleftedges, mean(Cfar(:,:,i))./max([mean(Cnear(:,:,i)) mean(Cfar(:,:,i))]),'r.-') % PSTH at "i"th ring size
    set(gca,'Xlim',[min(binleftedges) max(binleftedges)],'Ylim',[0 1]);
end

figure; % PSTHs as a function of ring diam and FP position (GRAY)
set(gcf,'Position',[440 70 240 730],'PaperPositionMode','auto'); % PSTHs as a function of ring diam and FP position (CORRIDOR)
for i = 1:size(Gnear,3)
    subplot(size(Gnear,3),1,i); hold on;
    plot(binleftedges, mean(Gnear(:,:,i))./max([mean(Gnear(:,:,i)) mean(Gfar(:,:,i))]),'b.-') % PSTH at "i"th ring size
    plot(binleftedges, mean(Gfar(:,:,i))./max([mean(Gnear(:,:,i)) mean(Gfar(:,:,i))]),'r.-') % PSTH at "i"th ring size
    set(gca,'Xlim',[min(binleftedges) max(binleftedges)],'Ylim',[0 1]);
end

% Differences between near and far
figure; % PSTHs as a function of ring diam and FP position (GRAY)
set(gcf,'Position',[440 70 240 730],'PaperPositionMode','auto'); % PSTHs as a function of ring diam and FP position (CORRIDOR)
for i = 1:size(Gnear,3)
    subplot(size(Gnear,3),1,i); hold on;
    normnear = mean(Gnear(:,:,i))./max([mean(Gnear(:,:,i)) mean(Gfar(:,:,i))]);
    normfar = mean(Gfar(:,:,i))./max([mean(Gnear(:,:,i)) mean(Gfar(:,:,i))]);
    plot(binleftedges, normnear-normfar,'k.-') % PSTH at "i"th ring size

    normnear = mean(Cnear(:,:,i))./max([mean(Cnear(:,:,i)) mean(Cfar(:,:,i))]);
    normfar = mean(Cfar(:,:,i))./max([mean(Cnear(:,:,i)) mean(Cfar(:,:,i))]);
    plot(binleftedges, normnear-normfar,'b.-') % PSTH at "i"th ring size

    set(gca,'Xlim',[min(binleftedges) max(binleftedges)],'Ylim',[-.5 .5]);
end



%%
% Ok, now trying to make a nice plot showing the temporal evolution of the
% RF shift. Pooling data across monkeys.
% Scott chucks neurons 1,11,43 from Freya

% loading the data from the manuscripts folder: Freya
Cnear1 = []; Cfar1 = []; Gnear1 = []; Gfar1 = [];
load('/Users/greghorwitz/Documents/Manuscripts/Size Illusion stuff/Matlab data/M1E1.mat');
if (size(Cfar,1) == 46)  % making sure this is M1E1
    Cfar1 = Cfar; Cnear1 = Cnear; Gfar1 = Gfar; Gnear1 = Gnear;
    skipcells = [1 11 43];
    Cfar1(skipcells,:,:) = [];
    Cnear1(skipcells,:,:) = [];
    Gfar1(skipcells,:,:) = [];
    Gnear1(skipcells,:,:) = [];
end

% loading the data from the manuscripts folder: Apollo
load('/Users/greghorwitz/Documents/Manuscripts/Size Illusion stuff/Matlab data/M2E1.mat');
if (size(Cfar,1) == 18)
    Cfar1 = Cfar; Cnear1 = Cnear; Gfar1 = Gfar; Gnear1 = Gnear;
    skipcells = [11 13];
    Cfar(skipcells,:,:) = [];
    Cnear(skipcells,:,:) = [];
    Gfar(skipcells,:,:) = [];
    Gnear(skipcells,:,:) = [];
end

Cnear1 = cat(1,Cnear,Cnear1);
Cfar1 = cat(1,Cfar,Cfar1);
Gnear1 = cat(1,Gnear,Gnear1);
Gfar1 = cat(1,Gfar,Gfar1);

load('/Users/greghorwitz/Documents/Manuscripts/Size Illusion stuff/Matlab data/M1E2.mat');
Cnear1 = cat(1,Cnear1,Cnear);
Cfar1 = cat(1,Cfar1,Cfar);
Gnear1 = cat(1,Gnear1,Gnear);
Gfar1 = cat(1,Gfar1,Gfar);

load('/Users/greghorwitz/Documents/Manuscripts/Size Illusion stuff/Matlab data/M2E2.mat');
Cnear1 = cat(1,Cnear1,Cnear);
Cfar1 = cat(1,Cfar1,Cfar);
Gnear1 = cat(1,Gnear1,Gnear);
Gfar1 = cat(1,Gfar1,Gfar);

Cnear = Cnear1;
Cfar = Cfar1;
Gnear = Gnear1;
Gfar = Gfar1;

% normalizing each neuron
% 4: condition, 3: ring size, 2: time.
normmat = max(max(max(cat(4,Gnear,Gfar,Cnear,Cfar),[],4),[],3),[],2);
normmat = mean(mean(mean(cat(4,Gnear,Gfar,Cnear,Cfar),4),3),2);
normmat = ones(size(normmat));

normmat = max(max(max(cat(4,Cnear,Cfar),[],4),[],3),[],2);
Cnormnear = Cnear./repmat(normmat,[1,size(Cnear,2),size(Cnear,3)]);

normmat = max(max(max(cat(4,Cfar,[]),[],4),[],3),[],2);
Cnormfar = Cfar./repmat(normmat,[1,size(Cfar,2),size(Cfar,3)]);

normmat = max(max(max(cat(4,Gnear,Gfar),[],4),[],3),[],2);
Gnormnear = Gnear./repmat(normmat,[1,size(Gnear,2),size(Gnear,3)]);

normmat = max(max(max(cat(4,Gfar,[]),[],4),[],3),[],2);
Gnormfar = Gfar./repmat(normmat,[1,size(Gfar,2),size(Gfar,3)]);

%Subtracting baseline, which appears to be different on the gray and
%corridor backgrounds
% baseline = mean(mean(Cnear(:,binleftedges < 0,:),2),3);
% Cnormnear = Cnear-repmat(baseline,[1,size(Cnear,2),size(Cnear,3)]);
% baseline = mean(mean(Cfar(:,binleftedges < 0,:),2),3);
% Cnormfar = Cfar-repmat(baseline,[1,size(Cfar,2),size(Cfar,3)]);
% baseline = mean(mean(Gnear(:,binleftedges < 0,:),2),3);
% Gnormnear = Gnear-repmat(baseline,[1,size(Gnear,2),size(Gnear,3)]);
% baseline = mean(mean(Gfar(:,binleftedges < 0,:),2),3);
% Gnormfar = Gfar-repmat(baseline,[1,size(Gfar,2),size(Gfar,3)]);



CfarM = squeeze(mean(Cnormfar));
CnearM = squeeze(mean(Cnormnear));
GfarM = squeeze(mean(Gnormfar));
GnearM = squeeze(mean(Gnormnear));
overallmin = min([CfarM(:); CnearM(:); GfarM(:); GnearM(:)]); % For plotting
overallmax = max([CfarM(:); CnearM(:); GfarM(:); GnearM(:)]); % For plotting

vars = {CnearM', CfarM', GnearM', GfarM'};

binwidth = binleftedges(2) - binleftedges(1);
titles = {'Corridor near','Corridor far','Gray near','Gray far'};
figure;colormap(jet(255));
for i = 1:4
    subplot(3,2,i);
    m = vars{i};
    image((m-overallmin)./(overallmax-overallmin)*255);
    ylabel('Ring size');
    set(gca,'Ytick',[1 4 7],'Yticklabel',{'small','medium','large'});
    xlabel('Time (ms)')
    set(gca,'Xtick',[1:3:length(binleftedges)]-0.5,'XtickLabel',num2str(binleftedges([1:3:length(binleftedges)])'*1000));
    title(titles{i});
end

% difference image for the corridor
C_Diff = CfarM - CnearM;

% difference image for gray
G_Diff = GfarM - GnearM;

%plot difference maps
% figure(1)
% imagesc(C_Diff)
% figure(2)
% imagesc(G_Diff)

%to get rid of that monitor artifact
subplot(3,2,5)
difference = C_Diff - G_Diff;
image((difference' - min(difference(:)))./(max(difference(:))-min(difference(:)))*255);% axis square;
ylabel('Ring size');
set(gca,'Ytick',[1 4 7],'Yticklabel',{'small','medium','large'});
xlabel('Time (ms)')
set(gca,'Xtick',[1:3:length(binleftedges)]-.5,'XtickLabel',num2str(binleftedges([1:3:length(binleftedges)])'*1000));
title('Diff_C - Diff_G');


% Just the corridor difference (ignoring Gray background data)
subplot(3,2,6)
difference = G_Diff;
image((difference' - min(difference(:)))./(max(difference(:))-min(difference(:)))*255);% axis square;
ylabel('Ring size');
set(gca,'Ytick',[1 4 7],'Yticklabel',{'small','medium','large'});
xlabel('Time (ms)')
set(gca,'Xtick',[1:3:length(binleftedges)]-.5,'XtickLabel',num2str(binleftedges([1:3:length(binleftedges)])'*1000));
title('Diff_C');



% % Trying a waterfall plot
% binwidth = binleftedges(2)-binleftedges(1)
% x = repmat(binleftedges+binwidth/2,7,1)
% y = repmat([1:7]',1,length(binleftedges))
% waterfall(x,y,difference')

% Trying a permutation test to test the hypothesis that there is a 
% % difference between GfarM and GnearM
% 
% niter = 2000;
% permdata = nan*ones(niter,size(Gnormfar,2),size(Gnormfar,3));
% Gnormtmp = cat(1,Gnormnear, Gnormfar);
% L = [true(size(Gnormfar,1),1); false(size(Gnormnear,1),1)];
% for i =1:niter
%     L = L(randperm(length(L)));
%     tmp1 = squeeze(mean(Gnormtmp(L,:,:)));
%     tmp2 = squeeze(mean(Gnormtmp(~L,:,:)));
%     permdata(i,:,:) = tmp1-tmp2;
% end
% 
% p = nan*ones(size(permdata,2),size(permdata,3));
% for i = 1:size(permdata,2)
%     for j = 1:size(permdata,3)
%        p(i,j) = sum(squeeze(permdata(:,i,j)) > (GfarM(i,j) - GnearM(i,j)))./niter;
%     end
% end

% grabbing a couple of slices
difference = C_Diff - G_Diff;
seofdiff = squeeze(sqrt((var(Cnormfar)./size(Cnormfar,1)^2)+...
    (var(Cnormnear)./size(Cnormnear,1)^2)+...
    (var(Gnormfar)./size(Gnormfar,1)^2)+...
    (var(Gnormnear)./size(Gnormnear,1)^2)));

% difference = C_Diff;
% seofdiff = squeeze(sqrt((var(Gnormfar)./size(Gnormfar,1)^2)+...
%     (var(Gnormnear)./size(Gnormnear,1)^2)));
% 

figure;
subplot(2,1,1); hold on;
idx = find((binleftedges -.024).^2 == min((binleftedges -.024).^2))
errorbar(difference(idx,:),seofdiff(idx,:));
%set(gca,'Ylim',[-.15 .15]);

subplot(2,1,2); hold on;
idx = find((binleftedges -.1).^2 == min((binleftedges -.1).^2))
errorbar(difference(idx,:),seofdiff(7,:));
%set(gca,'Ylim',[-.15 .15]);

figure;
for i = 1: size(CfarM,1)
    subplot(ceil(sqrt(size(CfarM,1))),ceil(sqrt(size(CfarM,1))),i);
    hold on;
    plot((CfarM(i,:)+CnearM(i,:)+GfarM(i,:)+GnearM(i,:))./5,'b-');
    plot(CfarM(i,:)-CnearM(i,:)-GfarM(i,:)+GnearM(i,:),'r-');
    set(gca,'Xlim',[1 7]); %,'Ylim',[-20 30]);
end

%%
% Working on a simple feedforward model that might explain the results in terms of a
% static RF shift that is relative to the vanishing point.

% Let's say the screen is 10x10 (-5 to 5 in both dimensions)
vanishingpoint = [-2 4];
RF = [-1 -1]; % relative to FP
fp1 = [-2 -2];
fp2 = [2 2];

fp3 = fp1-2*RF;
fp4 = fp2-2*RF;


% Let's say the "flow field" displaces RF centers toward the vanishing point
% a fixed amount
delta = 1;

figure; axes; hold on;
plot([-5 5 5 -5 -5],[-5 -5 5 5 -5],'k-')
set(gca,'Xlim',[-6 6],'ylim',[-6 6]);
plot(vanishingpoint(1),vanishingpoint(2),'m*','MarkerSize',15);

%experiment 1
plot(fp1(1),fp1(2),'k.');
plot(fp2(1),fp2(2),'k.');

plot(fp1(1)+RF(1),fp1(2)+RF(2),'ro');
v1 = vanishingpoint-fp1-RF; % getting a vector toward vanishingpoint
v1 = v1./norm(v1)*min(norm(v1),delta);
plot(fp1(1)+RF(1)+[0 v1(1)],fp1(2)+RF(2)+[0 v1(2)])

plot(fp2(1)+RF(1),fp2(2)+RF(2),'ro');
v2 = vanishingpoint-fp2-RF; % getting a vector toward vanishingpoint
v2 = v2./norm(v2)*min(norm(v2),delta);
plot(fp2(1)+RF(1)+[0 v2(1)],fp2(2)+RF(2)+[0 v2(2)])

%experiment 2
plot(fp3(1),fp3(2),'k.');
plot(fp4(1),fp4(2),'k.');

plot(fp3(1)+RF(1),fp3(2)+RF(2),'ro');
v3 = vanishingpoint-fp3-RF; % getting a vector toward vanishingpoint
v3 = v3./norm(v3)*min(norm(v3),delta);
plot(fp3(1)+RF(1)+[0 v3(1)],fp3(2)+RF(2)+[0 v3(2)])

plot(fp4(1)+RF(1),fp4(2)+RF(2),'ro');
v4 = vanishingpoint-fp4-RF; % getting a vector toward vanishingpoint
v4 = v4./norm(v4)*min(norm(v4),delta);
plot(fp4(1)+RF(1)+[0 v4(1)],fp4(2)+RF(2)+[0 v4(2)])


% Getting distances (experiment 1)
d1 = norm(fp1 - fp1+RF+v1)
d2 = norm(fp2 - fp2+RF+v2)

% Getting distances (experiment 2)
d3 = norm(fp3 - fp3+RF+v3)
d4 = norm(fp4 - fp4+RF+v4)

[d1-d2 d3-d4]
if (d1>d2)
    disp('Correctly predicts sign of experiment 1');
else
    disp('Fails to predict sign of experiment 1');    
end


if (d3>d4)
    disp('Correctly predicts sign of experiment 2');
else
    disp('Fails to predict sign of experiment 2');    
end


%%
% Trying to extract the sizes of the rings used in the psychophysical
% experiments.
% These numbers are in DVA*10
SUBJECTNUMBER = 2; % 1 = Freya, 2 = Apollo (these are the actual monkey number in the ms)
% 3 = Eric, 4 = Zack

% In the manuscript Freya is Monkey 1 (Apollo is Monkey 2)

if (SUBJECTNUMBER == 1)
%    pathname = '/VOLUMES/NO BACKUP/NexFiles/Amy/Illusion data/Data_Behavior_Freya';
    pathname = '/VOLUMES/NO BACKUP/NexFiles/Amy/Illusion data/DataFMFreyaUnits';
    pathname = '/VOLUMES/NO BACKUP/NexFiles/Amy/Illusion data/DataFSFreyaUnits';
    monkeyletter = 'F';
elseif (SUBJECTNUMBER == 2)
    pathname = '/VOLUMES/NO BACKUP/NexFiles/Amy/Illusion data/Data_Behavior_Apollo';
%    pathname = '/VOLUMES/NO BACKUP/NexFiles/Amy/Illusion data/DataFMApolloUnits';
%    pathname = '/VOLUMES/NO BACKUP/NexFiles/Amy/Illusion data/DataFSApolloUnits';
    monkeyletter = 'A';
elseif (SUBJECTNUMBER == 3)
    pathname = '/VOLUMES/NO BACKUP/NexFiles/Amy/Illusion data/Data_Behavior_E';
    monkeyletter = 'E';
elseif (SUBJECTNUMBER == 4)
    pathname = '/VOLUMES/NO BACKUP/NexFiles/Amy/Illusion data/Data_Behavior_Z';
    monkeyletter = 'Z';
end

a = dir(pathname);
data = [];
for i =1:length(a)
    name = a(i).name;
    if (name(1) == monkeyletter)
        stro = nex2stro([pathname,filesep,name]);
        if (all(isnan(stro.trial(:,[12 13]))))
            data = [data; mean(stro.trial(:,[19,20,21,22]))];
        else
            data = [data; mean(stro.trial(:,[12, 13]))];
        end
    end 
end

2*mean(data)/10  % 2 puts it in diameters, 10 puts it in degrees

%% 
% Regression analysis of size tuning curve peak shifts (corridor vs gray)
% to answer the reviewer's question: Is there a change in the slope?
% Need to first load the file "Monkey1_Correlation.mat" (or
% "Monkey2_Correlation.mat") from Scott that contain x,y coordinates in two
% variables: Corrid and Gray.

figure; axes; hold on;
plot(Gray(:,1),Gray(:,2),'bo');
plot(Corrid(:,1),Corrid(:,2),'r+');
x = [Gray(:,1); Corrid(:,1)];
n = size(Gray,1);
I = [ zeros(n,1); ones(n,1)];
X = [ones(2*n,1), x, x.*I, I];
y = [Gray(:,2); Corrid(:,2)];
[b, bint] = regress(y,X,.05);
% order of coefficients: (1)intercept, (2)slope, (3)additional slope for corridior,
% (4)additional intercept for corridor
b
bint

xplt = get(gca,'Xlim');
plot(xplt, b(1)+xplt.*b(2),'b-');
plot(xplt, b(1)+b(4)+xplt.*(b(2)+b(3)),'r-');

%%
% Comparing the slopes of psychometric functions (early vs late) files
% Most of this stuff was lifted from SMurray_p.m
% To analyze data from Apollo, go to NO BACKUP/NexFiles/Amy/Illusion
% Data/Apollo and select any file in this directory.
% I don't think I really need the gray background trials here, but
% gathering them anyway. 4/7/14 GH

[fname, pathname] = uigetfile('*.nex', 'Select a NeuroExplorer file in the right folder');
folder = pathname;
allfiles = dir(folder);
rowfiles = 1;
allblocks = [];
for i = 1:size(allfiles,1)
    if size(allfiles(i).name,2) == 14
        if allfiles(i).name(11:14) == fname(11:14)
            allblocks(rowfiles,:) = allfiles(i).name;
            rowfiles = rowfiles + 1;
        end
    end
end

% Making two big "stro" type structures
% Early files are early in the trial field of the new stro structures.
% Late files are at the end of the trial field.
row = 1;
rowC = 1;
for k = 1:size(allblocks,1)
    stro1 = nex2stro(strcat(folder, allblocks(k,:)));
    bkgnd = stro1.trial(:,strcmp(stro1.sum.trialFields(1,:),'bkgnd'));
    if bkgnd == 2 %gray
        %combining all trials with this background
        stro.sum = stro1.sum;
        stro.trial(row:row+length(stro1.trial)-1,:) = stro1.trial(1:length(stro1.trial),:);
        row = row+length(stro1.trial);
    elseif bkgnd == 0 %corridor
        %combining all trials with this background
        stroC.sum = stro1.sum;
        stroC.trial(rowC:rowC+length(stro1.trial)-1,:) = stro1.trial(1:length(stro1.trial),:);
        rowC = rowC+length(stro1.trial);
    end
end

stim_or_near = stroC.trial(:,strcmp(stroC.sum.trialFields(1,:),'stim_or_near'));
stim_or_far = stroC.trial(:,strcmp(stroC.sum.trialFields(1,:),'stim_or_far'));
%correct_t = stroC.trial(:,strcmp(stroC.sum.trialFields(1,:),'correct_resp'));
sacc_time_near = stroC.trial(:,strcmp(stro1.sum.trialFields(1,:),'sacc_time_near'));


nfdiff = stim_or_near - stim_or_far; %how much bigger the near stim is than the far stim
nfdiff = round(nfdiff*10) / 10; %rounded to hundredth of a degree
diffunique = unique(nfdiff);
nfdiff(:,2) = ~isnan(sacc_time_near);

% End of set up. Now doing the analysis
trialbreakpoints = round(linspace(0,size(nfdiff,1),4));
trial_begin_end = [];
for i = 1:length(trialbreakpoints)-1;
    trial_begin_end = [trial_begin_end; trialbreakpoints(i)+1 trialbreakpoints(i+1)];
end

data = []; % [ring diff, # of near choices, total # of trials]
for i = 1:size(trial_begin_end,1) % early, late, etc.
    Ltrials = zeros(size(nfdiff,1),1);
    Ltrials(trial_begin_end(i,1):trial_begin_end(i,2)) = 1;
    for j = 1:length(diffunique)
        Lring = nfdiff(:,1) == diffunique(j);
        data(j,:,i) = [diffunique(j) sum(nfdiff(Lring&Ltrials,2)) sum(Lring&Ltrials)];
    end
end

figure; axes; hold on;
colors = ['b','r','g','m'];
for i = 1:size(data,3)
    plot(data(:,1,i),data(:,2,i)./data(:,3,i),[colors(i),'o-'])
end
legend({'1st third','2nd third','3rd third'})
%legend({'1st quarter','2nd quarter','3rd quarter','4th quarter'})
set(gca,'Xlim',[-15 15]);

%%
% Looking at the distributions of sigmas from the Gaussian fits to the size
% tuning functions to try to get a feel for the "sizes of the RFs" (Reviewer
% 1).
% The files are called 
% One variable: DATA
% nx4 matrix: [Gray_Near Gray_Far Corr_Near Corr_Far]
% From the Gatass et al 1981 paper (Fig 13) I estimated that 
% sqrt(area) = .11*eccentricity+.813
% which gives
% Freya (Monkey 1) has RFs at ~2.25° = 1.06°
% Apollo (Monkey 2) has RFs at ~6° = 1.47° 
% (RF eccentricities estimated from poster)

for i = 1:2
    if i == 1
        load SigmaM1.mat 
    else
        load SigmaM2.mat
    end
    greg = log10(abs(DATA(:)));
    mn1 = 10^(mean(greg));
    mn2 = median(DATA(:));
    effectivediam = 5.5*mn1
    % area = (effectivediam/2).^2*pi
end



%%
% Mixed factors ANOVA looking at shifts on the corridor background.
% Need to have DATA_move and DATA_stay in the workspace
% columns of each are: gray near, gray far, corridor near, corridor far
n_stay = size(DATA_stay,1);
n_move = size(DATA_move,1);

X1 = [DATA_stay(:,3);DATA_stay(:,4);DATA_move(:,3);DATA_move(:,4)]; % data
X2 = [ones(2*n_stay,1); 2*ones(2*n_move,1)]; % Between subjects factor (stay vs move)
X3 = [ones(n_stay,1);2*ones(n_stay,1);ones(n_move,1);2*ones(n_move,1)]; % Within subjects factor (corridor near vs corridor far)
X4 = [[1:n_stay]'; [1:n_stay]'; n_stay+[1:n_move]'; n_stay+[1:n_move]']; % Cell labels

X = [X1,X2,X3,X4];

[SSQs, DFs, MSQs, Fs, Ps] = mixed_between_within_anova(X);
