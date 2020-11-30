% Population analyses SC_stimcue
% Author - Abhishek De, 6/18
close all; clearvars;
conn = database('Abhishek','horwitzlab','vector','Vendor','MySql','Server','128.95.153.12');
filename = fetch(conn,'SELECT filename FROM SCstimcue');
training = fetch(conn,'SELECT training FROM SCstimcue');
comments = fetch(conn,'SELECT comments FROM SCstimcue');
close(conn);
plot_counter = 1;
ind = find(strcmp(training,'no'));
figure(plot_counter); set(gcf,'Name','Laser');
figure(plot_counter); set(gcf,'Name','control');
L = ceil(sqrt(numel(ind)));
binwidth = .005;
bins = -0.4:binwidth:0.6;
for ii = 1:numel(ind)
    targetlocations = [];
    targethitsallNL = []; % non-laser trials
    targethitsallL = []; % laser trials
    RFloc = [];
    fileofinterest = char(filename(ind(ii),:));
    stro = nex2stro(findfile(fileofinterest));
    targ_x = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targ_x'));
    targ_y = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targ_y'));
    uniquetargxy = unique([targ_x targ_y],'rows');
    ntrials = size(stro.trial,1);
    Lcatchtrials = targ_x == 0 & targ_y == 0;
    laseron_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimon_t'));
    laseroff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimoff_t'));
    fpoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fpoff_t'));
    sacint_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'saccinit_t'));
    sacmade_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'saccmade_t'));
    targon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targon_t'));
    targoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targoff_t'));
    Llaser = ~isnan(laseron_t);
    samplerate = stro.sum.analog.storeRates{1};
    [C,ia,ib] = unique(stro.trial(:,12:13),'rows');
    targetlocations = [targetlocations; C];
    targethitsNL = zeros(size(C,1),1);
    targethitsL = zeros(size(C,1),1);
    RFloc = [RFloc; stro.sum.exptParams.rf_x stro.sum.exptParams.rf_y];
    if norm(RFloc)>100
        halfwidth = 5.0; % in degrees of visual angle
    else
        halfwidth = 3.0;
    end
    countL = 1; countNL = 1;
    PSTHlaser = zeros(1,length(bins));
    PSTHcontrol = zeros(1,length(bins));
    for i = 1:size(stro.trial,1)
        x = stro.ras{i,2}*4096/400;
        y = stro.ras{i,3}*4096/400;
        t = stro.ras{i,5}+[0:1:length(x)-1]/samplerate;
        Lt = t>fpoff_t(i) & t < fpoff_t(i)+0.30;
        analogstartime = stro.ras{i,5};
        spiketimes = stro.ras{i,1};
        if Llaser(i) & ~Lcatchtrials(i) & all(stro.trial(i,[12 13])==RFloc,2)
            laserontime = laseron_t(i);
            laserofftime = laseroff_t(i);
            timedurlaser = laserofftime - laserontime;
            PSTHlaser = PSTHlaser + hist(spiketimes-laserontime, bins);
            figure(plot_counter); subplot(L,L,ii); plot(spiketimes(spiketimes>laserontime-0.2 & spiketimes<laserofftime+0.2)-laserontime,countL*ones(size(spiketimes(spiketimes>laserontime-0.2 & spiketimes<laserofftime+0.2))),'k.'); hold on;
            countL = countL + 1;
        elseif ~Llaser(i) & ~Lcatchtrials(i) & all(stro.trial(i,[12 13])==RFloc,2)
            targont = targon_t(i);
            targofft = targoff_t(i);
            timedurcontrol = targofft - targont;
            PSTHcontrol = PSTHcontrol + hist(spiketimes-targont, bins);
            figure(plot_counter+1); subplot(L,L,ii); plot(spiketimes(spiketimes>targont-0.2 & spiketimes<targofft+0.2)-targont,countNL*ones(size(spiketimes(spiketimes>targont-0.2 & spiketimes<targofft+0.2))),'k.'); hold on;
            countNL = countNL + 1;
            
        end
    end
    figure(plot_counter); subplot(L,L,ii); plot(bins,PSTHlaser,'-','color',[0 0.5 1],'Linewidth',2); line([0 0],[0 countL]); line([timedurlaser timedurlaser],[0 countL]); xlabel('time'); ylabel('trials'); title(char(filename(ind(ii)))); set(gca,'Xlim',[-0.2 0.5]); hold off;
    figure(plot_counter+1); subplot(L,L,ii); plot(bins,PSTHcontrol,'-','color',[0.5 0.5 0.5],'Linewidth',2); line([0 0],[0 countNL]); line([timedurcontrol timedurcontrol],[0 countNL]); xlabel('time'); ylabel('trials'); title(char(filename(ind(ii)))); set(gca,'Xlim',[-0.2 0.5]); hold off;
end

plot_counter = plot_counter + 2;