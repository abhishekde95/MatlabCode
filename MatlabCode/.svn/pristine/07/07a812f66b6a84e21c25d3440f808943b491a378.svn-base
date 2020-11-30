% Population analyses SC_stimcue
% Author - Abhishek De, 6/18
close all; clearvars;
conn = database('Abhishek','horwitzlab','vector','Vendor','MySql','Server','128.95.153.12');
filename = fetch(conn,'SELECT filename FROM SCstimcue');
training = fetch(conn,'SELECT training FROM SCstimcue');
close(conn);
plot_counter = 1;
ind = find(strcmp(training,'no'));
L = ceil(sqrt(numel(ind)));
binwidth = .005;
bins = -0.4:binwidth:0.6;
withinRFlatencyL = cell(numel(ind),1);
withinRFlatencyNL = cell(numel(ind),1);
outsideRFlatencyL = cell(numel(ind),1);
outsideRFlatencyNL = cell(numel(ind),1);
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
    [C,ia,ib] = unique(stro.trial(:,12:13),'rows');
    RFloc = [RFloc; stro.sum.exptParams.rf_x stro.sum.exptParams.rf_y];
    latency = sacint_t - fpoff_t;
    latency(isnan(latency)) = 0.5; % Artificially allocating all the NaNs as 0.5 ms 
    withinRFlatencyL{ii} = latency(all(stro.trial(:,[12 13]) == RFloc,2) & Llaser & ~Lcatchtrials);
    withinRFlatencyNL{ii} = latency(all(stro.trial(:,[12 13]) == RFloc,2) & ~Llaser & ~Lcatchtrials);
    outsideRFlatencyL{ii} = latency(all(stro.trial(:,[12 13]) ~= RFloc,2) & Llaser & ~Lcatchtrials);
    outsideRFlatencyNL{ii} = latency(all(stro.trial(:,[12 13]) ~= RFloc,2) & ~Llaser & ~Lcatchtrials);
end
monkeyIdxs = cell2mat(filename(ind));
monkeyIdxs = monkeyIdxs(:,1);
idxs = monkeyIdxs == 'M'; % For Maui
bins = 0:0.02:0.5;
figure(plot_counter); set(gcf,'Name','Latency analysis');
subplot(221); histogram(cell2mat(withinRFlatencyL(idxs)),bins,'Normalization','probability','FaceColor',[0 0.5 1]); hold on; histogram(cell2mat(withinRFlatencyNL(idxs)),bins,'Normalization','probability','FaceColor',[0.5 0.5 0.5]);
plot(median(cell2mat(withinRFlatencyL(idxs))),0,'v','MarkerFaceColor',[0 0.5 1],'MarkerEdgeColor',[1 1 1]); plot(median(cell2mat(withinRFlatencyNL(idxs))),0,'v','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]);
xlabel('latency'); ylabel('pdf'); title('M: Within RF'); axis square; set(gca,'Tickdir','out'); hold off;
subplot(222); histogram(cell2mat(outsideRFlatencyL(idxs)),bins,'Normalization','probability','FaceColor',[0 0.5 1]); hold on; histogram(cell2mat(outsideRFlatencyNL(idxs)),bins,'Normalization','probability','FaceColor',[0.5 0.5 0.5]);
plot(median(cell2mat(outsideRFlatencyL(idxs))),0,'v','MarkerFaceColor',[0 0.5 1],'MarkerEdgeColor',[1 1 1]); plot(median(cell2mat(outsideRFlatencyNL(idxs))),0,'v','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]);
xlabel('latency'); ylabel('pdf'); title('M: Outside RF'); axis square; set(gca,'Tickdir','out'); hold off;
subplot(223); histogram(cell2mat(withinRFlatencyL(~idxs)),bins,'Normalization','probability','FaceColor',[0 0.5 1]); hold on; histogram(cell2mat(withinRFlatencyNL(~idxs)),bins,'Normalization','probability','FaceColor',[0.5 0.5 0.5]);
plot(median(cell2mat(withinRFlatencyL(~idxs))),0,'v','MarkerFaceColor',[0 0.5 1],'MarkerEdgeColor',[1 1 1]); plot(median(cell2mat(withinRFlatencyNL(~idxs))),0,'v','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]);
xlabel('latency'); ylabel('pdf'); title('A: Within RF'); axis square; set(gca,'Tickdir','out'); hold off;
subplot(224); histogram(cell2mat(outsideRFlatencyL(~idxs)),bins,'Normalization','probability','FaceColor',[0 0.5 1]); hold on; histogram(cell2mat(outsideRFlatencyNL(~idxs)),bins,'Normalization','probability','FaceColor',[0.5 0.5 0.5]);
plot(median(cell2mat(outsideRFlatencyL(~idxs))),0,'v','MarkerFaceColor',[0 0.5 1],'MarkerEdgeColor',[1 1 1]); plot(median(cell2mat(outsideRFlatencyNL(~idxs))),0,'v','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[1 1 1]);
xlabel('latency'); ylabel('pdf'); title('A: Outside RF'); axis square; set(gca,'Tickdir','out'); hold off;
plot_counter = plot_counter + 1;