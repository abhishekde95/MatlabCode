%analysis of data from Memsacc paradigm

stro = nex2stro;
fpon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fpon_t'));
fpacq_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fpacq_t'));
fpoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fpoff_t'));
stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimon_t'));
stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimoff_t'));
sacc_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'saccmade_t')); %time monkey left fix window
stim_xy = [stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_x')), ...
           stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_y'))]; %stimulus locations
rf_xy = [stro.sum.exptParams.rf_x, stro.sum.exptParams.rf_y]; %receptive field location
unique_xy = unique(stim_xy, 'rows'); %unique stimulus locations

%plot raster plot and PSTH
offset = [-.1 .8]; %pre and post sync-time to plot
binwidth = .005;
bins = offset(1):binwidth:offset(2);
PSTH = zeros(size(unique_xy,1),length(bins));
Lspikechans = strncmp(stro.sum.rasterCells,'sig0',4);
if (any(Lspikechans))
    for whichspike = find(Lspikechans)
        figure; 
        for i = 1:size(unique_xy,1)
            L = stim_xy == repmat(unique_xy(i,:),size(stim_xy,1), 1);
            LL = find(L(:,1) & L(:,2));
            clear spikes
            spikes = stro.ras(LL,whichspike);
            subplot(2, size(unique_xy,1), i); hold on;
            title(['XY: ' num2str(unique_xy(i,1)/10) ', ' num2str(unique_xy(i,2)/10)]);
            for j = 1:length(LL)
                clear tmpspikes
                tmpspikes = spikes{j}-stimon_t(LL(j));
                tmpspikes(tmpspikes < offset(1) | tmpspikes > offset(2)) = [];
                nspikestot = length(tmpspikes);
                plot([tmpspikes tmpspikes]',[zeros(nspikestot,1) .5*ones(nspikestot,1)]'+j,'k-');
                PSTH(i,:) = PSTH(i,:) + hist(tmpspikes, bins);
            end
            PSTH(i,:) = PSTH(i,:)/binwidth/length(LL);
            
            %plotting the time course of stimulation
            dur = mode(stimoff_t-stimon_t);
            plot([0 dur], [-3 -3], 'k-', 'linewidth', 4);
            text(0, -1, 'stimulus', 'FontSize', 10)
            set(gca,'XLim', offset, 'Ytick', [], 'YLim', [-4 length(LL)+5]);
            xlabel('Time (s)')
            %plotting time of saccade (when saccade entered stim window)
            saccmade = mode(sacc_t-stimon_t);
            plot([saccmade saccmade], [-4 length(LL)+5], 'k-');
            text(saccmade, -3, ' \leftarrow saccade', 'FontSize', 10)
        end
        
        %PSTH
        for i = 1:size(unique_xy,1)
            subplot(2, size(unique_xy,1), i+size(unique_xy,1)); hold on;
            plot(bins, PSTH(i,:), 'k-', 'LineWidth', 2);
            plot([0 dur], [-2 -2], 'k-', 'linewidth', 4); %stim duration
            plot([saccmade saccmade], [0 10*ceil(max(PSTH(:)/10))], 'k-'); %saccade left fix window
            set(gca,'YLim',[-4 10*ceil(max(PSTH(:)/10))]);
            set(gca,'Xlim',offset);
        end
    end
end
