%SCStimcueAnalyses

%% Yasmine El-Shamayleh (5/21/2015)
% UNPACK SCSTIMCUE DATA
% PLOT RASTERS

%% NEX2 STRO BASICS

clear all; close all;

filename = 'F051915007.nex';

addpath(genpath('-GENERAL-')); 
stro = nex2stro(filename);

% synctimestr = 'FPOFF';
fpacq_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fpacq_t'));
fpon_t  = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fpon_t'));
fpoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fpoff_t'));
stimon_t  = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimon_t'));
stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimoff_t'));
targon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targon_t'));
targoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targoff_t'));
saccinit_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'saccinit_t'));
saccmade_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'saccmade_t'));
optstim = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'optstim'));            %THINKS STIM IS ALWAYS ON; EXCEPT BLANKS!
optstim_freq = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'optstim_freq'));  %NOT USED
targ_x = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targ_x'));
targ_y = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targ_y'));
cta = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'cta'));                    %??
visstim = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'visstim'));            %??
direction = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'direction'));        %??
correct= stro.trial(:,strcmp(stro.sum.trialFields(1,:),'correct'));             %??

%%

% synctimestr = 'SACCADE ON'; 
% sync_t = saccinit_t;
% offset = [-.3 0.45];              % pre and post time

synctimestr = 'FP OFF';
sync_t = fpoff_t;
offset = [-.3 0.7];              % pre and post time

Lspikechans = strncmp(stro.sum.rasterCells,'sig01',4);

% parse by 8 target locations
all_loc = [targ_x, targ_y];
loc = unique(all_loc,'rows');
loc = sortrows(loc);
%foo= loc./max(loc(:));
loc_order = [3 5 7 1 8 2 4 6];          %this orders the rasters by x/y coords
loc = loc(loc_order,:);


if (any(Lspikechans))
    for whichspike = find(Lspikechans)
        whichspike;
        spikes = stro.ras(:,whichspike);
        binwidth = .01;
        bins = offset(1):binwidth:offset(2);
        PSTH = zeros(length(loc),length(bins));
        
        fig1 = figure;
        annotation('textbox',[.4 .87 .2 .1],'string',['SC STIMCUE   ', filename(1:end-4)],'HorizontalAlignment','Left','FitBoxToText','on','EdgeColor','none')
        
        for i = 1:length(loc)             % 8 target directions
            for j = 1  %unique(optstim)   % just taking target trials (not blank)
                
                L = targ_x == loc(i,1) & targ_y == loc(i,2) & optstim == j;
                trlidxs = find(L);
                
                if i<5
                    subplot(6,3,i); hold on;
                else
                    subplot(6,3,i+1); hold on;
                end


                % PLOT RASTERS
                for counter = 1:sum(L)
                    trlidx = trlidxs(counter);
                    tmpspikes = spikes{trlidx}-sync_t(trlidx);
                    tmpspikes(tmpspikes < offset(1) | tmpspikes > offset(2)) = [];
                    nspikestot = length(tmpspikes);
                    plot([tmpspikes tmpspikes]',[zeros(nspikestot,1) .5*ones(nspikestot,1)]'+counter,'k-');
                    PSTH(i,:) = PSTH(i,:) + hist(tmpspikes, bins);
                end
                PSTH(i,:) = PSTH(i,:)./binwidth/sum(L);   
                set(gca,'XLim', offset,'Ytick',[],'YLim',[-1 sum(L)+2]);
                if (i == 1)
                    title(['Trials aligned to ',synctimestr]);
                end
            end
            
            % PLOT PSTHs
            for i = 1:length(loc)
                if i<5
                    subplot(6,3,i+9); hold on;
                else
                    subplot(6,3,i+1+9); hold on;
                end


                plot(bins,PSTH(i,:),'k-','LineWidth',2);
                set(gca,'YLim',[0 10*ceil(max(PSTH(:)/10))]);
                set(gca,'Xlim',offset);
            end
            
        end
    end
end

