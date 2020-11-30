%% YASMINE EL-SHAMAYLEH (11/8/2015)


%% LOAD

clear all; close all;
%fname = 'F100115002';
%fname = 'F100415002';
%fname = 'F101215001';
fname = 'F101215005';

nexfilepath = ['/Users/yasmine/Data-horwitzlab/nex/',fname,'.nex'];

%% PARSE

% Section 5
% Hacked together analysis of Fixstim_interfreq data
stro = nex2stro(nexfilepath);
offset = [-.1 .1];  % pre and post time wrt stimon/stimoff
stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimon_t'));
fpoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fpoff_t'));
fpacq_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fpacq_t'));
stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimoff_t'));
targon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targon_t'));
targoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targoff_t'));
opt_stimfreq = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'opt_stimfreq'));
targ_shown = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targ_shown'));
optfreq = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'elec_stimfreq'));
uniqfreqs = unique(optfreq);

stimon_t(isnan(stimon_t)) = fpacq_t(isnan(stimon_t))+nanmean(stimon_t-fpacq_t); % 

% There was a brief period of time when you couldn't determine if the modulator was used or not
% (the modulator header parameter didn't exist). Account for that below.
% [~, toks] = isvalidnexfilename(stro.sum.fileName);
% file_date = datenum(str2double([['20' toks{1}{4}] toks{1}(2) toks{1}(3)]));
% if file_date >= datenum([2015 5 7]) && file_date <= datenum([2015 5 14])
%     modulator = 1;
% elseif file_date < datenum([2015 5 7])
%     modulator = 0;
% else % after May 14, 2015
%     modulator = stro.sum.exptParams.modulator;
% end
modulator = stro.sum.exptParams.modulator;

% Hack below 
if (isnan(targ_shown(1)))
    targ_shown(:,1) = 0;
end
dur = mode(stimoff_t-stimon_t);

sync_t = stimon_t;  % Ecode for alignment
Lspikechans = strncmp(stro.sum.rasterCells,'sig0',4);

%% PLOT

fig1 =figure('position',[1400 1600 1200 1300]);
annotation('textbox',[.4 .89 .2 .1],'string',['FIXSTIM  CELL  ',fname],'HorizontalAlignment','Center','FitBoxToText','on','EdgeColor','none')


if (any(Lspikechans))
    for whichspike = find(Lspikechans)
        spikes = stro.ras(:,whichspike);
        binwidth = .002; % s
        bins = offset(1):binwidth:dur+offset(2);
        PSTH = zeros(1,length(bins));
        for j = uniqfreqs'
            PSTH = zeros(1,length(bins));
            %figure; subplot(2,1,1); hold on;
            subplot(length(uniqfreqs), 2, 2*(find(uniqfreqs==j))-1); hold on;
            L = optfreq == j;
            trlidxs = find(L);
            for counter = 1:sum(L)
                trlidx = trlidxs(counter);
                tmpspikes = spikes{trlidx}-sync_t(trlidx);
                tmpspikes(tmpspikes < offset(1) | tmpspikes > dur+offset(2)) = [];
                nspikestot = length(tmpspikes);
                plot([tmpspikes tmpspikes]',[zeros(nspikestot,1) .5*ones(nspikestot,1)]'+counter,'k-');
                PSTH = PSTH + hist(tmpspikes, bins);
            end
            PSTH = PSTH./(sum(L).*binwidth);
            if (j > 0) % Plotting the time course of optical stimulation
                if modulator
                    f = 0.5*1000^(1/255)^j;
                    t = linspace(bins(1),bins(end),5e4);
                    y = sin(2*pi*f*t-pi/2);
                    y(t<0) = -1;
                    y(t>dur) = -1;
                    y = y-.5; % for plotting below the spikes
                    plot(t,y,'k-','linewidth',2);
                else
                    secspercycle = 1/unique(j(j > 0));                    
                    transitions = 0:secspercycle/2:dur;
                    if (ceil(length(transitions)/2) ~= floor(length(transitions)/2))
                        transitions(end+1) = dur; %automatic shutoff
                    end
                    x = [transitions; transitions];
                    x = [offset(1); x(:); max(x(:))+offset(2)];
                    y = [repmat([0 1],1,length(transitions)/2) 0]*2-2;
                    y = [y;y];
                    if (length(x(:)) == length(y(:)))
                        plot(x,y(:)','k-','linewidth',2);
                    end
                end
            end
            set(gca,'XLim', [offset(1) dur+offset(2)],'Ytick',[],'YLim',[-4 sum(L)+5]);
            if modulator && j > 0
                title(['Frequency: ',num2str(round(f)),' Hz']);
            else
                title(['Frequency: ',num2str(j)]);
            end
            
            % PSTH
            %subplot(2,1,2); hold on;
            subplot(length(uniqfreqs), 2, 2*(find(uniqfreqs==j))-0); hold on;
            plot(bins,PSTH,'k-','LineWidth',2);
            set(gca,'YLim',[0 10*ceil(max(PSTH(:)/10))+1]);
            set(gca,'Xlim',[offset(1) dur+offset(2)]);
            xlabel('Time (s)','FontSize',12);
            ylabel('Response (sp/s)','FontSize',12);
        end
    end
end

%%

% FIG
fpath = (['/Volumes/hamlet/Users/yasmine/Data-horwitzlab/Figs/FixStim/','FixStim_',fname,'.eps']);
exportfig(fig1,fpath,'Color','cmyk');

clear all; close all;