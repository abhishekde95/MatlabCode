% IsoSamp/LMTF analyses 
%
% Contents:
% 1) Setting stuff up.
% 2) Plotting rasters and heatmap of stimulus modulation
% 3) Looking at raw spike count as a function of color direction/TF (useful
% for low frequency conditions in which < 1 cycle of the stimulus is shown)
% 4) Trial-by-trial calculation of signal-to-noise (several variants)
% 5) Comparing luminance SNR (as calculated above) to color SNR on a TF-by-TF basis
% 6) Estimating how many midget and parasol cell RFs overlap the stimulus.
% 6.5) Estimating how many midget and parasol cell RFs overlap the
% stimulus + estmating d' from RF overlaps.
% 7) Looking at efficiency of encoding (cone model to LGN cell)
% and decoding (population of LGN cells to behavior) as a function of
% TF/color dir.
% 8) Cross-correlations functions (in case two neurons were recorded
% simultaneously).
% 9) Comparing contrasts used in an IsoSamp experiment to the most recent
% psychophysical model fits (LMTF.mat)
% 10) Trying spike train classification with the Daniel Reich & Jonathan
% Victor spike distance metric
% 11) Comparing ideal observer performance (of cones and LGN cells),
% allowing the counting window to move around
% 12) Simulation of cone current and photon absorption ideal observers for
% a range of contrasts. Cone current SNR should increase linearly with
% contrast and photon absorption SNR should increase with the square root
% of contrast. (Trying to figure out why photon absorption ideal observer
% is always ~1.6 more sensitive than the cone current ideal observer.)

%% Section 1: Setting stuff up.
stro = nex2stro;
if stro.sum.paradigmID ~= 107
    error('Not an IsoSamp file');
end
rfx = stro.sum.exptParams.rf_x/10;
rfy = stro.sum.exptParams.rf_y/10;
Lcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_l'));
Mcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_m'));
TF = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'tf'));
stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimon_t'));
stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimoff_t'));
fpoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'rew_t'));

dur = mean(stimoff_t-stimon_t);
spikeidxs = strncmp(stro.sum.rasterCells(1,:),'sig',3);
%spikes = stro.ras(:,spikeidx);
uniquestim = sortrows(unique([Lcc Mcc TF],'rows'),3); % sorting by TF
totaltrials = size(stro.trial,1);

%% Section 2: Stimulus and rasters
bins = linspace(0,dur,100);
data = zeros(size(uniquestim,1), length(bins));
temporalenvelope = ones(size(bins));
temporalenvelope(1:round((1/4)*length(bins))) = linspace(0,1,round((1/4)*length(bins)));
temporalenvelope(end:-1:round((3/4)*length(bins))+1) = linspace(0,1,round((1/4)*length(bins)));

for i = 1:size(uniquestim,1)
    contrast = sqrt(uniquestim(i,1).^2+uniquestim(i,2).^2);
    data(i,:) = contrast*temporalenvelope.*sin(2*pi*uniquestim(i,3)*bins);
end

for spikeidx = find(spikeidxs)
    spikes = stro.ras(:,spikeidx);
    for stimtype = 1:2 % non-opponent, opponent
        if (stimtype == 1)
            Lstimtype = sign(uniquestim(:,1)) == sign(uniquestim(:,2));
        else
            Lstimtype = sign(uniquestim(:,1)) ~= sign(uniquestim(:,2));
        end
        if (any(Lstimtype))
            
            figure;
            subplot(1,2,2);
            im = 255*(data - min(data(:)))/(max(data(:)) - min(data(:)));
            image(flipud(im(Lstimtype,:)));
            xtick = 0:.2:dur;
            set(gca,'XTick',interp1([bins(1) bins(end)],[1 length(bins)],xtick));
            set(gca,'Xticklabel',xtick);
            
            ytick = round(uniquestim(1,3)):5:round(uniquestim(end,3));
            set(gca,'YTick',interp1([ytick(1) ytick(end)],[1 size(data,1)],ytick));
            set(gca,'Yticklabel',fliplr(ytick));
            colormap(jet(255));
            
            % Plotting rasters
            offset = [-.1 .1];  % pre and post time wrt stimon/stimoff
            subplot(1,2,1); hold on; counter = 0;
            for j = find(Lstimtype)'
                L = Lcc == uniquestim(j,1) & Mcc == uniquestim(j,2) & TF == uniquestim(j,3);
                for i = find(L)'
                    tmpspikes = spikes{i}-stimon_t(i);
                    tmpspikes(tmpspikes < offset(1) | tmpspikes > dur+offset(2)) = [];
                    nspikestot = length(tmpspikes);
                    plot([tmpspikes tmpspikes]',[zeros(nspikestot,1) 1*ones(nspikestot,1)]'+counter,'k-','linewidth',1);
                    counter = counter + 1;
                end
            end
            set(gca,'Xlim',[0+offset(1) dur+offset(2)],'Ylim',[0 counter],'Ytick',[],'Xtick',[0 .2 .4 .6 .8],'Box','on');
            if stimtype ==1
                title('LUM');
            else
                title('RG');
            end
        end
        set(gcf,'name',stro.sum.rasterCells{spikeidx})
    end
end

%% Section 3
% ----------------------------
% Heat map based on spike counts
% Only makes sense for low frequencies

for spikeidx = find(spikeidxs)
    % Getting raw spike rates
    spikerates = [];
    for i = 1:size(stro.trial,1)
        spiketimes = stro.ras{i,spikeidx};
        if isempty(spiketimes)
            nspikes = 0;
        else
            nspikes = sum(spiketimes > stimon_t(i) & spiketimes < stimon_t(i)+dur);
        end
        spikerates = [spikerates; nspikes./(stimoff_t(i)-stimon_t(i))];
    end
    
    % Getting normalized spikerates per condition
    normfr = nan(size(uniquestim,1),1);
    for j = 1:size(uniquestim,1)
        L = Lcc == uniquestim(j,1) & Mcc == uniquestim(j,2) & TF == uniquestim(j,3);
        normfr(j)= mean(spikerates(L));
    end
    normfr = normfr-min(normfr);
    normfr = normfr./max(normfr);
    
    % Plotting
    cmap = jet(255);
    figure; axes; hold on;
    Lblank = Lcc == 0 & Mcc == 0 & TF == 0;
    
    for j = 1:size(uniquestim,1)
        L = Lcc == uniquestim(j,1) & Mcc == uniquestim(j,2) & TF == uniquestim(j,3);
        if (uniquestim(j,3) < 5)
            if (uniquestim(j,3) == 0) % blank trials
                TF_for_plotting = 1;
            else
                TF_for_plotting = uniquestim(j,3);
            end
            h = plot3(uniquestim(j,1),uniquestim(j,2),TF_for_plotting,'ko','MarkerSize',20*normfr(j)+10);
            set(h,'MarkerEdgeColor','none','MarkerFaceColor',cmap(round((size(cmap,1)-1)*normfr(j)+1),:));
            
            % T-tests
            L = Lcc == uniquestim(j,1) & Mcc == uniquestim(j,2) & TF == uniquestim(j,3);
            [h,p] = ttest2(spikerates(Lblank),spikerates(L));
            if p < 0.05
                h = plot3(uniquestim(j,1),uniquestim(j,2),TF_for_plotting,'y*','MarkerSize',20*normfr(j)+10);
            end
        end
    end
    set(gca,'View',[50 20], 'Zscale','log');
    xlabel('L-cone contrast'); ylabel('M-cone contrast'); zlabel('TF');
    title('Symbol size is spike rate');
    set(gcf,'name',stro.sum.rasterCells{spikeidx})
end
%%
% Section 4
% Individual tests for each trial. Null hypothesis depends on
% the number of spikes.
% See "StatsStuff Section 7" for Poisson spike trains projected onto basis
% vectors.
% This code should be replaced by "IsoSampGetDPrime.m".

for spikeidx = find(spikeidxs)
    signaloffsets = [0.05 dur];
    spiketimes = {}; % First, getting all the spike times
    nspikes = [];
    spikes = stro.ras(:,spikeidx);
    for i = 1:totaltrials
        current_TF = TF(i);
        tmpspikes = spikes{i}-stimon_t(i);
        spiketimes{i} = tmpspikes(tmpspikes > signaloffsets(1) & tmpspikes < signaloffsets(end));
        spiketimes{i} = spiketimes{i} - signaloffsets(1); % So time = 0 is begining of counting window
        nspikes(i) = length(spiketimes{i});
    end
    Lblank = Lcc==0 & Mcc==0 & TF==0;
    blank_trl_idxs = find(Lblank);
    
    % Getting dot products
    signal_data = nan(totaltrials,2); % dot products onto basis vectors
    for i = 1:totaltrials
        current_TF = TF(i);
        signal_data(i,:) = [sum(cos(2*pi*current_TF.*spiketimes{i})) sum(sin(2*pi*current_TF.*spiketimes{i}))];
        nspikes(i) = length(spiketimes{i});
    end
    
    ps = [];
    T2s = [];
    for i = 1:totaltrials
        current_TF = TF(i);
        mu_noise = nspikes(i)*[sin(2*pi*current_TF*dur)+0 -cos(2*pi*current_TF*dur)+1]/(2*pi*current_TF*dur); % Integral of cos(ax) is sin(ax)/a. Integral of sin(ax) is -cos(ax)/a
        crossprods(1,1) = (1/dur)*((dur/2)+sin(4*pi*current_TF*dur)/(8*pi*current_TF)) - (mu_noise(1)/nspikes(i))^2;
        crossprods(1,2) = (1/dur)*(-1*cos(4*pi*current_TF*dur)+1)/(8*pi*current_TF)-mu_noise(1)*mu_noise(2)/nspikes(i)^2;
        crossprods(2,1) = crossprods(1,2);
        crossprods(2,2) = (1/dur)*((dur/2)-sin(4*pi*current_TF*dur)/(8*pi*current_TF)) - (mu_noise(2)/nspikes(i))^2;
        S2 = crossprods*nspikes(i);
        if isnan(rcond(S2))
            T2s(i) = 0;
            ps(i) = 1;
        else
            T2s(i) = (signal_data(i,:)-mu_noise)*inv(S2)*(signal_data(i,:)-mu_noise)';
            ps(i) = 1-chi2cdf(T2s(i),2);
        end
    end
    
    % Debugging
    %figure; axes; hold on;
    %hist(T2s,linspace(0,20,100))
    %plot(linspace(0,20,100),chi2pdf(linspace(0,20,100),2)*length(T2s)/5,'y-')
    
    % Getting an empirical distribution of T2 values (per TF) for the zero
    % contrast trials.
    empirical_noise_T2 = [];
    for i = 1:size(uniquestim,1)
        current_TF = uniquestim(i,3);
        tmp = [];
        for j = find(Lblank)'
            mu_noise = nspikes(j)*[sin(2*pi*current_TF*dur)+0 -cos(2*pi*current_TF*dur)+1]/(2*pi*current_TF*dur); % Integral of cos(ax) is sin(ax)/a. Integral of sin(ax) is -cos(ax)/a
            crossprods(1,1) = (1/dur)*((dur/2)+sin(4*pi*current_TF*dur)/(8*pi*current_TF)) - (mu_noise(1)/nspikes(j))^2;
            crossprods(1,2) = (1/dur)*(-1*cos(4*pi*current_TF*dur)+1)/(8*pi*current_TF)-mu_noise(1)*mu_noise(2)/nspikes(j)^2;
            crossprods(2,1) = crossprods(1,2);
            crossprods(2,2) = (1/dur)*((dur/2)-sin(4*pi*current_TF*dur)/(8*pi*current_TF)) - (mu_noise(2)/nspikes(j))^2;
            S2 = crossprods*nspikes(j);
            if isnan(rcond(S2))
                tmp = [tmp; 0];
            else
                noiseprojections = [sum(cos(2*pi*current_TF.*spiketimes{j})) sum(sin(2*pi*current_TF.*spiketimes{j}))];
                tmp = [tmp; (noiseprojections-mu_noise)*inv(S2)*(noiseprojections-mu_noise)'];
            end
        end
        empirical_noise_T2{i} = tmp;
    end
    
    % Plotting rasters
    alpha = .1;
    offset = [-.1 .1];  % pre and post time wrt stimon/stimoff
    figure;
    for stimtype = 1:2 % non-opponent, opponent
        if (stimtype == 1)
            Lstimtype = sign(uniquestim(:,1)) == sign(uniquestim(:,2));
        else
            Lstimtype = sign(uniquestim(:,1)) ~= sign(uniquestim(:,2));
        end
        if (any(Lstimtype))
            % Plotting rasters
            subplot(1,2,stimtype); hold on; counter = 0;
            for j = find(Lstimtype)'
                L = Lcc == uniquestim(j,1) & Mcc == uniquestim(j,2) & TF == uniquestim(j,3);
                for i = find(L)'
                    tmpspikes = spikes{i}-stimon_t(i);
                    tmpspikes(tmpspikes < offset(1) | tmpspikes > dur+offset(2)) = [];
                    nspikestot = length(tmpspikes);
                    h = plot([tmpspikes tmpspikes]',[zeros(nspikestot,1) 1*ones(nspikestot,1)]'+counter,'-','linewidth',1,'color',[.5 .5 .5]);
                    sig = ps(i) < alpha;
                    if sig
                        set(h,'Color','black','Linewidth',2);
                    end
                    counter = counter + 1;
                end
            end
            set(gca,'Xlim',[0 dur],'Ylim',[0 counter],'Ytick',[],'Xtick',[0 .2 .4 .6 .8],'Box','on');
            if stimtype ==1
                title('LUM');
            else
                title('RG');
            end
        end
    end
    
    % "Signal to noise" as a function of condition
    snr = {};
    for j = 1:size(uniquestim,1)
        L = Lcc == uniquestim(j,1) & Mcc == uniquestim(j,2) & TF == uniquestim(j,3);
        snr{j} = log10(T2s(L));
        snr{j}(isinf(snr{j})) = nan; % T2s = 0 converted to snr = nan. Ignored by nanmean, etc.
    end
    set(gcf,'name',stro.sum.rasterCells{spikeidx})

    % Plotting raw "SNRs"
    figure; subplot(2,2,1); hold on;
    for j = 1:size(uniquestim,1)
        if (sign(uniquestim(j,1)) ~= sign(uniquestim(j,2)))
            plot(uniquestim(j,3),snr{j},'r.')
            h = errorbar(uniquestim(j,3),nanmean(snr{j}),1.96*nanstd(snr{j})./sum(~isnan(snr{j})));
            set(h,'LineWidth',2,'Color','red');
            hsig = plot(uniquestim(j,3),nanmean(snr{j}),'ro','MarkerFaceColor','red');
        else
            plot(uniquestim(j,3),snr{j},'k.')
            h = errorbar(uniquestim(j,3),nanmean(snr{j}),1.96*nanstd(snr{j})./sum(~isnan(snr{j})));
            set(h,'LineWidth',2,'Color','black');
            hsig = plot(uniquestim(j,3),nanmean(snr{j}),'ko','MarkerFaceColor','black');
        end
        h = plot(uniquestim(j,3)+.1,nanmean(log10(empirical_noise_T2{j})),'bs','MarkerFaceColor','blue','MarkerSize',10);
        plot(uniquestim(j,3)+.1,log10(empirical_noise_T2{j}),'b.')
        [~,p] = ttest2(log10(empirical_noise_T2{j}), snr{j}');
        if p < 0.05
            set(hsig,'MarkerFaceColor','yellow');
        end
    end
    set(gca,'Xscale','log','Yscale','linear');
    ylabel('log10(T^2)');
    xlabel('TF(Hz)');
    
    % ROCs
    subplot(2,2,2); hold on;
    for j = 1:size(uniquestim,1)
        if (sign(uniquestim(j,1)) ~= sign(uniquestim(j,2)))
            plot(uniquestim(j,3),roc(log10(empirical_noise_T2{j}),snr{j}),'ro','MarkerFaceColor','red');
        else
            plot(uniquestim(j,3),roc(log10(empirical_noise_T2{j}),snr{j}),'ko','MarkerFaceColor','black');
        end
    end
    plot([1 40],[.5 .5],'k--');
    ylabel('AUROC');xlabel('TF (Hz)'); set(gca,'Xscale','log')
    
    % D-primes (assuming different noise distributions for each TF)
    subplot(2,2,3); hold on;
    d_prime=[];
    for j = 1:size(uniquestim,1)
        signal = snr{j};
        noise = log10(empirical_noise_T2{j});
        noise(isinf(noise)) = nan;
        n_signal = sum(~isnan(signal));
        n_noise = sum(~isnan(noise));
        pooled_var = ((n_signal-1)*nanvar(signal)+(n_noise-1)*nanvar(noise))/(n_signal+n_noise-2);
        d_prime(j) = (nanmean(signal)-nanmean(noise))/sqrt(pooled_var);
        if (sign(uniquestim(j,1)) ~= sign(uniquestim(j,2)))
            plot(uniquestim(j,3),d_prime(j),'ro','MarkerFaceColor','red');
        else
            plot(uniquestim(j,3),d_prime(j),'ko','MarkerFaceColor','black');
        end
    end
    plot([1 40],[0 0],'k:'); % Complete insensitivity
    plot([1 40],[1.27 1.27],'k--'); % Psychophysical d'
    ylabel('d''');xlabel('TF (Hz)'); set(gca,'Xscale','log')
    title('individual noise distns');
    
    % D-primes (same noise distribution for each TF)
    noise_means = [];
    noise_vars = [];
    for j = 1:size(uniquestim,1)
        tmp = empirical_noise_T2{j};
        if all(tmp == 0)
            continue
        else
            tmp(tmp == 0) = [];
            noise_means = [noise_means; nanmean(log10(tmp))];
            noise_vars = [noise_vars; nanvar(log10(tmp))];
        end
    end
    noise_mean = mean(noise_means);
    noise_std = sqrt(mean(noise_vars));
    
    subplot(2,2,4); hold on;
    d_prime=[];
    for j = 1:size(uniquestim,1)
        signal = snr{j};
        d_prime(j) = (nanmean(snr{j})-noise_mean)/noise_std;
        if (sign(uniquestim(j,1)) ~= sign(uniquestim(j,2)))
            plot(uniquestim(j,3),d_prime(j),'ro','MarkerFaceColor','red');
        else
            plot(uniquestim(j,3),d_prime(j),'ko','MarkerFaceColor','black');
        end
    end
    plot([1 40],[0 0],'k:'); % Complete insensitivity
    plot([1 40],[1.27 1.27],'k--'); % Psychophysical d'
    ylabel('d''');xlabel('TF (Hz)'); set(gca,'Xscale','log')
    title('single noise distn');
    set(gcf,'name',stro.sum.rasterCells{spikeidx})
end
%%
% Section 4.5
% Calling IsoSampGetDPrime with fake data to test the sensitivity/bias of
% the d' prime estimate as a function of temporal frequency and phase. What
% happens to neurons that are suppressed by the stimulus (haven't answered
% this question yet).

TESTING = 0;
if TESTING == 1
    niter = 100;
    phi = pi/4;
    figure; axes; hold on;
    for iter = 1:niter
        tmpstro = stro;
        for i = 1:size(tmpstro.trial,1)
            nspikes = length(tmpstro.ras{i});
            if nspikes == 0
                tmpstro.ras{i} = [];
                continue
            elseif length(tmpstro.ras{i}) == 1
                continue
            else
                fakespiketimes = unifrnd(tmpstro.ras{i}(1),tmpstro.ras{i}(end),length(tmpstro.ras{i}),1);
                %fakespiketimes = linspace(tmpstro.ras{i}(1),tmpstro.ras{i}(end),length(tmpstro.ras{i}));
                
                % What is the integral of (sin(ax+b)+1)/2?
                % It's (x-cos(ax+b)/a)/2
                if TF(i) > 0
                    c = TF(i).*2.*pi;
                    z = (fakespiketimes-cos(fakespiketimes.*c+phi)/c);
                    tmpstro.ras{i} = sort(z)';
                else
                    tmpstro.ras{i} = fakespiketimes;
                end
            end
        end
        [tmpstim, tmpd] = IsoSampGetDPrime(tmpstro);
        Lrg = sign(tmpstim(:,1)) ~= sign(tmpstim(:,2));
        Llum = sign(tmpstim(:,1)) == sign(tmpstim(:,2)) & sign(tmpstim(:,3)) ~= 0;
        
        plot(tmpstim(Lrg,3),tmpd(Lrg),'ro-');
        plot(tmpstim(Llum,3),tmpd(Llum),'ko-');
        set(gca,'Xscale','log');
        drawnow;
    end
    
end

%%
% Section 5
% Comparing luminance to color on a TF by TF basis

[~, d_prime] = IsoSampGetDPrime(stro);

sigcompare = [];
for currentTF = unique(TF)'
    L = TF == currentTF;
    LuniqueTF = uniquestim(:,3) == currentTF;
    Lchrom = LuniqueTF & sign(uniquestim(:,1)) ~= sign(uniquestim(:,2));
    Llum = LuniqueTF & sign(uniquestim(:,1)) == sign(uniquestim(:,2));
    if sum(Lchrom) == 1 & sum(Llum) == 1
        sigcompare = [sigcompare; currentTF d_prime(Lchrom) d_prime(Llum)];
    end
end

cmap = colormap(jet(64));
figure; axes; hold on;
for i = 1:size(sigcompare,1)
    h = plot(sigcompare(i,2),sigcompare(i,3),'o','MarkerEdgeColor','none','MarkerSize',10);
    cidx = (log10(sigcompare(i,1))-log10(min(sigcompare(:,1))))./(log10(max(sigcompare(:,1)))-log10(min(sigcompare(:,1))));
    set(h,'MarkerFaceColor',cmap(floor(cidx*63+1),:));
end
xlabel('Chrom d-prime');
ylabel('Lum d-prime');

ylims = ylim;
xlims = xlim;
lims = [max(xlims(1),ylims(1)),min(xlims(2),ylims(2))];
plot(lims,lims,'k:'); % Identity line


%% 
% Section 6
% Estimating the number of midget and parasol retinal ganglion cell RFs
% that overlap the stimulus (crudely) and estimating the d' of each one if
% the exactly matched the monkeys d'=1.27.

% Key reference: Dacey and Petersen (1992) Dendritic field size and
% morphology of midget and parasol ganglion cells of the human retina.
MMPERDEG = 0.223; % mm/deg (Perry and Cowey 1985)
DEGPERMM = 1/MMPERDEG; % deg/mm
sigma = unique(stro.trial(:,strcmp(stro.sum.trialFields(1,:),'sigma'))); % DVA
sigmas_n = unique(stro.trial(:,strcmp(stro.sum.trialFields(1,:),'sigmas_n')));
% Converting RF location to mm on retina
rfy_mm = rfy*MMPERDEG;
rfx_mm = rfx*MMPERDEG;
rf_theta = atan2(rfy_mm,rfx_mm);
rf_r_mm = sqrt(rfx_mm^2 + rfy_mm^2);
% Perry and Cowey method <--- not needed. Linear approx. is fine (better?)
% with central 40�. Watson 2014.
% deg = 0.1+4.21(mm)+0.038(mm^2) 
%rf_r_deg = sqrt(rfx^2+rfy^2);
%rf_r_mm = (5/38)*(sqrt(1520*rf_r_deg + 177089) - 421); % gives similar answers

% For starters, how many RFs fit within 'n' stds 
% y is dendritic field diameter in microns
% THESE EQUATIONS ARE FOR HUMAN DATA (which, for parasol cells, is different from monkey)
% Dacey and Petersen (1992)
RF_diam_midget_mm = 8.64*rf_r_mm^1.04/1000; % midgets
%RF_diam_parasol_mm = 70.2*rf_r_mm^0.65/1000; % human parasol cells
RF_diam_midget_deg = RF_diam_midget_mm*DEGPERMM;
n_midgets = (2*sigma/RF_diam_midget_deg)^2;
pred_d_prime_midget = 1.27/sqrt(n_midgets);

% Now using the Watanabe and Rodieck (1989) equations
% for monkey parasol cells. This seems like a gross overestimate at the RF
% locations I'm interested in, near the fovea.
%RF_diam_parasol_mm = 0.0206*rf_r_mm+0.0518;
%RF_diam_parasol_deg = RF_diam_parasol_mm*DEGPERMM;
%n_parasol = ((2*sigma)/RF_diam_parasol_deg)^2

% A rough estimate of parasol RF sizes based on figure 6A of Perry et al.
% (1984)
RF_diam_parasol_mm = .03; 
RF_diam_parasol_deg = RF_diam_parasol_mm*DEGPERMM;
n_parasol = ((2*sigma)/RF_diam_parasol_deg)^2;
pred_d_prime_parasol = 1.27/sqrt(n_parasol);

% Plotting D-primes with predictions from RGC counting
figure; axes; hold on;
for j = 1:size(uniquestim,1)
    signal = snr{j};
    noise = log10(empirical_noise_T2{j});
    noise(isinf(noise)) = nan;
    n_signal = sum(~isnan(signal));
    n_noise = sum(~isnan(noise));
    pooled_var = ((n_signal-1)*nanvar(signal)+(n_noise-1)*nanvar(noise))/(n_signal+n_noise-2);
    d_prime(j) = (nanmean(signal)-nanmean(noise))/sqrt(pooled_var);
    if (sign(uniquestim(j,1)) ~= sign(uniquestim(j,2)))
        plot(uniquestim(j,3),d_prime(j),'ro','MarkerFaceColor','red');
    else
        plot(uniquestim(j,3),d_prime(j),'ko','MarkerFaceColor','black');
    end
end
h = [];
plot([1 40],[0 0],'k:','linewidth',1); % Complete insensitivity
h(1) = plot([1 max(TF)],[1.27 1.27],'b-','linewidth',1); % Psychophysical d'
h(2) = plot([1 max(TF)],[1 1]*pred_d_prime_midget,'r--','linewidth',1); 
h(3) = plot([1 max(TF)],[1 1]*pred_d_prime_parasol,'k--','linewidth',1); 

legend(h,{'Monkey','Single indep. midget','Single indep. parasol'},'location','southeast');
ylabel('d''');xlabel('TF (Hz)'); set(gca,'Xscale','log','Xlim',[min(TF) max(TF)]);

%%
% Section 6.5
% A slightly more involved "simulation" using a mosaic of midget ganglion
% cells overlapping the stimulus. First doing it in discrete space.
% According to Chichilnisky, ON and OFF midgets form anatomical mosaics
% with dendrites that do not overlap and this corresponds, physiologically,
% to RFs that abut at the 1 STD point of a Gaussian fit.

% Things to improve:
% 1) Yoke grain of simulated screen to cone mosaic or pixels on
% real screen.
% 2) Average nasal and temporal retinae RF sizes?
sigma_gabor = unique(stro.trial(:,strcmp(stro.sum.trialFields(1,:),'sigma'))); % DVA
sigmas_n = unique(stro.trial(:,strcmp(stro.sum.trialFields(1,:),'sigmas_n')));
[x_deg,y_deg] = meshgrid(linspace(-sigma_gabor*2*sigmas_n,sigma_gabor*2*sigmas_n,40));
stim = normpdf(x_deg,0,sigma_gabor).*normpdf(y_deg,0,sigma_gabor);

% Plotting the stimulus as it appears in DVA
figure; axes; hold on;
imagesc(x_deg(1,:)',y_deg(:,1),stim);
axis equal;
axis tight;
ylabel('DVA');
xlabel('DVA');
RF_diam_midget_mm = 8.64*rf_r_mm^1.04/1000; % midget dendritic field from Dacey and Petersen.
RF_diam_midget_deg = RF_diam_midget_mm*DEGPERMM;

% Estimating parasol RF size
% perrydata = xlsread ('PerryEtAl1984ParasolData'); % x = eccentricity (mm), y = dendritic field diam. (microns)
% b = regress(perrydata(:,2)*DEGPERMM/1000,[ones(size(perrydata,1),1) perrydata(:,1)*DEGPERMM]); % intercept and slope in degrees (DF diam) as a fucntion of degrees (eccentrcity)
RF_diam_parasol_deg = 0.0673+0.0165*rf_r_mm*DEGPERMM; % Regression coefficient from Perry et al. 1984

if (RF_diam_midget_deg < 0.01)
    disp('Too many midget cells for brute force simulation');
    return
end
% Plotting the central midget and parasol RFs
%plot(RF_diam_midget_deg/2*cos(linspace(0,2*pi,100)), RF_diam_midget_deg/2*sin(linspace(0,2*pi,100)),'r-')
%plot(RF_diam_parasol_deg/2*cos(linspace(0,2*pi,100)), RF_diam_parasol_deg/2*sin(linspace(0,2*pi,100)),'k-')

for i = 1:2
    if i == 1
        RF_diam_deg = RF_diam_midget_deg;
        celltypestr = 'midget';
        linetype = 'r-';
    else
        RF_diam_deg = RF_diam_parasol_deg;
        celltypestr = 'parasol';
        linetype = 'k';
    end
    % Below code adapted from http://matlabdatamining.blogspot.com/2008/04/generating-hexagonal-grids-for-fun-and.html
    % Generate hexagonal grid
    Rad3Over2 = sqrt(3)/2;
    x_centers = [x_deg(1,1):RF_diam_deg:x_deg(1,end)];
    closest_to_zero = find(abs(x_centers) == min(abs(x_centers)),1);
    x_centers = x_centers-x_centers(closest_to_zero);
    x_centers_mat = repmat(x_centers,size(x_centers,2),1);
    y_centers = x_centers;
    x_centers_mat = x_centers_mat*Rad3Over2;
    y_centers_mat = repmat(y_centers',1,size(y_centers,2));
    y_centers_mat(:,rem(find(y_centers == 0),2)+1:2:end) = y_centers_mat(:,rem(find(y_centers == 0),2)+1:2:end)+.5*RF_diam_deg;
    tmp = RF_diam_deg/2*[cos(linspace(0,2*pi,200))', sin(linspace(0,2*pi,200))'];
    weights = []; 
    for i = 1:numel(x_centers_mat)
        plot(x_centers_mat(i)+tmp(:,1),y_centers_mat(i)+tmp(:,2),linetype);
        % Formulas used below are from https://www.cs.nyu.edu/~roweis/notes/gaussid.pdf
        % See below for sanity check. Ignoring terms that are common to all
        % RF locations. Taking advantage of fact that all RFs have the same
        % sigma and mean of stimulus is (0,0).
        b = [x_centers_mat(i); y_centers_mat(i)];
        weights(i) = exp(-.5*(b'*b)/(RF_diam_midget_deg^2+sigma_gabor^2)); %https://math.stackexchange.com/questions/1260015/normalizing-factor-for-product-of-gaussian-densities-interpretation-with-bayes
    end
    weights = weights./max(weights);
%     % Sanity checking normalization
%     [x,y] = meshgrid(linspace(-10,10,200),linspace(-10,10,200));
%     tmpdata = [];
%     for tmpi = 1:100
%         mu1 = [0 0]'; mu2 = [normrnd(0,2) normrnd(0,2)]'; sigma1 = 2; sigma2 = 1;
%         y1 = normpdf(x,mu1(1),sigma1).*normpdf(y,mu1(2),sigma1);
%         y2 = normpdf(x,mu2(1),sigma2).*normpdf(y,mu2(2),sigma2);
%         sigma2_12 = 1/((1/sigma1^2)+(1/sigma2^2));
%         mu_product = (sigma2_12/sigma2^2)*mu2;
%         %mu_12 = ((mu1/sigma1^2)+(mu2/sigma2^2))* sigma2_12;
% %        tmpdata = [tmpdata; max(y1.*y2) exp(-.5*((mu2^2/sigma2^2)-mu_12^2/sigma2_12))];
%         tmpdata = [tmpdata; max(y1(:).*y2(:)) exp(-.5*(mu2'*mu2)/(sigma1^2+sigma2^2))];
%     end
%     figure; axes; hold on;
%     plot(tmpdata(:,1),tmpdata(:,2),'.'); % "normalization constant" is linearly related to peak of product pdf 
    % figure; axes; hold on;
    % plot3(x_centers_mat(:),y_centers_mat(:),weights,'o');
    % axis equal
    % axis tight
    
    % Assuming that I'm recording from the most sensitive neuron (centered on
    % the Gabor). For now, collapsing across conditions
    
    L_opp = sign(uniquestim(:,1))~=sign(uniquestim(:,2));
    L_nopp = sign(uniquestim(:,1))==sign(uniquestim(:,2)) & uniquestim(:,1) ~= 0;
    
    pred_chrom_d_prime = nanmedian(d_prime(L_opp))*(sum(weights.^2)./sqrt(sum(weights.^2))); % square root of sum of squared weights
    pred_lum_d_prime = nanmedian(d_prime(L_nopp))*(sum(weights.^2)./sqrt(sum(weights.^2)));
    
    disp(['Assuming this is a ',celltypestr,' cell, the predicted L-M d'' is ',num2str(pred_chrom_d_prime)]);
    disp(['Assuming this is a ',celltypestr,' cell, the predicted L+M d'' is ',num2str(pred_lum_d_prime)]);
end

% Is the above equation correct for the d' of a weighted sum? Doing some
% simulations as a sanity check:
TESTING = 0;
if TESTING
    tmpdata =[];
    signalquality = [.1 .3 .1 10]; % number of elements is number of neurons
    signal = 1;
    for i = 1:100
        simulated_spike_counts = normrnd(0,1,100,length(signalquality));
        simulated_spike_counts = simulated_spike_counts+repmat(signalquality,size(simulated_spike_counts,1),1); % adding signal
        simulated_weighted_spike_counts = simulated_spike_counts.*repmat(signalquality,size(simulated_spike_counts,1),1); % applying weights
        num = sum(mean(simulated_weighted_spike_counts));
        den = sqrt(sum(var(simulated_weighted_spike_counts)));
        tmpdata(i) = num/den;
    end
    figure; axes; hold on;
    plot(tmpdata)
    plot([1, length(tmpdata)],[1 1]*sum(signalquality.^2)./sqrt(sum(signalquality.^2)),'k-');
end
%%
% Section 7
% Encoding and decoding efficiency using cone model and behavior as the
% benchmarks for the input and output of the system
% Using stuff that was calculated in previous section (6.5).
for spikeidx = find(spikeidxs)
    [uniquestim, dprime] = IsoSampGetDPrime(stro,2,spikeidx);
    
    fundamentals = load('T_cones_smj10');
    MMPERDEG = 0.223; % mm/deg
    DEGPERMM = 1/MMPERDEG; % deg/mm
    
    % Setting up for cone model
    params.runType = 'isosamp';
    params.obsMethod = 'obsMethod_filteredWtFxn';
    params.impulseResponse = 'rieke';
    params.DTV1_fname = [];
    params.DTNT_fname = [];
    params.unitTest = false;
    params.eqMosaic = false;
    params.saveDir = '/Users/greghorwitz/Documents/MATLAB';
    params.notes = 'IsoSamp test';
    params.parallelOperations = false;
    params.eyeType = 'monkey';
    params.eyeNumber = 1;
    params.coneSampRate = 2400;
    params.flatPowerSpect = false;
    params.enableScones = false;
    params.sacamp_deg = 0;
    params.sacdur_s = 0;
    params.stro = stro;
    %cal.gamma = params.stro.sum.exptParams.gamma_table; % Only used for "dtv1". Can ignore.
    spds = params.stro.sum.exptParams.mon_spd;
    spds = SplineSpd([380:4:780]',reshape(spds,length(spds)/3,3),[380:5:780]');
    cal.monSpect = spds(:);
    M = fundamentals.T_cones_smj10*spds;
    cal.Mmtx = M(:);
    cal.frameRate = params.stro.sum.exptParams.framerate;
    cal.bkgndrgb = params.stro.sum.exptParams.bkgndrgb'; % intensities
    cal.fname = 'test';
    cal.monSpectWavelengths = linspace(380,780,101);
    cal.pixperdeg = params.stro.sum.exptParams.pixperdeg;
    params.monCalFile = cal;
    % End of cone model set up
    
    rfx = stro.sum.exptParams.rf_x/10;
    rfy = stro.sum.exptParams.rf_y/10;
    rf_r_deg = sqrt(rfx^2+rfy^2);
    rf_r_mm = (5/38)*(sqrt(1520*rf_r_deg + 177089) - 421);
    RF_diam_midget_mm = 8.64*rf_r_mm^1.04/1000; % midgets
    RF_diam_midget_deg = RF_diam_midget_mm*DEGPERMM;
    RF_diam_parasol_deg = 0.0673+0.0165*rf_r_deg; % Regression coefficients from Perry et al. 1984
    
    % Should make this a query to the user
    choice = questdlg('What type of cell to simulate?','CellMenu','Magno','Parvo','Cancel','Cancel');
    if strcmp(choice,'Parvo')
        params.gab.sd = RF_diam_midget_deg;
    elseif strcmp(choice,'Magno')
        params.gab.sd = RF_diam_parasol_deg;
    else
        return
    end
    [gab, cones, mon, idlob, params] = DTcones_gh(params,0);
    
    conemodeldata = nan*ones(size(uniquestim,1),1);
    for i = 1:size(idlob.analyticMean,1) % looping over color direction
        for j = 1:size(idlob.analyticMean(i,:),2) % looping over contrast/TF
            if ~isempty(idlob.analyticMean{i,j})
                tmp_lm_mu = idlob.analyticMean{i,j}([1 2]);
                tmp_lm_var = idlob.analyticVar{i,j}([1 2]);
                tf = gab.driftRates{i}(j);
                L = uniquestim(:,3) == gab.driftRates{i}(j) & sign(gab.colorDirs(i,1)) == sign(uniquestim(:,1));
                if (sum(L) == 0)
                    error('cannot find condition');
                end
                conemodeldata(L) = [sum(abs(tmp_lm_mu))/sqrt(sum(tmp_lm_var))];
            end
        end
    end
    
    Lchrom = sign(uniquestim(:,1)) ~= sign(uniquestim(:,2));
    Llum = sign(uniquestim(:,1)) == sign(uniquestim(:,2)) & uniquestim(:,1)~= 0;
    
    figure; axes; hold on;
    h = [];
    h(1) = plot(uniquestim(Lchrom,3), dprime(Lchrom),'ro-','MarkerFaceColor','red');
    h(2) = plot(uniquestim(Llum,3), dprime(Llum),'ko-','MarkerFaceColor','black');
    h(3) = plot(uniquestim(Lchrom,3), conemodeldata(Lchrom),'ro:');
    h(4) = plot(uniquestim(Llum,3), conemodeldata(Llum),'ko:');
    set(gca,'Xscale','log');
    ylabel('d-prime');
    xlabel('TF (hz)');
    legend(h,{'Neuron L-M','Neuron L+M','Cones L-M','Cones L+M'});
    title(['Assuming neuron is ',choice]);
end
%%
% Section 8
% In case two cells were recorded simultaneously compute the cross
% correlogram.
rfx = stro.sum.exptParams.rf_x/10;
rfy = stro.sum.exptParams.rf_y/10;

if any(strcmp(stro.sum.rasterCells, 'sig001a')) & any(strcmp(stro.sum.rasterCells, 'sig001b'))
    spikes1 = stro.ras(:,strcmp(stro.sum.rasterCells, 'sig001a'));
    spikes2 = stro.ras(:,strcmp(stro.sum.rasterCells, 'sig001b'));
    bins = linspace(0,dur,dur*5000);
    data = zeros(1,2*(length(bins)-1)-1);
    for i = 1:size(stimon_t,1)
        tmp1 = histcounts(spikes1{i}-stimon_t(i) ,bins);
        tmp2 = histcounts(spikes2{i}-stimon_t(i) ,bins);
        if sum(tmp1 > 0) & sum(tmp2 > 0)
            data = data+xcorr(tmp1, tmp2,'coeff');
        end
    end
    figure; axes;
    binwidth = bins(2)-bins(1);
    timebins = [-1*floor(size(data,2)/2):1:floor(size(data,2)/2)]*binwidth;
    plot(timebins,data./size(stimon_t,1),'ko-');
    xlabel('lag (s)');
    title('negative numbers mean spike 2 leads spike 1');
    
    % repeating the analysis on a condition-by-condition basis
    data = zeros(size(uniquestim,1),2*(length(bins)-1)-1); % reusing this variable
    for i = 1:size(uniquestim,1)
        L = Lcc == uniquestim(i,1) & Mcc == uniquestim(i,2) & TF == uniquestim(i,3);
        for j = find(L)'
            tmp1 = histcounts(spikes1{i}-stimon_t(i),bins);
            tmp2 = histcounts(spikes2{i}-stimon_t(i),bins);
            data(i,:) = data(i,:)+xcorr(tmp1, tmp2,'coeff');
        end
    end
    Ltime = timebins > -0.0015 & timebins < 0;
    for stimtype = 1:2
        if (stimtype == 1)
            L = sign(uniquestim(:,1)) == sign(uniquestim(:,2));
            L = L & uniquestim(:,1) > 0 & uniquestim(:,2) > 0; % avoiding the blanks
        else
            L = sign(uniquestim(:,1)) ~= sign(uniquestim(:,2));
        end
        figure;
        idxs = find(L);
        for i=1:length(idxs)
            subplot(sum(L),1,i);
            plot(timebins(Ltime),data(idxs(i),Ltime),'ko-');
            ylabel(num2str(uniquestim(idxs(i),3)))
        end
        subplot(sum(L),1,1); % Adding a title to the top
        if (stimtype == 1)
            title('LUM');
        else
            title('RG');
        end
        equatesubplotaxeslims;
    end
    % Looking at trial-by-trial correlations in spike counts (within
    % condition)
    data = zeros(size(stro.trial,1),2);
    rs = [];
    %figure; axes; hold on;
    for i = 1:size(uniquestim,1)
        L = Lcc == uniquestim(i,1) & Mcc == uniquestim(i,2) & TF == uniquestim(i,3);
        data = [];
        for j = find(L)'
            data = [data; sum(spikes1{j} > stimon_t(j) & spikes1{j} < stimoff_t(j)) sum(spikes2{j} > stimon_t(j) & spikes2{j} < stimoff_t(j))];
        end
        %h=plot(data(:,1),data(:,2),'ko');
        tmp = corrcoef(data);
        rs(i) = tmp(1,2);
    end
    figure; axes; hold on;
    L = sign(uniquestim(:,1)) == sign(uniquestim(:,2));
    plot(uniquestim(L,3),rs(L),'ko-');
    plot(uniquestim(~L,3),rs(~L),'ro-');
    set(gca,'Xscale','log');
end


%% Section 9
% Comparing contrasts used to model fits
LMTFdata = load(fullfile(fileparts(which('IsoSampOnline')), 'private', 'data', 'LMTF.mat'));
[~,fname] = fileparts(stro.sum.fileName);
sid = fname(1);
LMTFstruct = getfield(LMTFdata,sid);
global_params = LMTFstruct.legacy.mode5params;
local_params_from_repository = LMTF_global_to_local_model(global_params, rfx, rfy, 5);
local_params_from_file = stro.sum.exptParams.localparams;

[xx,yy,zz] = meshgrid(linspace(-max(abs(Lcc)),max(abs(Lcc)),50),...
     linspace(-max(abs(Mcc)),max(abs(Mcc)),50),...
     linspace(1,max(TF),50));
surfcolors = {'green','blue'};
figure; axes; hold on;
for i = 1:2
    if i == 1
        fpar = local_params_from_repository;
    else 
        fpar = local_params_from_file;
    end
    xi_1 = fpar(1);
    zeta_1 = fpar(2);
    n1_1 = fpar(3);
    n2_1 = fpar(3)+fpar(4); % convention: n2 = n1+delta_n
    tau1_1 = 10^fpar(5);
    tau2_1 = 10^(fpar(5)+fpar(6)); % convention: tau2 = kappa*tau1
    xi_2 = fpar(7);
    zeta_2 = fpar(8);
    n1_2 = fpar(9);
    n2_2 = fpar(9)+fpar(10);
    tau1_2 = 10^fpar(11);
    tau2_2 = 10^(fpar(11)+fpar(12));
    theta = fpar(13);
    
    f1 = @(omega)xi_1*abs(((1i*2*pi*tau1_1.*omega+1).^-n1_1)-zeta_1*((1i*2*pi*tau2_1.*omega+1).^-n2_1));
    f2 = @(omega)xi_2*abs(((1i*2*pi*tau1_2.*omega+1).^-n1_2)-zeta_2*((1i*2*pi*tau2_2.*omega+1).^-n2_2));
    
    a = abs(f1(zz)); % luminance sensitivity
    b = abs(f2(zz)); % chromatic sensitivity
    mechs = [cos(theta) 1/sqrt(2); sin(theta) -1/sqrt(2)];
    V = zeros(size(xx));
    for j = 1:numel(xx)
        stimvector = [xx(j) yy(j)];
        V(j) = sqrt(sum((stimvector*mechs*diag([a(j) b(j)])).^2));
    end
    
    FV = isosurface(xx,yy,zz,V,1);
    h = patch(FV);
    set(h,'FaceColor',surfcolors{i},'FaceAlpha',.2','EdgeColor','none','EdgeAlpha',0);
end

% Now plotting a representation of the data
[uniquestim, dprime] = IsoSampGetDPrime(stro,2, find(spikeidx));

set(gcf,'Renderer','painters');
cmap = jet(64);
for i = 1:size(uniquestim,1)
    h = plot3(uniquestim(i,1),uniquestim(i,2),uniquestim(i,3),'ko','MarkerSize',15);
    set(h,'MarkerFaceColor',cmap(max(1,round(dprime(i)/max(dprime)*64)),:));
end
set(gca,'ZScale','log','View',[180 20])

%%
% Section 10 Trying the Victor spike train distance metric to classify
% spike trains a signal or noise.
% See StatsStuff.m section 14.

Q = .1; % Cost for moving a spike by 1 s (adding or removing a spike costs '1')
signaloffset = .05;

DEBUG = 0;
if DEBUG % Randomizing the Lcc, Mcc, and TF values. This is destructive!
    idxs = randperm(length(Lcc));
    Lcc = Lcc(idxs);
    Mcc = Mcc(idxs);
    TF = TF(idxs);
end
stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimon_t'));
stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimoff_t'));
dur = mean(stimoff_t-stimon_t);
spikeidx = ismember(stro.sum.rasterCells(1,:),{'sig001a'});

% Setting up the spikes so that time 0 is stim on and cutting off
% spikes after stimoff
ntrials = length(stro.ras(:,spikeidx));
spikes = {};
startidxs = nan*ones(ntrials,1);
endidxs = nan*ones(ntrials,1);

for i = 1:ntrials
    tmpspikes = stro.ras{i,spikeidx}-stimon_t(i);
    tmpspikes(tmpspikes<signaloffset) = [];
    tmpspikes(tmpspikes>dur) = [];
    if i == 1
        if isempty(tmpspikes)
            startidxs(i) = 0;
            endidxs(i) = -1;
        end
        startidxs(i) = 1;
        endidxs(i) = length(tmpspikes);
    else
        if isempty(tmpspikes)
            startidxs(i) = startidxs(i-1);
            endidxs(i) = startidxs(i)-1;
        end
        startidxs(i) = numel(vertcat(spikes{:}))+1;
        endidxs(i) = startidxs(i)+length(tmpspikes)-1;
    end
    spikes{i} = tmpspikes;
end
spiketimes = vertcat(spikes{:});

%figure; axes; hold on;
%for i = 1:length(startidxs)
%    plot(spiketimes(startidxs(i):endidxs(i)))
%end
% spkdl is blazingly fast so calculating *all* pairwise distances
d = spkdl(spiketimes,startidxs,endidxs,Q); % 10 ms
d_mat = reshape(d,length(startidxs),length(startidxs));
% imagesc(d_mat);

% ---------------------------------------------------------
% Doing the k-nearest neighbors calculation
% Getting the blanks in the correct format for spkdl
Lblank = Lcc == 0 & Mcc == 0;
nblanktrials = sum(Lblank);
data = [];
for i = 1:size(uniquestim,1)
    Lstim = Lcc == uniquestim(i,1) & Mcc == uniquestim(i,2) & TF == uniquestim(i,3);
    nsignaltrials = sum(Lstim);
    
    ntrials_tot = nblanktrials+nsignaltrials;
    sub_d_mat = [d_mat(Lblank,Lblank), d_mat(Lblank,Lstim);...
        d_mat(Lstim,Lblank), d_mat(Lstim,Lstim)];
    sub_d_mat = sub_d_mat+diag(nan(ntrials_tot,1));
    sorted_d_mat = sort(sub_d_mat); % sorts each column in ascending order (so first rows correspond to nearby spike trains)
    k = 2*floor((nsignaltrials-1)/2)+1; % an odd number of nearest neighbors to avoid ties
    critical_values = sorted_d_mat(k,:); % K nearest neighbor
    k_close_points = sub_d_mat<=repmat(critical_values,size(sub_d_mat,1),1);
    if ~all(sum(k_close_points) == k)
        error('bug in nearest neighbor selection');
    end
    
    mask = [-1*ones(nblanktrials,ntrials_tot); ones(nsignaltrials,ntrials_tot)]; % -1 means blank, 1 means stim
    whichcategoryassigned = sum(k_close_points.*mask)>0; % false (0) means classified as a blank, true (1) means classified as a stim
    % For blank-against-blank classification *every* trial is
    % misclassified because every trial has a duplicate of itself with
    % distance 0
    if any(sum(k_close_points.*mask) == 0)
        disp('got a tie');
        if i ~= 1 & ~DEBUG
            keyboard
        end
        tieidxs = sum(k_close_points.*mask) == 0;
        whichcategoryassigned(tieidxs) = unidrnd(2,1,sum(tieidxs))-1; % flipping a coin
    end
    
    nmisclassified_blank = sum(whichcategoryassigned(1:nblanktrials) ~= 0);
    nmisclassified_signal = sum(whichcategoryassigned(nblanktrials+1:ntrials_tot) == 0);
    % Each term in the nmisclassified sum, above, is biased upward because there are more "unlike" trials than "like" trials
    % similarly, ncorrectlyclassified is biased downward
    ncorrectlyclassified_blank = sum(whichcategoryassigned(1:nblanktrials) == 0);
    ncorrectlyclassified_signal = sum(whichcategoryassigned(nblanktrials+1:ntrials_tot) ~= 0);
    
    % Prob. that a blank is correctly classified is:
    % P(# nearby blanks > # nearby stims) =
    % P(# nearby blanks > k/2) =
    % 1-hygecdf(k/2,ntrials_tot-1, nblanktrials-1, k)
    % probability that a blank is *incorrectly* classified is
    % hygecdf(k/2,ntrials_tot-1, nblanktrials-1, k)
    % So E(# of correctly classified blanks) = nblanktrials*P(a blank is correctly classified)
    E_nmisclassified_blank = nblanktrials*hygecdf(k/2,ntrials_tot-1, nblanktrials-1, k);
    E_nmisclassified_signal = nsignaltrials*hygecdf(k/2,ntrials_tot-1, nsignaltrials-1, k);
    E_correctlyclassified_blank = nblanktrials-E_nmisclassified_blank;
    E_correctlyclassified_signal = nsignaltrials-E_nmisclassified_signal;
    
    % Gives same answer as: nblanktrials*(1-hygecdf(k/2,ntrials_tot-1, nblanktrials-1, k)) + nsignaltrials*(1-hygecdf(k/2,ntrials_tot-1, nsignaltrials-1, k))
    
    data = [data; nmisclassified_blank nmisclassified_signal nblanktrials nsignaltrials E_nmisclassified_blank E_nmisclassified_signal];
    % figure; subplot(2,2,1); hist(data(2:end,2))
    % data is in # correct classified beyond what is expected by chance.
end

% Correction for zeros and ones in hits and FAs
% https://stats.stackexchange.com/questions/134779/d-prime-with-100-hit-rate-probability-and-0-false-alarm-probability
if any(data(:,1) == 0)
    data(data(:,1) == 0,1) = 0.5;
end
if any(data(:,2) == 0)
    data(data(:,2) == 0,2) = 0.5;
end
if any(data(:,1) == data(:,3))
    data(data(:,1) == data(:,3),1) = data(data(:,1) == data(:,3),3)-0.5;
end
if any(data(:,2) == data(:,4))
    data(data(:,2) == data(:,4),2) = data(data(:,2) == data(:,4),4)-0.5;
end

H = (data(:,4)-data(:,2))./data(:,4);
FA = data(:,1)./data(:,3); 
dprime_raw = (1-norminv(FA,0,1))-(1-norminv(H,0,1));

E_H = (data(:,4)-data(:,6))./data(:,4);
E_FA = data(:,5)./data(:,3); 
dprime_correction = (1-norminv(E_FA,0,1))-(1-norminv(E_H,0,1));
dprime_corrected = dprime_raw-dprime_correction;

% Comparing to parametric dprime
[~, dprime_parametric] = IsoSampGetDPrime(stro,1,1); % No offset

figure; axes; hold on;
Lblank = all(uniquestim == 0,2);
Lrg = sign(uniquestim(:,1)) ~= sign(uniquestim(:,2)) & ~Lblank;
Llum = sign(uniquestim(:,1)) == sign(uniquestim(:,2)) & ~Lblank;
plot(uniquestim(Lrg,3), dprime_corrected(Lrg),'ro--');
plot(uniquestim(Lrg,3), dprime_parametric(Lrg),'r-.');

plot(uniquestim(Llum,3), dprime_corrected(Llum),'ko--');
plot(uniquestim(Llum,3), dprime_parametric(Llum),'k-.');
set(gca,'Xscale','log');
xlabel('TF'); ylabel('"non-parametric" d-prime');
pause
[~,dprime] = IsoSampGetKNNDPrime(stro,Q,1,signaloffset,1); % Testing IsoSampGetKNNDPrime function (should overlie -- lines)

%%
% Section 11
% Comparing sensitivity of ideal observer of cone signals and LGN signals
% as a function of counting window. For each time window (2D) we have cone
% d'/LGN d' as a function of TF and color dir. Based on section 7.

flashtime_idx = strcmp(stro.sum.trialFields(1,:),'flash_time');
flashtime = unique(stro.trial(:,flashtime_idx))/1000; % in sec

%for spikeidx = find(spikeidxs)
    timewindoffset = [-.2 .1]; % relative to stimon, when to look at cone responses
    % Do not go more than 200 ms before stimon nor more than 200 ms after
    % stimoff because spikes were not collected during that time.
    
    [uniquestim, dprime] = IsoSampGetDPrime(stro,3,spikeidx,timewindoffset,1);
    
    fundamentals = load('T_cones_smj10');
    MMPERDEG = 0.223; % mm/deg
    DEGPERMM = 1/MMPERDEG; % deg/mm
    
    % Setting up for cone model
    params.runType = 'isosamp';
    params.obsMethod = 'obsMethod_filteredWtFxn';
    params.impulseResponse = 'rieke';
    params.DTV1_fname = [];
    params.DTNT_fname = [];
    params.unitTest = false;
    params.eqMosaic = false;
    params.saveDir = '/Users/greghorwitz/Documents/MATLAB';
    params.notes = 'IsoSamp test';
    params.parallelOperations = false;
    params.eyeType = 'monkey';
    params.eyeNumber = 1; % comparing sensitivity of a monocular neuron to cones
    params.coneSampRate = 2400;
    params.flatPowerSpect = false;
    params.enableScones = false;
    params.sacamp_deg = 0;
    params.sacdur_s = 0;
    params.stro = stro;
    params.timewindoffset = timewindoffset; % relative to stimon, when to look at cone responses
    % "gab.length" in DTcones_gh.m is .666 and taken from flashtime in stro structure.
    % cal.gamma = params.stro.sum.exptParams.gamma_table; % Only used for "dtv1". Can ignore.
    spds = params.stro.sum.exptParams.mon_spd;
    spds = SplineSpd([380:4:780]',reshape(spds,length(spds)/3,3),[380:5:780]');
    cal.monSpect = spds(:);
    M = fundamentals.T_cones_smj10*spds;
    cal.Mmtx = M(:);
    cal.frameRate = params.stro.sum.exptParams.framerate;
    cal.bkgndrgb = params.stro.sum.exptParams.bkgndrgb'; % intensities
    cal.fname = 'test';
    cal.monSpectWavelengths = linspace(380,780,101);
    cal.pixperdeg = params.stro.sum.exptParams.pixperdeg;
    params.monCalFile = cal;
    % End of cone model set up
    
    rfx = stro.sum.exptParams.rf_x/10;
    rfy = stro.sum.exptParams.rf_y/10;
    rf_r_deg = sqrt(rfx^2+rfy^2);
    rf_r_mm = (5/38)*(sqrt(1520*rf_r_deg + 177089) - 421);
    RF_diam_midget_mm = 8.64*rf_r_mm^1.04/1000; % midgets
    RF_diam_midget_deg = RF_diam_midget_mm*DEGPERMM;
    RF_diam_parasol_deg = 0.0683+0.0165*rf_r_deg; % Regression coefficients from Perry et al. 1984
    
    % query to the user
    choice = questdlg('What type of cell to simulate?','CellMenu','Magno','Parvo','Cancel','Cancel');
    if strcmp(choice,'Parvo')
        params.gab.sd = RF_diam_midget_deg;
    elseif strcmp(choice,'Magno')
        params.gab.sd = RF_diam_parasol_deg;
    else
        return
    end
    
    [gab, cones, mon, idlob, params] = DTcones_gh(params,0); 
    
    conemodeldata = nan*ones(size(uniquestim,1),1);
    for i = 1:size(idlob.analyticMean,1) % looping over color direction
        for j = 1:size(idlob.analyticMean(i,:),2) % looping over contrast/TF
            if ~isempty(idlob.analyticMean{i,j})
                tmp_lm_mu = idlob.analyticMean{i,j}([1 2]);
                tmp_lm_var = idlob.analyticVar{i,j}([1 2]);
                tf = gab.driftRates{i}(j);
                L = uniquestim(:,3) == gab.driftRates{i}(j) & sign(gab.colorDirs(i,1)) == sign(uniquestim(:,1));
                if (sum(L) == 0)
                    error('cannot find condition');
                end
                conemodeldata(L) = [sum(abs(tmp_lm_mu))/sqrt(sum(tmp_lm_var))];
            end
        end
    end
    
    Lchrom = sign(uniquestim(:,1)) ~= sign(uniquestim(:,2));
    Llum = sign(uniquestim(:,1)) == sign(uniquestim(:,2)) & uniquestim(:,1)~= 0;
    
    figure; axes; hold on;
    h = [];
    h(1) = plot(uniquestim(Lchrom,3), dprime(Lchrom),'ro-','MarkerFaceColor','red');
    h(2) = plot(uniquestim(Llum,3), dprime(Llum),'ko-','MarkerFaceColor','black');
    h(3) = plot(uniquestim(Lchrom,3), conemodeldata(Lchrom),'ro:');
    h(4) = plot(uniquestim(Llum,3), conemodeldata(Llum),'ko:');
    set(gca,'Xscale','log');
    ylabel('d-prime');
    xlabel('TF (hz)');
    legend(h,{'Neuron L-M','Neuron L+M','Cones L-M','Cones L+M'});
    title(['Assuming neuron is ',choice]);
%end



%% OUTDATED CODE BELOW

%%
% Compute "F1" from activity during and after the stimulus. This will provide 
% a reasonable signal to noise ratio that should be independent of incomplete 
% numbers of cycles or periodicity in the neuron's baseline firing rate.

% offsets = [.05 dur]; % relative to stimulus on (in s)
% spiketimes = {}; % First, getting all the spike times
% for i = 1:size(uniquestim,1)
%     current_TF = uniquestim(i,3);
%     L = Lcc == uniquestim(i,1) & Mcc == uniquestim(i,2) & TF == current_TF;
%     tmpspikes = [];
%     for j = find(L)'
%         tmpspikes = [tmpspikes; spikes{j}-stimon_t(j)];
%     end
%     spiketimes{i} = tmpspikes(tmpspikes > offsets(1) & tmpspikes < offsets(2));
% end
% 
% % For each frequency, projecting the spike times onto (brief) sin and
% % cos basis vectors (pooling spikes across trials within a condition).
% % Then taking the amplitude of the vector (cos^2+sin^2). This is done 
% % for both 0 contrast trials and trials at the particular frequency being
% % examined.
% % Unclear whether the projections should be normalized by the number of
% % spikes or the number of trials. 
% data = [];
% noiseidx = find(all(uniquestim == 0,2));
% for i = 1:size(uniquestim,1)
%     current_TF = uniquestim(i,3);
%     ntrials = sum(Lcc == uniquestim(i,1) & Mcc == uniquestim(i,2) & TF == uniquestim(i,3));
%     cos_n_sin = [cos(2*pi*current_TF.*spiketimes{i}) sin(2*pi*current_TF.*spiketimes{i})];
%     F1_sig = sum(sum(cos_n_sin).^2)./length(spiketimes{i})^2;
%     %F1_sig = sum(sum(cos_n_sin).^2)./ntrials;
%     cos_n_sin = [cos(2*pi*current_TF.*spiketimes{noiseidx}) sin(2*pi*current_TF.*spiketimes{noiseidx})];
%     F1_noise = sum(sum(cos_n_sin).^2)./length(spiketimes{noiseidx})^2;
%     %F1_noise = sum(sum(cos_n_sin).^2)./ntrials;
% 
%     data = [data; F1_sig F1_noise];
% end
% % ---------------
% if (any(data(:) == 0))
%     error('zero spikes in some condition');
% end
% 
% % One sample t-test on SNRs across TFs
% % Not sure how legit this is because the noise is not independent across
% % freqencies (same spike trains projected onto similar basis vectors)
% figure; cols = [1 0 0; 0 0 0];
% hist(log(data(:,1)./data(:,2)));
% [h,p] = ttest(log(data(:,1)./data(:,2)));
% xlabel('SNR'); ylabel('Stimulus conditions'); title(['p = ',num2str(p)]);
% 
% figure;
% subplot(2,2,1); title('LUM'); xlabel('TF (Hz)'); 
% subplot(2,2,2); xlabel('TF (Hz)'); ylabel('SNR');
% subplot(2,2,3); title('RG'); xlabel('TF (Hz)');
% subplot(2,2,4); xlabel('TF (Hz)'); ylabel('SNR');
% 
% for stimtype = 1:2 % lum, color
%     if (stimtype == 1)
%         L = sign(uniquestim(:,1)) == sign(uniquestim(:,2));
%         L = L & uniquestim(:,1) > 0 & uniquestim(:,2) > 0; % avoiding the blanks
%     else
%         L = sign(uniquestim(:,1)) ~= sign(uniquestim(:,2));
%     end
%     for j = 1:size(uniquestim,1)
%         if (L(j))
%             subplot(2,2,1+2*(stimtype-1)); hold on;       
%             plot(uniquestim(j,3),data(j,1),'o','MarkerSize',6,'MarkerFaceColor',[0 .75 0],'MarkerEdgeColor','none');
%             plot(uniquestim(j,3),data(j,2),'o','MarkerSize',6,'MarkerFaceColor',[.5 .5 .5],'MarkerEdgeColor','none');
%             set(gca,'Xlim',[min(uniquestim(L,3)) max(uniquestim(L,3))],'Xscale','log','Yscale','log');
%             subplot(2,2,2+2*(stimtype-1)); hold on;      
%             plot(uniquestim(j,3),data(j,1)./data(j,2),'o','MarkerSize',6,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor','none');
%             set(gca,'Xlim',[min(uniquestim(L,3)) max(uniquestim(L,3))],'Xscale','log');
%             set(gca,'Ylim',[min(data(:,1)./data(:,2)) max(data(:,1)./data(:,2))],'Yscale','log');
%         end
%     end
%     if any(L)
%         plot([min(uniquestim(L,3)) max(uniquestim(L,3))],[1 1],'k--');
%     end
% end
% 
% % Bubble plot of SNR
% SNR = data(:,1)./data(:,2);
% normSNR = 127/max(abs(log(SNR)))*log(SNR)+127;
% cmap = jet(255);
% figure; axes; hold on;
% for j = 1:size(uniquestim,1)
%     L = Lcc == uniquestim(j,1) & Mcc == uniquestim(j,2) & TF == uniquestim(j,3);
%     h = plot3(uniquestim(j,1),uniquestim(j,2),uniquestim(j,3),'ko','MarkerSize',normSNR(j)/10+1);
%     set(h,'MarkerEdgeColor','none','MarkerFaceColor',cmap(round(normSNR(j))+1,:));
% end
% standard = plot3(0,0,1,'ko','MarkerSize',128/10+1);
% set(standard,'MarkerEdgeColor','none','MarkerFaceColor',cmap(128,:));
% 
% set(gca,'View',[50 20], 'Zscale','log','Zlim',[1 max(uniquestim(:,3))]);
% xlabel('L-cone contrast'); ylabel('M-cone contrast'); zlabel('TF'); 
% title('Symbol size is SNR');
% axis vis3d;

%%
% Section 12 
% Comparing the rate of SNR increase for the photon absorption ideal
% observer and the cone current ideal observer.

% The answer is that the SNR, in both cases, increases *linearly* with the
% contrast. This makes sense for the cone current observer because the noise
% is fixed and additive. It's true for the photon absorption observer
% because the noise is basically identical for the noise and all levels of 
% contrast signal (remember, the noise is not Poisson, its the sum of Poisson
% RVs, each of which has large, similar variance).
% The tiny change in the ratio of SNRs is due to the fact
% that, for the photon absorption ideal observer, the variance does change
% a tiny bit with the amplitude of the signal.

filename = 'A061517005.nex'; % example magnocell
stro = nex2stro(findfile(filename,[nexfilepath,filesep,'Greg',filesep,'Apollo']));
spikeNum = 'sig001a';

contrasts = linspace(0,1,100);
stro.trial([1:length(contrasts)],20) = contrasts; % Lcc
stro.trial([1:length(contrasts)],21) = contrasts; % Mccs
stro.trial([1:length(contrasts)],17) = 8; % TF

stro.trial(length(contrasts)+1:end,:) = [];
stro.ras(length(contrasts)+1:end,:) = [];

fundamentals = load('T_cones_smj10');
MMPERDEG = 0.223; % mm/deg
DEGPERMM = 1/MMPERDEG; % deg/mm

% Setting up for cone model
params.runType = 'isosamp';
params.obsMethod = 'obsMethod_filteredWtFxn';
params.impulseResponse = 'rieke';
params.DTV1_fname = [];
params.DTNT_fname = [];
params.unitTest = false;
params.eqMosaic = false;
params.saveDir = '/Users/greghorwitz/Documents/MATLAB';
params.notes = 'IsoSamp test';
params.parallelOperations = false;
params.eyeType = 'monkey';
params.eyeNumber = 1;
params.coneSampRate = 2400;
params.flatPowerSpect = false;
params.enableScones = false;
params.sacamp_deg = 0;
params.sacdur_s = 0;
params.stro = stro;
%cal.gamma = params.stro.sum.exptParams.gamma_table; % Only used for "dtv1". Can ignore.
spds = params.stro.sum.exptParams.mon_spd;
spds = SplineSpd([380:4:780]',reshape(spds,length(spds)/3,3),[380:5:780]');
cal.monSpect = spds(:);
M = fundamentals.T_cones_smj10*spds;
cal.Mmtx = M(:);
cal.frameRate = params.stro.sum.exptParams.framerate;
cal.bkgndrgb = params.stro.sum.exptParams.bkgndrgb'; % intensities
cal.fname = 'test';
cal.monSpectWavelengths = linspace(380,780,101);
cal.pixperdeg = params.stro.sum.exptParams.pixperdeg;
params.monCalFile = cal;
% End of cone model set up

uniquestim = IsoSampGetDPrime(stro);

[gab, cones, mon, idlob, params] = DTcones_gh(params,0);
conemodeldata = nan*ones(length(idlob.analyticMean),1);
for i = 1:length(idlob.analyticMean) % looping over contrast/TF (only one color direction)
    if ~isempty(idlob.analyticMean{i})
        tmp_lm_mu = idlob.analyticMean{i}([1 2]);
        tmp_lm_var = idlob.analyticVar{i}([1 2]);
        conemodeldata(i) = [sum(abs(tmp_lm_mu))/sqrt(sum(tmp_lm_var))];
    end
end


t = 0:1/mon.frameRate:gab.length;
nframes = length(t);
flashTimeProfile = ones(1,nframes);
ramp = linspace(0,1,nframes/4);
flashTimeProfile(1:length(ramp)) = ramp;
flashTimeProfile(end:-1:end-length(ramp)+1) = ramp;
photon_dprimes = IsoSampGetPhotonDPrime (flashTimeProfile, mon.frameRate, mon.bkgndlms_Rstar, params.gab.sd*cal.pixperdeg, cat(3,cones.num_L,cones.num_M), uniquestim);

figure; 
rms = sqrt(uniquestim(:,1).^2+uniquestim(:,2).^2);
subplot(2,1,1); hold on;
plot(rms(rms>0), conemodeldata,'kv');
plot(rms(rms>0), photon_dprimes(rms>0),'k*');
subplot(2,1,2);
plot(rms(rms>0), conemodeldata./photon_dprimes(rms>0),'rv');