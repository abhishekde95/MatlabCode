%%
% A couple of analyses for Patrick's experiment

stro = nex2stro;
Lcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'Lcc'));
Mcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'Mcc'));
fpacq_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fpacq_t'));
stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimon_t'));
stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fpoff_t'));
direction = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'Theta'));
ntrials = size(stro.trial,1);
LM = Lcc+Mcc;
LvM = (Lcc-Mcc);

bins = linspace(0,2,50);
PSTH = zeros(1,length(bins));
% First PSTHs to find good values for "offset"
for i = 1:ntrials
    spiketimes = stro.ras{i,1}-stimon_t(i);
    PSTH = PSTH +histc(spiketimes',bins);
end

% Getting baseline firing rate
dt = zeros(ntrials,1);
baselinerates = zeros(ntrials,1);
for i = 1:ntrials
    nsp =  sum(stro.ras{i,1} > fpacq_t(i)+.05 & stro.ras{i,1} < stimon_t(i));
    baselinerates(i) = nsp./(stimon_t(i) - fpacq_t(i) -.05);
    dt(i) = stimon_t(i) - fpacq_t(i);
end

figure; axes; hold on;
bar(bins,PSTH);
gauss = @(x,params)(params(2)+params(1)*exp((-(x-params(3)).^2)/params(4).^2));
bguess = mean(PSTH(bins<.25));
aguess = max(PSTH)-bguess;
muguess = bins(find(PSTH == max(PSTH),1,'first'));
sigmaguess = .1; % terrible. Need something better.
%plot(bins,gauss(bins,[aguess,bguess,muguess,sigmaguess]),'m.')
fittedparams = fminsearch(@(params)sum((gauss(bins,params)-PSTH).^2),[aguess,bguess,muguess,sigmaguess]);
plot(bins,gauss(bins,fittedparams),'b-.','LineWidth',3)
offset = [fittedparams(3)-fittedparams(4) fittedparams(3)+fittedparams(4)];
plot(offset,[fittedparams(1) fittedparams(1)],'k-','linewidth',3);

% Now counting up spikes in a window
ntrials = size(stro.trial,1);
nspikes = zeros(ntrials,1);
spikerates = zeros(ntrials,1);
for i = 1:ntrials
    nspikes(i) = sum(stro.ras{i,1} > stimon_t(i)+offset(1) & stro.ras{i,1}<stimon_t(i)+offset(2));
    dt = offset(2)-offset(1);
    spikerates(i) = nspikes(i)./dt;
end

%%
% Generating an initial guess 

vlb = [0    0 0.001 0.001 1 1 0];
vub = [200 200    1     1 6 6 100];
options = optimset('MaxFunEvals',50000,'MaxIter',50000,'TolFun',10^-6,'TolX',10^-6);

sigmaguess = 0.1;
Ldir = logical(direction == 1);
if sum(Ldir) == 0; Ldir = true(length(Ldir),1); end
maxpred = []; f = [];
titles = {'L+M','L-M'};
for i = 1:2
    if (i == 1)
        L = Ldir & LvM == 0;
        projs = LM(L);
    else
        L = Ldir & LM == 0;
        projs = LvM(L);
    end
    
    topfr = max(nspikes(L));
    if (topfr == 0)
        f(i,:) = [0 0 1 1 2 2 0]
    else
        params0 = [topfr, topfr, sigmaguess, sigmaguess, 2, 2, 0];  % need to constrain B across color directions
        f(i,:) = fmincon('MyFitNakaRushtonFun',params0,[],[],[],[],vlb,vub,[],options,projs,nspikes(L),'asymmetric');
    end
    pred = MyComputeNakaRushton(f(i,:),projs,'asymmetric');
    subplot(2,1,i); hold on;
    plot(pred,nspikes(L),'.');
    plot([0 max(pred)],[0 max(pred)],'k-');
    xlabel('pred'); ylabel('actual response');
end
a = mean(f(1,[1 2]));
b = mean(f(2,[1 2]));
%params0(1) = max(max(f([1 2],[1 2])));
params0 = max(nspikes);
params0(2) = a./sum([a b]); %this doesn't work so well
params0(2) = .5;
params0(3) = .25;
params0(4) = 2; 
params0(5) = (f(1,2)./f(1,4))./(f(1,1)./f(1,3)); if (isnan(params0(5))) params0(5) = 1; end;
params0(6) = (f(2,2)./f(2,4))./(f(2,1)./f(2,3)); if (isnan(params0(6))) params0(6) = 1; end;%
%params0(5) = 1;
%params0(6) = 1;
params0(7) = mean(f(:,end))

%%
vlb = [10   0 0.0001 1  0  0  0];
vub = [1000 1     1  6 10 10 20];
options = optimset('MaxFunEvals',50000,'MaxIter',50000,'TolFun',10^-6,'TolX',10^-6,'Display','on');


Ldir = logical(direction == 1);
if sum(Ldir) == 0; Ldir = true(length(Ldir),1); end

f0 = fmincon('MyFitNakaRushtonFun',params0,[],[],[],[],vlb,vub,[],options,[LM(Ldir) LvM(Ldir)],nspikes(Ldir),'surface1');
[x,y] = meshgrid(linspace(min(LM),max(LM),50), linspace(min(LvM),max(LvM),50));
surface = MyComputeNakaRushton(f0,[x(:) y(:)],'surface1');
figure; axes; hold on;
surf(x,y,reshape(surface,size(x)))
plot3(LM(Ldir), LvM(Ldir), nspikes(Ldir),'ko');
title(['Fit to all the data (n = ',num2str(sum(Ldir)),')']);

% Let's try fitting to the L+M and L-M directions only
% Minor cheat - using result from one analysis to feed the other
L = LM == 0 | LvM == 0;
f1 = fmincon('MyFitNakaRushtonFun',f0,[],[],[],[],vlb,vub,[],options,[LM(L&Ldir) LvM(L&Ldir)],nspikes(L&Ldir),'surface1');
[x,y] = meshgrid(linspace(min(LM),max(LM),50), linspace(min(LvM),max(LvM),50));
surface = MyComputeNakaRushton(f1,[x(:) y(:)],'surface1');
figure; axes; hold on;
surf(x,y,reshape(surface,size(x)))
plot3(LM(Ldir), LvM(Ldir), nspikes(Ldir),'ko');
title(['Fit to L+M and L-M only (n = ',num2str(sum(L&Ldir)),')']);

% Looking at residuals, first for one model, then the other.
% Predicting data points (~L) because we're looking at the off-cardinal data points
for i = 1:2
    if i == 1
        pred = MyComputeNakaRushton(f0,[LM(~L&Ldir) LvM(~L&Ldir)],'surface1');
    else
        pred = MyComputeNakaRushton(f1,[LM(~L&Ldir) LvM(~L&Ldir)],'surface1');
    end
    
    % making predictions
    figure; subplot(2,1,1); hold on;
    plot(pred,nspikes(~L&Ldir),'.');
    plot([0 max(pred)],[0 max(pred)],'k-')
    xlabel('pred'); ylabel('actual response');
    subplot(2,1,2); hold on;
    resid = nspikes(~L&Ldir)-pred;
    plot(pred,resid,'.');
    plot([min(pred) max(pred)],[0 0],'k-')
    xlabel('pred'); ylabel('residual');
    
    % residuals in the plane
    figure; hax(1) = subplot(1,2,1); hold on;
    plot3(LM(~L&Ldir), LvM(~L&Ldir), nspikes(~L&Ldir)-pred,'k.')
    out2d = tpaps([LM(~L&Ldir), LvM(~L&Ldir)]', nspikes(~L&Ldir)'-pred');
    h = fnplt(out2d);
    surf(h{1},h{2},h{3});
    f = sum(pred)-(log(pred')*nspikes(~L&Ldir));
    axis tight; axis square; title(['residuals: ',num2str(round(f))]);
    
    hax(2) = subplot(1,2,2); hold on;
    plot3(LM(Ldir), LvM(Ldir), nspikes(Ldir),'r.')
    out2d = tpaps([LM(Ldir), LvM(Ldir)]', nspikes(Ldir)');
    h = fnplt(out2d);
    surf(h{1},h{2},h{3});
    axis tight; axis square; title('TPS to raw data');
    
    xlims = [get(hax(1),'Xlim'); get(hax(2),'Xlim')];
    set(hax,'XLim',[max(xlims(:,1)), min(xlims(:,2))]);
    ylims = [get(hax(1),'Ylim'); get(hax(2),'Ylim')];
    set(hax,'YLim',[max(ylims(:,1)), min(ylims(:,2))]);
    zlims = [get(hax(1),'Zlim'); get(hax(2),'Zlim')];
    set(hax,'ZLim',[min(zlims(:,1)), max(zlims(:,2))]);
end