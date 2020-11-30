%%
% Contents
% Section 1) Comparing actual neurometric thresholds from those predicted by a
% cardinal, linear model.
% Section 2) Fitting a half-squaring model to neurophysiological data.
% Section 3) Showing that we can get cone weight from Gaussina RGB noise
% (for the reviewers of the Journal of Vision paper)
% Section 4) Using the Zohary et al. formalism for computing the signal-to-noise
% ratio of a pool of correlated neurons
% Section 5) Working on the "potpourri" figure for the JOV paper 
% Section 6) Hacking around with the sum of Poisson RVs.
% Section 7) Sensitivity of ideal observer with access to Poisson
% distributed photon isomerizations (Figure 3).
% Section 8) Ratio of sensitivities between Charlie's cone model and the
% ideal observer photon absorptions (Figure 3).
% Section 9) Rendering nicely shaded isodetection surface figures (and an
% icon for the JOV paper)
% Section 10) Making an artificial cone mosaic
% Section 11) Counting cones inside the Gabor
% Section 12) Looking at cone-opponent and non-opponent signals inside the
% Gabor
% Section 13) More isodetection surface rendering stuff
% Section 14) Examining how errors in cone fundamentals would manifest in
% cone model sensitivity
% Section 15) Examining how errors in L:M cone ratios would manifest in
% cone model sensitivity
% Section 16) Getting data for Xiaomao, David Brainard's student. Seeing
% whether asymmetry in orange-cyan/lime-magenta detection might be due to
% errors in pre-retinal filtering estimates
%% LINEAR MECHANISMS MODEL (Section 1)
%cycle through the valid neurons and compare the actual and predicted TRs
%of the intermediate color direction based on a linear model
l_validConds = ~(out.errors(:,cardTieInd) | out.errors(:,intTieInd) | out.errors(:, lt10spikesInd) | out.errors(:, lt24DTtrialsInd) | out.errors(:,lt5GTtrialsInd));
validCells = find(l_validConds);
for a = 1:length(validCells)
    cellInd = validCells(a);
    
    cardNT = out.dat(cellInd).c.alpha(out.maskCardTR(cellInd,:)');
    cardColor = out.dat(cellInd).expt.standColors(out.maskCardTR(cellInd,:)', :);
    cardUnit(a,:) = cardColor ./ norm(cardColor);
    cardVect = cardUnit(a,:)*cardNT;
    
    measIntNT(a) = out.dat(cellInd).c.alpha(out.maskIntTR(cellInd,:)');
    intColor = out.dat(cellInd).expt.standColors(out.maskIntTR(cellInd,:)', :);
    intUnit(a,:) = intColor ./ norm(intColor);
    intVect = intUnit(a,:)*measIntNT(a);
     
    % Under the linear model, a light is detected when the dot product of any light 
    % onto the mechanism, cardUnit(a,:), equals cardNT
    
    % It works for the cardinal light at threshold...
    [cardVect*cardUnit(a,:)'  cardNT]
    
    % ...but it doesn't work for the intermediate light at threshold.
    [abs(intVect*cardUnit(a,:)')  cardNT]
    
    % What do we have to multiply intVect by so that its
    % dot product onto the mechanism, cardUnit(a,:), is cardNT?
    % (x*light2*mech) = light1*mech
    % x = (light1*mech)/(light2*mech)
    intScaleFactor = (cardVect*cardUnit(a,:)') / abs(intVect*cardUnit(a,:)');
    
    % intScaleFactor is how much do we have to multiply intVect by so that
    % the dotproduct between intVect and cardUnit is same as dot prod
    % between cardVect and cardUnit.
    % The two expressions below should evaluate to the same number.
    % They are (1) the projection of the cardinal light onto the
    % cardinal mechanism and (2) the projection of the scaled intermediate
    % light onto the cardinal mechanism.
    [cardVect*cardUnit(a,:)' abs(intScaleFactor*intVect*cardUnit(a,:)')]
    
    % We can use this to derive a prediction about the threshold in the
    % intermediate color direction.
  	predIntNT(a) = intScaleFactor*norm(intVect)
    
    % To make sure this analysis is invariant to how we represent the data,
    % we are going to linearly transformation of the space by multiplying 
    % everything by some random matrix, A.
    A = normrnd(0,1,3,3);
    mech = inv(A')*(cardUnit(a,:))';
    light1 = A*cardVect';
    light2 = A*intVect';

    % Now just repeating the above code with the new mechanism and new
    % lights.
    newScaleFactor(a) = (mech'*light1) / abs(mech' * light2);
    newpredIntNT(a) = intScaleFactor*norm(light2);
    newmeasIntNT(a) = norm(light2);
end


% Hopefully, the ratio of predicted to actual neurometric thresholds
% is the same with and without the random transformation of the space.  
% Otherwise, the analysis is lame.
newpredIntNT./newmeasIntNT
predIntNT./measIntNT
% On the other hand, we expect the absolute values of the thresholds to
% depend on the linear transformation (like converting from
% inches to centimeters - the absolute values of the numbers change but the
% ratios don't).

%plot the data
figure
subplot(1,2,1), hold on
plot(measIntNT,predIntNT , 'b.')
xlabel('measured Int NT')
ylabel('predicted Int NT')
maxVal = max([measIntNT(:) ; predIntNT(:)]) .* 1.1;
minVal = min([measIntNT(measIntNT>0)' ; predIntNT(predIntNT>0)']) .* 0.90;
plot([minVal maxVal], [minVal, maxVal], 'k-')
nNans = sum(isnan(measIntNT(:)+predIntNT(:))); %either the pred or meas TR is a nan
title(sprintf('n = %d, nNan = %d', sum(l_validConds), nNans));
axis tight
hold off

subplot(1,2,2), hold on
plot(predIntNT./measIntNT, newpredIntNT./newmeasIntNT, 'b.')
xlabel('Pred/Meas NT in cone space');
ylabel('Pred/Meas NT in ???? space');
x = get(gca,'XLim');
plot(x.*[.9 1.1], x.*[.9 1.1],'k-');
axis tight


%%
% Section 2) 
% Fitting a half-squaring (with a threshold) to simulated data
b = [3 1 2];
contrasts = linspace(0,4,6);
lambda = b(1)+b(3)*(max(contrasts-b(2),0)).^2;

ntrials = 7;
nspikes = nan*ones(ntrials,length(contrasts));
for i = 1:ntrials
    nspikes(i,:) = poissrnd(lambda,1,length(lambda));
end

% OK, I've got my fake data. Now I just need to fit it.
x = repmat(contrasts,ntrials,1); % making a design matrix
%beta = glmfit(x(:),nspikes(:),'poisson','link','identity');
beta = regress(nspikes(:),[ones(length(x(:)),1) x(:).^2]);
%

options = optimset;
options = optimset(options,'Diagnostics','off','Display','off');
options = optimset(options,'LargeScale','off');
vlb = [0 0 0];
vub = [50 50 50];

params0 = [beta(1) .1 beta(2)];
c = repmat(contrasts,ntrials,1);
CAHfit(params0,c(:),nspikes(:))

f1 = fmincon('CAHfit',params0,[],[],[],[],vlb,vub,[],options,c(:),nspikes(:));

figure; axes; hold on;
plot(contrasts,nspikes,'k.');
x = linspace(contrasts(1),contrasts(end),30);
% plot(x,beta(1)+beta(2)*x.^2,'b-');  % Initial fit
plot(x,b(1)+b(3)*(max(x-b(2),0).^2),'k:'); % truth
plot(x,f1(1)+f1(3)*(max(x-f1(2),0).^2),'m-'); % fmincon fit

% Here's the -llik  of the fit
CAHfit(f1,c(:),nspikes(:))

% ---------------------------------------
% Now trying a *pair* of color directions
% ---------------------------------------
b_int = [b(1) b(2) b(3)]
contrast_scaling = 1;  % bounded between .1 and 10
contrasts_int = [0:1:4]/contrast_scaling;
lambda_int = b_int(1)+b_int(3)*(max(contrasts_int.*contrast_scaling-b_int(2),0).^2);

ntrials_int = 7;
nspikes_int = nan*ones(ntrials_int,length(contrasts_int));
for i = 1:ntrials_int
    nspikes_int(i,:) = poissrnd(lambda_int,1,length(lambda_int));
end

% --------------------
% Fitting them individually
x = repmat(contrasts,ntrials,1); % making a design matrix
beta1 = regress(nspikes(:),[ones(length(x(:)),1) x(:).^2]);
params0 = [beta1(1) 0 beta1(2)];
c = repmat(contrasts,ntrials,1);
f1 = fmincon('CAHfit',params0,[],[],[],[],vlb,vub,[],options,c(:),nspikes(:));
dev1 = CAHfit(f1,c(:),nspikes(:));

x = repmat(contrasts_int,ntrials_int,1); % making a design matrix
beta2 = regress(nspikes_int(:),[ones(length(x(:)),1) x(:).^2]);
params0 = [beta2(1) 0 beta2(2)];
c = repmat(contrasts_int,ntrials_int,1);
f2 = fmincon('CAHfit',params0,[],[],[],[],vlb,vub,[],options,c(:),nspikes_int(:));
dev1(2)  = CAHfit(f2,c(:),nspikes_int(:));
dev1 = sum(dev1);

% plotting
figure; axes; hold on;
plot(contrasts,nspikes','k.');
plot(contrasts_int,nspikes_int','m.')
x = linspace(0,max([contrasts, contrasts_int]),100);
plot(x,f1(1)+f1(3)*(max(x-f1(2),0).^2),'k-');
plot(x,f2(1)+f2(3)*(max(x-f2(2),0).^2),'m-');
set(gca,'Ylim',[0 max([nspikes_int(:);nspikes(:)])]);

% --------------------
% Trying to fit both yoked
vlb = [0 0 0 .1];
vub = [50 50 50 10];

params1 = [mean([f1(1); f2(1)]) f1(2) f1(3) min([sqrt(f2(3)./f1(3)), f1(2)./f2(2)])];
params1 = [mean([f1(1); f2(1)]) f1(2) f1(3) sqrt(f2(3)./f1(3))];
params1 = [mean([f1(1); f2(1)]) f1(2) f1(3) contrast_scaling]; % cheating
c1 = repmat(contrasts,ntrials,1);
c2 = repmat(contrasts_int,ntrials_int,1);
X = [c1(:);c2(:)];
X(:,2) = [zeros(1,numel(c1)) ones(1,numel(c2))];
y = [nspikes(:); nspikes_int(:)];

f = fmincon('CAHfit',params1,[],[],[],[],vlb,vub,[],options,X,y);
dev2 = CAHfit(f,X,y);

% plotting
figure; axes; hold on;
plot(contrasts,nspikes','k.');
plot(contrasts_int,nspikes_int','m.')
x = linspace(0,max([contrasts, contrasts_int]),30);
plot(x,f(1)+f(3)*(max(x-f(2),0).^2),'k-'); % fmincon fit
plot(x,f(1)+f(3)*(max(x*f(4)-f(2),0).^2),'m-'); % fmincon fit
set(gca,'Ylim',[0 max([nspikes(:);nspikes_int(:)])]);

devstat = 2*(dev2-dev1);  % "full model" goes second. dev1 and dev2 are -1*llik
p = 1- chi2cdf(devstat,2);
title(['p = ',num2str(p)]);

%[b(1) b(2) b(3) b_int(3)]
%params1


%%
% Section 3)
% L, M, and S cone weights for an example neuron stimulated with Gaussian
% RGB noise. This is for a figure to give to the DTspot manuscript reviewers
whichframes = [5 6 7];

WN=nex2stro;
framerate = WN.sum.exptParams.framerate;
nstixperside = WN.sum.exptParams.nstixperside;
ntrials = length(WN.sum.absTrialNum);
stimonidx = find(strcmp(WN.sum.trialFields(1,:),'stim_on'));
stimoffidx = find(strcmp(WN.sum.trialFields(1,:),'all_off'));
nframesidx = find(strcmp(WN.sum.trialFields(1,:),'num_frames'));
noisetypeidx = find(strcmp(WN.sum.trialFields(1,:),'noise_type'));
sigmaidxs = strmatch('sigma',WN.sum.trialFields(1,:));

hepidx = find(strcmp(WN.sum.rasterCells(1,:),'AD11'));
vepidx = find(strcmp(WN.sum.rasterCells(1,:),'AD12'));
anlgStartTimeidx = find(strcmp(WN.sum.rasterCells(1,:),'anlgStartTime'));
eyestart_t = [WN.ras{:,anlgStartTimeidx}]';
eyesampperiod = 1/WN.sum.analog.storeRates{1};
gammaTable = WN.sum.exptParams.gamma_table;
gammaTable = reshape(gammaTable, length(gammaTable)/3, 3);
gammaTable1 = interp1(linspace(0,255,256),gammaTable,linspace(0,255,65536), 'spline');
invgamma = InvertGamma(gammaTable, 0);

% Reconstructing the M matrix
fundamentals = WN.sum.exptParams.fundamentals;
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
mon_spd = WN.sum.exptParams.mon_spd;
mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
M = fundamentals'*mon_spd;

% Getting the background rgb/lms
ridx = find(strcmp(WN.sum.trialFields(1,:),'bkgnd_r'));
gidx = find(strcmp(WN.sum.trialFields(1,:),'bkgnd_g'));
bidx = find(strcmp(WN.sum.trialFields(1,:),'bkgnd_b'));
bkgndRGB = [mode(WN.trial(:,ridx)), mode(WN.trial(:,gidx)), mode(WN.trial(:,bidx))];
bkgndrgb = [gammaTable(bkgndRGB(1)+1,1); gammaTable(bkgndRGB(2)+1,2); gammaTable(bkgndRGB(3)+1,3)];
bkgndlms = M*bkgndrgb;
spikename = getSpikenum(WN);
spikeidx = find(strcmp(WN.sum.rasterCells(1,:),spikename));
maxT = 9;
Lgunnoise = WN.trial(:,noisetypeidx) == 1;
Lconenoise = WN.trial(:,noisetypeidx) == 2;
WN.ras(Lconenoise,:) = [];
WN.trial(Lconenoise,:) = [];
out = getWhtnsStats(WN,maxT,'STCOVmex', {nstixperside^2, 3, maxT}, spikename);


nspikes = out{3};
STAs = [];
for i = 1:length(whichframes)
    lms = inv(M')*reshape(out{1}(:,whichframes(i)),100,3)';
    lms = lms';
  %  lms = reshape(out{1}(:,whichframes(i)),100,3); % debugging
    STAs(:,i) = lms(:);
end
STAs = mean(STAs,2);
maxes = max(STAs(:));
mins = min(STAs(:));

potentialnormfactors = [(1-[.5; .5; .5]-eps)./(maxes-[.5; .5; .5]); (-[.5; .5; .5]+eps)./(mins-[.5; .5; .5])];
% 'eps' in above line is a kludge that is required for avoiding
% out of bounds errors.
potentialnormfactors(potentialnormfactors < 0) = []; % if min > mu or max < mu
normfactor = min(potentialnormfactors);

%muvect = reshape(repmat(bkgndrgb',nstixperside^2,1),nstixperside^2*3,1);
muvect = reshape(repmat([.5 .5 .5],nstixperside^2,1),nstixperside^2*3,1);

% Plotting
figure;
for i = 1:size(STAs,2)
    STA = normfactor*(STAs(:,i)-muvect)+muvect;
    STA = reshape(STA,[nstixperside nstixperside 3]);
    for j = 1:3
        subplot(3,size(STAs,2),j);
        image(255*STA(:,:,j)); colormap(gray(255));
        set(gca,'XTick',[],'YTick',[]); axis square;
    end
end


%%
% Section 4)
% The signal to noise ratio for a population of V1 neurons
% using the formalism of Zohary et al 1994.

% If mu = 1.25 and sigma = 1, this implies d'= 1.25 =>  ROC area 0.81
% Average neuron is 2/3 as sensitive at the monkey, so monkey's d' is 1.875.

% How much do we need to crank M (number of neurons) by to get 1.875
% Cones are 3x more sensitive than the monkey and so 4.5x more sensitive
% than the V1 neurons

mu = 1.25;
sigma = 1;
M = 10;
r =0.2;


pooledSNR = (M*mu)/sqrt(M*sigma.^2 +M*(M-1)*r*sigma.^2) % signal to noise
pooledSNR/(mu/sigma) % <--- should be 1.5 to match monkey, should be 4.5 to match cones

% Assymptotic value of pooledSNR
pooledSNR_asym = mu./sqrt(r*sigma.^2)

% Do these calculations really depend on the specific values selected for
% mu and sigma?
% Let's say each V1 neuron has a d' of 5, then the monkey's d' is 1.6*5 = 8
mu = .5;
sigma = 1;
M = 10;
r =0.2;

pooledSNR = (M*mu)/sqrt(M*sigma.^2 +M*(M-1)*r*sigma.^2) % signal to noise

[pooledSNR 1.6*(mu/sigma)]
% The pool of 10 neurons (with r = .2) is always 1.9x more sensitive than the
% single neuron and 1.2x more sensitive than the photocurrent ideal
% observer.

%%
% Section 5)
% Working on Figure 4 (Potpourri figure) for the modeling paper with Juan and Fred

load('/Users/greghorwitz/Documents/Manuscripts/Charlie''s model/Figure 4 Material/fig4.mat');
%%
% First panel
ylims = [0.0035 0.6097];
mn = fig4.a.kali.threshold;
sem = fig4.a.kali.sem;
modeldat = fig4.a.model.threshold;
sf = fig4.a.SF;
colordirs = fig4.a.colors;
plotcols = [241 32 42; 0 0 0; 0 0 0; 53 85 160; 112 207 222; 181 83 156]/255; 
% L-M, L+M, Ach, S, L-M-S, L-M+S
% rows = colors, cols = SFs
figure;
axes; hold on;
for i = [1 3 4 5 6] % skipping L+M
    for j = 1:length(sf)
        h = plot([sf(j) sf(j)], mn(i,j)+sem(i,j)*[-1 1],'-','LineWidth',2);
        set(h,'MarkerFaceColor',plotcols(i,:),'MarkerEdgeColor',plotcols(i,:),'Color',plotcols(i,:))
    end  
    h(1) = plot(sf, mn(i,:),'-','LineWidth',2);
    h(2) = plot(sf, mn(i,:),'o','LineWidth',1);
    set(h,'MarkerFaceColor',plotcols(i,:),'MarkerEdgeColor',plotcols(i,:),'Color',plotcols(i,:))
end
set(gca,'Xscale','log','YScale','log')
set(gca,'Ylim',ylims);
% Model
figure;
axes; hold on;
for i = [1 3 4 5 6] % skipping L+M
    h(1) = plot(sf, modeldat(i,:),'-','LineWidth',2);
    h(2) = plot(sf, modeldat(i,:),'o','LineWidth',1);
    set(h,'MarkerFaceColor',plotcols(i,:),'MarkerEdgeColor',plotcols(i,:),'Color',plotcols(i,:))   
end
set(gca,'Xscale','log','YScale','log')
set(gca,'Ylim',10.^(log10(ylims)-log10(ylims(1))-3));

% Now a plot with both sets of data on the same axes
figure; axes; hold on;
for i = [1 3 4 5 6] % skipping L+M
    h(1) = plot(sf, mn(i,:),'-','LineWidth',2);
    h(2) = plot(sf, mn(i,:),'o','LineWidth',1);
    set(h([1 2]),'MarkerFaceColor',plotcols(i,:),'MarkerEdgeColor',plotcols(i,:),'Color',plotcols(i,:))   

    h(3) = plot(sf, modeldat(i,:),':','LineWidth',2);
    h(4) = plot(sf, modeldat(i,:),'o','LineWidth',1);
    set(h([3 4]),'MarkerFaceColor',plotcols(i,:),'MarkerEdgeColor',plotcols(i,:),'Color',plotcols(i,:))   

end
set(gca,'Xscale','log','YScale','log')
set(gca,'Ylim',[.002, .5]);

% What is the S-cone contrast at threshold for the intermediates at 2 cpd?
fig4.a.kali.threshold(end-1,3)*fig4.a.colors(end-1,:)
fig4.a.kali.threshold(end,3)*fig4.a.colors(end,:)

%%
% Second panel
ylims = [0.009  0.2258];
tf = fig4.b.TF;
modeldat = fig4.b.model.thresh;
mn = fig4.b.nut.thresh;
sem = fig4.b.nut.sem;
whichtfs = [2:length(tf)-1];

figure; axes; hold on;
for i = 1:size(mn,1)
    for j = whichtfs % Skipping first two TFs
        h = plot([tf(j) tf(j)], mn(i,j)+sem(i,j)*[-1 1],'-','LineWidth',2);
        set(h,'MarkerFaceColor',plotcols(i,:),'MarkerEdgeColor',plotcols(i,:),'Color',plotcols(i,:))
    end  
    h(1) = plot(tf(whichtfs), mn(i,whichtfs),'-','LineWidth',2);
    h(2) = plot(tf(whichtfs), mn(i,whichtfs),'o','LineWidth',1);
    set(h,'MarkerFaceColor',plotcols(i,:),'MarkerEdgeColor',plotcols(i,:),'Color',plotcols(i,:));
end
set(gca,'Xscale','log','YScale','log')
set(gca,'Ylim',ylims);
set(gca,'Xlim',[.9*tf(whichtfs(1)) tf(whichtfs(end))*1.1])

% % Plotting DC
% figure; axes; hold on;
% for i = 1:size(mn,1)
%     h = plot([tf(1) tf(1)], mn(i,1)+sem(i,1)*[-1 1],'-','LineWidth',2);
%     set(h,'MarkerFaceColor',plotcols(i,:),'MarkerEdgeColor',plotcols(i,:),'Color',plotcols(i,:))
%     h = plot(tf(1), mn(i,1),'o','LineWidth',1);
%     set(h,'MarkerFaceColor',plotcols(i,:),'MarkerEdgeColor',plotcols(i,:),'Color',plotcols(i,:));
% end
% set(gca,'YScale','log');
% set(gca,'Ylim',ylims);
% 
% % Plotting model thresholds
% figure; axes; hold on;
% h(1) = plot(tf(whichtfs), modeldat(whichtfs),'-','LineWidth',2);
% h(2) = plot(tf(whichtfs), modeldat(whichtfs),'o','LineWidth',1);
% set(h,'MarkerFaceColor','blue','MarkerEdgeColor','blue','Color','blue')
% set(gca,'YScale','log','Xscale','log')
% set(gca,'Ylim',10.^(log10(ylims)-log10(ylims(1))-2.7));
% set(gca,'Xlim',[.9*tf(whichtfs(1)) tf(whichtfs(end))*1.1])
% 
% 
% % model
% figure; axes; hold on;
% h = plot(tf(1), modeldat(1),'o','LineWidth',1);
% set(h,'MarkerFaceColor','blue','MarkerEdgeColor','blue','Color','blue');
% set(gca,'YScale','log')
% set(gca,'Ylim',10.^(log10(ylims)-log10(ylims(1))-2.7));

% Now plotting it all on common axes


figure; axes; hold on;
for i = 1:size(mn,1)
    for j = whichtfs % Skipping first two TFs
        h = plot([tf(j) tf(j)], mn(i,j)+sem(i,j)*[-1 1],'-','LineWidth',2);
        set(h,'MarkerFaceColor',plotcols(i,:),'MarkerEdgeColor',plotcols(i,:),'Color',plotcols(i,:))
    end  
    h(1) = plot(tf(whichtfs), mn(i,whichtfs),'-','LineWidth',2);
    h(2) = plot(tf(whichtfs), mn(i,whichtfs),'o','LineWidth',1);
    set(h,'MarkerFaceColor',plotcols(i,:),'MarkerEdgeColor',plotcols(i,:),'Color',plotcols(i,:));
end
set(gca,'Xscale','log','YScale','log')
set(gca,'ylim',[.002 .1])
set(gca,'Xlim',[.9*tf(whichtfs(1)) tf(whichtfs(end))*1.1])


% Plotting DC
figure; axes; hold on;
for i = 1:size(mn,1)
    h = plot([tf(1) tf(1)], mn(i,1)+sem(i,1)*[-1 1],'-','LineWidth',2);
    set(h,'MarkerFaceColor',plotcols(i,:),'MarkerEdgeColor',plotcols(i,:),'Color',plotcols(i,:))
    h = plot(tf(1), mn(i,1),'o','LineWidth',1);
    set(h,'MarkerFaceColor',plotcols(i,:),'MarkerEdgeColor',plotcols(i,:),'Color',plotcols(i,:));
end
set(gca,'YScale','log');
set(gca,'ylim',[.002 .1])

% Plotting model thresholds
figure; axes; hold on;
h(1) = plot(tf(whichtfs), modeldat(whichtfs),'-','LineWidth',2);
h(2) = plot(tf(whichtfs), modeldat(whichtfs),'o','LineWidth',1);
set(h,'MarkerFaceColor','blue','MarkerEdgeColor','blue','Color','blue')
set(gca,'YScale','log','Xscale','log')
set(gca,'ylim',[.002 .1])
set(gca,'Xlim',[.9*tf(whichtfs(1)) tf(whichtfs(end))*1.1])

% model
figure; axes; hold on;
h = plot(tf(1), modeldat(1),'o','LineWidth',1);
set(h,'MarkerFaceColor','blue','MarkerEdgeColor','blue','Color','blue');
set(gca,'YScale','log')
set(gca,'ylim',[.002 .1])



%%
% Panel 3
%ylims = [0.02 0.7];
ylims = [0.006 0.7];

modeldat = [fig4.c.model_15Hz_green; fig4.c.model_15Hz_blue; fig4.c.model_25Hz_green; fig4.c.model_25Hz_blue];
mn = [fig4.c.apollo_15Hz_green; fig4.c.apollo_15Hz_blue; fig4.c.apollo_25Hz_green; fig4.c.apollo_25Hz_blue];
ecc = fig4.c.apollo_eccentricity;
ecc_mod = fig4.c.model_eccentricity;
plotcols = [13 129 71; 53 85 160; 140 204 140; 205 202 227]/255;

figure; axes; hold on;
for i = 1:size(mn,1)
    h = plot(ecc, mn(i,:),'-','LineWidth',2);
    set(h,'MarkerFaceColor',plotcols(i,:),'MarkerEdgeColor',plotcols(i,:),'Color',plotcols(i,:));
end
set(gca,'YScale','log');
set(gca,'Ylim',ylims);
set(gca,'Xlim',[0 7.5]);

% model
%figure; axes; hold on;
for i = 1:size(mn,1)
    h = plot(ecc_mod, modeldat(i,:),'--','LineWidth',2);
    set(h,'MarkerFaceColor',plotcols(i,:),'MarkerEdgeColor',plotcols(i,:),'Color',plotcols(i,:));
end
set(gca,'YScale','log');
%set(gca,'Ylim',10.^(log10(ylims)-log10(ylims(1))-2.2));
set(gca,'Xlim',[0 7.5]);

%%
% Section 6)
% Is it true that the dot product of a modulated poisson process onto a
% function is a Poisson random variable? No, duh, but its approximately 
% Gaussian with mean and variance that can be computed and a bin by bin
% basis and combined

niter = 1000;
x = linspace(0,10,100);
f1 = (sin(x)+1)/100;
f1= abs(normrnd(0,1,1,length(f1)))
f2 = sin(x*10)+1;
f2 = normrnd(0,1,1,length(f2))
%f2 = ones(1,length(f1));
%f2 = f2./sum(f2).*length(f1);

metadata = [];
for j = 1:30
    data = [];
    for i = 1:niter
        spikecounts = poissrnd(f1);
        data = [data; f2*spikecounts'];
    end
%    figure; axes; hist(data);
%    title(num2str([mean(data) var(data)]));
    metadata(j,:) = [mean(data) var(data)];
end
mn = f2*f1'
v = f2.^2*f1';

figure; axes; hold on;
plot(metadata(:,1),metadata(:,2),'k.')
plot(mn,v,'m*')
%plot([min(metadata(:)), max(metadata(:))], [min(metadata(:)), max(metadata(:))])

%%
% Section 7
% Trying to generate an isodetection surface for an ideal observer with
% access to Poisson photoisomerizations.

% cones per mm^2
ecc = sqrt((-5)^2 + (-3.5)^2); % in deg
mmperdeg = .22;
conespermm2 = 150.9*10^3*exp(-1.2*ecc)+35.9*10^3*exp(-.16*ecc)+9.9*10^3*exp(-.03*ecc); % actually cones/mm^2 ~27,000
Sconespermm2 = 2.5*10^3 *exp(-.2*ecc)+1.8*10^3*exp(-.05*ecc);

conespermm2 = conespermm2*2; % Two eyes
Sconespermm2 = Sconespermm2*2;
gaborsddeg = 0.4;
gaborsdmm = gaborsddeg*mmperdeg;
pixelsizeindeg = 0.0194; % This doesn't matter very much - make it 0.2
pixelsizeindeg = 0.2;
gaborlim = 3; % sds
npix = round(2*gaborlim*gaborsddeg/pixelsizeindeg);
sperframe = 1/75;

pixelsizeinmm = pixelsizeindeg*mmperdeg;
conesperpixel = pixelsizeinmm.^2*conespermm2;
Sconesperpixel = conesperpixel*(Sconespermm2/conespermm2);
Lconesperpixel = (conesperpixel-Sconesperpixel)/2;
Mconesperpixel = Lconesperpixel;

% How many cones inside 1 SD?
gaborsddeg = 0.15;
gaborsdmm = gaborsddeg*mmperdeg;
A = pi*gaborsdmm^2 % mm^2
S = A*Sconespermm2
M = ((A*conespermm2)-S)/2
totalcones = round(S+M*2)
sqrt(totalcones*M/totalcones*(1-M/totalcones))

sqrt(totalcones*S/totalcones*(1-S/totalcones))

% making the stimulus/weighting function
envelope = [linspace(0,1,12) ones(1,26) linspace(1,0,12)]; % in frames
nframes = length(envelope);
sinusoid = sin(2*pi*3*linspace(0,nframes*sperframe,nframes)); % sin phase provides zero DC
temporalwtframes = envelope.*sinusoid;
spatialwtpix1D = normpdf(linspace(-gaborlim*gaborsdmm,gaborlim*gaborsdmm,npix),0,gaborsdmm);
spatialwtpix1D = spatialwtpix1D./max(spatialwtpix1D);
spatialwtpix2D = spatialwtpix1D'*spatialwtpix1D;

% assuming identical temporal phase at each location
% Right now stimulus is specified in pixels/frames
% Now computing mean R*s

bkgndRstarLMS = [7131 6017 1973]; % per cone per second
spacetimewt = repmat(spatialwtpix2D,[1 1 nframes]).*repmat(permute(temporalwtframes,[1 3 2]),[npix npix 1]);

% For debugging
%spacetimewt = ones(size(spacetimewt));
% spacetimewt(12,12,25) = 1;
% End of debugging code

% Getting the noise (no Gabor) distribution
LMScc = [0 0 0]; % cone contrast at peak
L = bkgndRstarLMS(1)*(1+LMScc(1).*spacetimewt)*Lconesperpixel*sperframe; % L isomerizations each pixel, each frame 
M = bkgndRstarLMS(2)*(1+LMScc(2).*spacetimewt)*Mconesperpixel*sperframe; % M isomerizations each pixel, each frame 
S = bkgndRstarLMS(3)*(1+LMScc(3).*spacetimewt)*Sconesperpixel*sperframe; % S isomerizations each pixel, each frame 

Lmn = sum(sum(sum(L.*spacetimewt)));
Lsd = sqrt(sum(sum(sum(L.*spacetimewt.^2)))); % Poisson assumption: var(spacetimewt*X(lambda)) = spacetimewt^2*lambda
Mmn = sum(sum(sum(M.*spacetimewt)));
Msd = sqrt(sum(sum(sum(M.*spacetimewt.^2)))); 
Smn = sum(sum(sum(S.*spacetimewt)));
Ssd = sqrt(sum(sum(sum(S.*spacetimewt.^2))));
noisedist = [Lmn Mmn Smn; Lsd Msd Ssd];
% Finished with noise distribution

% brute force discriminabilty measurements
[lcc,mcc,scc] = meshgrid(linspace(-.004,.004,80),linspace(-.004,.004,80),linspace(-.01,.01,80));
LMScc = [lcc(:), mcc(:), scc(:)]; % cone contrast at peak
data = zeros(size(LMScc,1),1);
for i = 1:size(LMScc,1)
    L = bkgndRstarLMS(1)*(1+LMScc(i,1).*spacetimewt)*Lconesperpixel*sperframe; % L isomerizations each pixel, each frame
    M = bkgndRstarLMS(2)*(1+LMScc(i,2).*spacetimewt)*Mconesperpixel*sperframe; % M isomerizations each pixel, each frame
    S = bkgndRstarLMS(3)*(1+LMScc(i,3).*spacetimewt)*Sconesperpixel*sperframe; % S isomerizations each pixel, each frame
    
    Lmn = sum(sum(sum(L.*spacetimewt)));
    Lsd = sqrt(sum(sum(sum(L.*spacetimewt.^2)))); % Poisson assumption: var(spacetimewt*X(lambda)) = spacetimewt^2*lambda
    Mmn = sum(sum(sum(M.*spacetimewt)));
    Msd = sqrt(sum(sum(sum(M.*spacetimewt.^2))));
    Smn = sum(sum(sum(S.*spacetimewt)));
    Ssd = sqrt(sum(sum(sum(S.*spacetimewt.^2))));
    
    signaldist = [Lmn Mmn Smn; Lsd Msd Ssd];
    
    v = (signaldist(1,:)-noisedist(1,:))./signaldist(2,:).^2; % linear discriminant vector (is cone contrast vector ?!)
    v = v./norm(v); % <--- This step is important and isn't mentioned in the text
    mu1 = signaldist(1,:)*v';
    mu2 = noisedist(1,:)*v';
    sigmasq = signaldist(2,:).^2*(v.^2)';
    data(i) = (mu1-mu2)./sqrt(sigmasq);
    
    % This is equivalent to above 
    % dprime = (signaldist(1,:)-noisedist(1,:))./signaldist(2,:); % d-prime
    % data(i) = sqrt(dprime*dprime');
end
fn = reshape(data,size(lcc));
figure;
isosurfstruct = isosurface(lcc,mcc,scc,fn,1.25);
h = patch(isosurfstruct);
set(h,'FaceColor',[0 .5 0],'edgecolor','none');
set(h,'FaceLighting','gouraud');
set(gca,'View',[30 36])
h_light = lightangle(45,45);
set(h_light,'position',[4 -1.5 -.5])
axis equal
set(gca,'Visible','off');
set(gcf,'Color',[1 1 1]);
export_fig junk3 -tiff -m4

%%
% d-prime of what correponds to AUC of 0.816? A: 1.25
% dprime = 1.25;
% x = linspace(-5,5,100);
% a = normpdf(x,0,1);
% a = a./sum(a)
% b = 1-normcdf(x,dprime,1);
% auc = sum(a.*b)

%%
% Lcc = 0.0018;
% deltaL = Lcc*bkgndRstarLMS(1); % is R*/cone/sec
% elementaldprime = ((bkgndRstarLMS(1)*sperframe)-((bkgndRstarLMS(1)+deltaL)*sperframe))/sqrt(bkgndRstarLMS(1)*sperframe);
% sqrt(sum(sum(sum(Lconesperpixel.*(spacetimewt.*elementaldprime).^2))))
% % Adding up all the tiny d's
% % Justification for multipling by spacetimewt in the pooling of d's:
% % Scaling the signal by a factor of 'x' scales d' by a factor of 'x'.
% 
% % Now considering each pixel to be "elemental"
% elementaldprime = ((bkgndRstarLMS(1)*sperframe*Lconesperpixel)-((bkgndRstarLMS(1)+deltaL)*sperframe*Lconesperpixel))/sqrt(bkgndRstarLMS(1)*sperframe*Lconesperpixel);
% sqrt(sum(sum(sum((spacetimewt.*elementaldprime).^2))))
% 
% % % A very simple (but illuminating) toy problem
% % mu1 = 0
% % mu2 = 1;
% % sigmasq = 1;
% % n = 2;
% % % way 1:
% % elementaldprime = (mu2-mu1)/sqrt(sigmasq);
% % sqrt(n*elementaldprime.^2)
% % % way 2:
% % ((n*mu2)-(n*mu1))/sqrt(n*sigmasq)
% 

%%
% Section 8)
% Direct comparison between Charlie's model and the Poisson model
% need to load "modIsoSurf.mat" and run the cell (two) above.
modIsoSurf.thresholds
viewSetting = [172.2 10.161]; % GDLH: Want to make surface same orientation as heat map
lightattenuationfactor = .5;

scales = log10(logspace(0,.2,20));
% I need to find thresholds in these directions
thresholds = nan*zeros(size(modIsoSurf.colorDirs,1),1);
for i = 1:size(modIsoSurf.colorDirs,1)
    tmpdata = [];
    for j = 1:length(scales)
        
        L = bkgndRstarLMS(1)*(1+modIsoSurf.colorDirs(i,1)*scales(j).*spacetimewt)*Lconesperpixel*sperframe;
        M = bkgndRstarLMS(2)*(1+modIsoSurf.colorDirs(i,2)*scales(j).*spacetimewt)*Mconesperpixel*sperframe;
        S = bkgndRstarLMS(3)*(1+modIsoSurf.colorDirs(i,3)*scales(j).*spacetimewt)*Sconesperpixel*sperframe;
        
        Lmn = sum(sum(sum(L.*spacetimewt)));
        Lsd = sqrt(sum(sum(sum(L.*spacetimewt.^2)))); % Poisson assumption: var(spacetimewt*X(lambda)) = spacetimewt^2*lambda
        Mmn = sum(sum(sum(M.*spacetimewt)));
        Msd = sqrt(sum(sum(sum(M.*spacetimewt.^2))));
        Smn = sum(sum(sum(S.*spacetimewt)));
        Ssd = sqrt(sum(sum(sum(S.*spacetimewt.^2))));
        
        signaldist = [Lmn Mmn Smn; Lsd Msd Ssd];
        if (j > 1)
            v = (signaldist(1,:)-noisedist(1,:))./signaldist(2,:).^2; % linear discriminant vector (is cone contrast vector ?!)
            v = v./norm(v);
            mu1 = signaldist(1,:)*v';
            mu2 = noisedist(1,:)*v';
            sigmasq = signaldist(2,:).^2*(v.^2)';
            tmpdata(j) = (mu1-mu2)./sqrt(sigmasq);
        else
            tmpdata(j) = 0;
            noisedist = signaldist;
        end
    end
    thresholds(i) = interp1(tmpdata, scales, 1.25);
end
[modIsoSurf.thresholds thresholds modIsoSurf.thresholds./thresholds]

% Plotting isothreshold contour for Poisson catch model
tmp = modIsoSurf.colorDirs.*repmat(thresholds,[1 3]);
% getting radii
D = [tmp(:,1) .* tmp(:,1),...
    tmp(:,2) .* tmp(:,2),...
    tmp(:,3) .* tmp(:,3),...
    2*tmp(:,1) .* tmp(:,2),...
    2*tmp(:,1) .* tmp(:,3),...
    2*tmp(:,2) .* tmp(:,3)];
v = (D' * D) \(D' * ones(size(tmp,1),1));
A = [v(1) v(4) v(5);...
    v(4) v(2) v(6);...
    v(5) v(6) v(3)];
[evecs, evals] = eig(A);
evals = diag(evals);
evals = evals([3 2 1]);  % back in LMS order. Fragile
radii(:,1) = sqrt(1./evals);
[lcc,mcc,scc] = meshgrid(linspace(-radii(1,1),radii(1,1),40),...
    linspace(-radii(2,1),radii(2,1),40),...
    linspace(-radii(3,1),radii(3,1),40));
fn = (lcc./radii(1,1)).^2+(mcc./radii(2,1)).^2+(scc./radii(3,1)).^2;

% Plotting isothreshold contour for photon catch model
figure('position', [256 5 1150 824]);
subplot(1,2,1); hold on; % behavioral data
p = patch(isosurface(lcc,mcc,scc,fn,1));
set(p,'edgecolor','none','facecolor',[0 .6 0]);
set(p,'SpecularExponent',1)
lighting gouraud;
set(gca,'View',viewSetting); axis equal;
set(gca,'Xlim',[-.0025 .0025],'Ylim',[-.0025 .0025],'Zlim',[-.01 .01]);
xlabel('L'); ylabel('M'); zlabel('S'); 
title('Poisson quantal catch model');
h_light = camlight(20, 20);
set(h_light,'Color',[1 1 1]*lightattenuationfactor);
%h_light = light('Position',[-1 1 .6]);
%h_light = lightangle(viewSetting(1),viewSetting(2));
%plot3(tmp(:,1),tmp(:,2),tmp(:,3),'ko','MarkerFaceColor',[0 .4 0],'MarkerEdgeColor','none','markersize', 4);
%plot3(-tmp(:,1),-tmp(:,2),-tmp(:,3),'ko','MarkerFaceColor',[0 .4 0],'MarkerEdgeColor','none','markersize', 4);
%set(gcf,'renderer','zbuffer'); 
set(gca,'Visible','off');
set(gcf,'Color',[1 1 1]);
%export_fig junk_poisson -tiff -m3 -zbuffer


% Plotting isothreshold contour for cone noise model
figure('position', [256 5 1150 824]);
subplot(1,2,1); hold on;% behavioral data
tmp = modIsoSurf.colorDirs.*repmat(modIsoSurf.thresholds,[1 3]);
% getting radii
D = [tmp(:,1) .* tmp(:,1),...
    tmp(:,2) .* tmp(:,2),...
    tmp(:,3) .* tmp(:,3),...
    2*tmp(:,1) .* tmp(:,2),...
    2*tmp(:,1) .* tmp(:,3),...
    2*tmp(:,2) .* tmp(:,3)];
v = (D' * D) \(D' * ones(size(tmp,1),1));
A = [v(1) v(4) v(5);...
    v(4) v(2) v(6);...
    v(5) v(6) v(3)];
[evecs, evals] = eig(A);
evals = diag(evals);
evals = evals([3 2 1]);  % back in LMS order. Fragile
radii(:,2) = sqrt(1./evals);
[lcc,mcc,scc] = meshgrid(linspace(-radii(1,2),radii(1,2),40),...
    linspace(-radii(2,2),radii(2,2),40),...
    linspace(-radii(3,2),radii(3,2),40));
fn = (lcc./radii(1,2)).^2+(mcc./radii(2,2)).^2+(scc./radii(3,2)).^2;
p = patch(isosurface(lcc,mcc,scc,fn,1));
set(p,'edgecolor','none','facecolor',[.0 .6 .0])
set(p,'SpecularExponent',1)
lighting gouraud;
set(gca,'View',viewSetting);
axis equal;
set(gca,'Xlim',[-.0025 .0025],'Ylim',[-.0025 .0025],'Zlim',[-.01 .01]);
title('Cone noise model');
xlabel('L'); ylabel('M'); zlabel('S');
h_light = camlight(20, 20);
set(h_light,'Color',[1 1 1]*lightattenuationfactor)
% h_light = light('Position',[-1 1 .6]);

%h_light = lightangle(viewSetting(1),viewSetting(2));
%plot3(tmp(:,1),tmp(:,2),tmp(:,3),'ko','MarkerFaceColor',[0 .4 0],'MarkerEdgeColor','none','markersize', 4);
%plot3(-tmp(:,1),-tmp(:,2),-tmp(:,3),'ko','MarkerFaceColor',[0 .4 0],'MarkerEdgeColor','none','markersize', 4);
% for surface and points
%set(gcf,'renderer','zbuffer'); set(gca,'Visible','off'); export_fig junk_cone -tiff -m2
set(gca,'Visible','off');
set(gcf,'Color',[1 1 1]);
%export_fig junk_cone -tiff -m3 -zbuffer

%set(gcf,'renderer','painters') % for axes

% Plotting a heat map of the ratio
figure; axes; hold on;
[X,Y,Z] = sphere(800);
colors = [X(:), Y(:), Z(:)];
colors = bsxfun(@rdivide, colors, sqrt(sum(colors.^2, 2)));
poissThresh = coleThresh(diag(1./radii(:,1)), 2, colors);
conesThresh = coleThresh(diag(1./radii(:,2)), 2, colors);
C = conesThresh./poissThresh;
C = reshape(C, size(X,1), size(X,2)); % for plotting
set(gcf, 'position', [256 5 1150 824]);
hand = surf(X,Y,Z,C);
set(hand, 'edgealpha', 0);
set(gca,'View',viewSetting);

axis square;
title('Cone noise thresh/Poiss thresh');
%colorbar;
colormap(jet(255));
caxis(gca,[1,2.5]);
xlabel('L'); ylabel('M'); zlabel('S');
set(gca,'Visible','off');
%export_fig junk -tiff -m3 -zbuffer

%%
% Section 9
% Trying to render a DTNT surface with surface shading
% Now using a file generated by "coneNoiseAnalysis.m" that contains the
% threshold measurments. The structure "dtnt" is in the files
% sedna_DTNT_022815.mat and kali_DTNT_022815.mat
INHERIT_fpar_cones = 1;
MAKEICON  = 1;

defaultViewSetting = [172.2 10.161]; % Want to make surface same orientation as heat map
% load('/Users/greghorwitz/Documents/Manuscripts/Charlie''s model/For Greg/kali_DTNT_022815.mat'); % Might need to do this by hand if file is not in the path
monkeyAlphas = cellfun(@(x,y) x./y, dtnt.alpha',  repmat({100}, numel(dtnt.alpha),1), 'uniformoutput', 0);

sfIdx = logical([1;0;0;0])  % which SF should be analyzed? [0.5 1 2 4]
sfs = unique(cat(1,dtnt.sfs{:}))

% make a logical vector to weed out conditions (if need be)
l_questBeta = dtnt.l_beta > 0;
l_bkgndrgb = dtnt.l_bkgndrgb > 0;

% pull out the spatial frequncy condition
l_sfs = cat(1,dtnt.sfs{:}) == sfs(sfIdx);

% grab the appropriate files
l_valid = [l_sfs & l_questBeta & l_bkgndrgb];
colors = cat(1, dtnt.colorDirs{l_valid});
alphas = cat(1, monkeyAlphas{l_valid});

% default init params on the monkey
[fpar_monkey_default, fval_monkey_default] = fitDetectionSurface(colors, alphas, 'ellipsoid');

% init params from kali at 0.5 cpd (on the monkey)
initparams_kali = [1.15; -15.82; 38.29; -14.36; -17.29; 6.62; -8.74; -81.82; 68.40; 9.36];  % optional input argument for threshSurfPlot for 0.5cpd
[fpar_monkey_kali, fval_monkey_kali] = fitDetectionSurface(colors, alphas, initparams_kali);

% particle swarm on the monkey
%npsoiters = 0;
%[fpar_monkey_pso, fval_monkey_pso] = fitDetectionSurface(colors, alphas, 'pso', npsoiters);

% find the best method
%method = {'default'; 'kali'; 'pso'};
%fvals = [fval_monkey_default; fval_monkey_kali; fval_monkey_pso]
%fpars = [fpar_monkey_default'; fpar_monkey_kali'; fpar_monkey_pso];

method = {'default'; 'kali'};
fvals = [fval_monkey_default; fval_monkey_kali]
fpars = [fpar_monkey_default'; fpar_monkey_kali'];

[~, idx] = min(fvals);
bestMethod = method{idx}
fpar_monkey = fpars(idx,:);

% fit the cone noise model's data
if (~INHERIT_fpar_cones)
    fpar_cones = fitDetectionSurface(modIsoSurf.colorDirs, modIsoSurf.thresholds, 'ellipsoid');
end

% plot the raw data
data = bsxfun(@times, colors, alphas);
lims = max(abs(data(:,[1 2 3])))*1.1;
npts = 150;
[xx yy zz] = meshgrid(linspace(-lims(1),lims(1),npts),...
    linspace(-lims(2),lims(2),npts),...
    linspace(-lims(3),lims(3),npts));
xformedxyz = [xx(:) yy(:) zz(:)];
Loog = logical(data(:,end) == 1);

fr = sum(abs(xformedxyz *reshape(fpar_monkey(2:end),3,3)).^fpar_monkey(1),2);
surfstruct = isosurface(xx,yy,zz,reshape(fr,size(xx)), 1);

% Now doing the plotting
% 3-D renderings with specularities
for i = 1:2 % two orientations
    if i == 1
        viewSetting = defaultViewSetting;
    else
        viewSetting = defaultViewSetting+[-90 0];
    end
    
    figure; axes; hold on;
    p = patch(surfstruct);
    set(p,'EdgeColor', 'none', 'FaceAlpha',.8,'FaceColor',[0 .6 0],'Edgealpha',0);
    set(p,'SpecularExponent',1)
    lighting gouraud;
    h = [];
    for k = [-1, 1]
        h = [h; plot3(k*data(~Loog,1),k*data(~Loog,2),k*data(~Loog,3),'ko')];
    end
    set(h,'MarkerFaceColor',[0 .25 0],'MarkerSize',4,'MarkerEdgeColor','none');
    view(viewSetting);
    %h_light=light('Position',lightposition);
    h_light = camlight(20, 20);
   % h_light=lightangle(viewSetting(1)-14, viewSetting(2)-5);
    
    xlabel('L'); ylabel('M'); zlabel('S');
    axis equal;
    axis vis3d
    set(gca,'Xlim',[-lims(1) lims(1)],'Ylim',[-lims(2) lims(2)],'Zlim',[-lims(3) lims(3)]);
    
   % cd /Volumes/Macintosh' HD'/Users/greghorwitz/Desktop/
   % set(gcf,'renderer','painters') % for axes
   % eval(['print -depsc ',['junk',num2str(i)]]);

    set(gca,'Visible','off');% for surface and points
    set(gcf,'Color',[1 1 1]);
%    eval(['export_fig ',['junk',num2str(i),'1'],' -tiff -m5 -zbuffer;']);
    eval(['export_fig ',['junk',num2str(i),'1'],' -tiff -m3;']);
end

%
% make a map of the ratio b/w Monkey and Retina. Where the ratio is large,
% the monkey does worse than the retina....
%
%%%%%%

[X,Y,Z] = sphere(400);
spherecolors = [X(:), Y(:), Z(:)];
spherecolors = bsxfun(@rdivide, spherecolors, sqrt(sum(spherecolors.^2, 2)));
monkeyThresh = coleThresh(reshape(fpar_monkey(2:end),3,3), fpar_monkey(1), spherecolors);
conesThresh = coleThresh(reshape(fpar_cones(2:end),3,3), fpar_cones(1), spherecolors);
C = monkeyThresh ./ conesThresh;
C = log10(C);
C = reshape(C, size(X,1), size(X,2)); % for plotting
minVal = min(C(:))
minVec = [X(C(:) == minVal), Y(C(:) == minVal), Z(C(:) == minVal)] .* 1.5;
if size(minVec,1) ==1; minVec = [minVec; -minVec]; end
maxVal = max(C(:))
maxVec = [X(C(:) == maxVal), Y(C(:) == maxVal), Z(C(:) == maxVal)] .* 1.5;
if size(maxVec,1) ==1; maxVec = [maxVec; -maxVec]; end

figure, hold on,
set(gcf, 'position', [256 5 1150 824]);
hand = surf(X,Y,Z,C);
set(hand, 'edgealpha', 0);
%plot3([-1.4, 1.4], [-1.4 1.4], [0,0], '-', 'linewidth', 18, 'color', 'k') % the L+M direction
%plot3([-1, 1], [-1 1], [-1 1], '-', 'linewidth', 5, 'color', 'k')           % the achromatic direction
%plot3([0, 0], [0 0], [-1.4 1.4], '-', 'linewidth', 18, 'color', 'b')         % S-iso color direction
%plot3([1.2, -1.2], [-1.2 1.2], [0 0], '-', 'linewidth', 18, 'color', 'r')            % L-M color direction
%plot3(minVec(:,1), minVec(:,2), minVec(:,3), '--', 'linewidth', 5, 'color', 'm')
%plot3(maxVec(:,1), maxVec(:,2), maxVec(:,3), '--', 'linewidth', 18, 'color', [.3 .3 .3])
xlabel('\DeltaL/L'); ylabel('\DeltaM/M'); zlabel('\DeltaS/S')
xlim([-1.3, 1.3]); ylim([-1.3, 1.3]); zlim([-1.3, 1.3]);
set(gca, 'view', [172.2 10.161],...
    'visible', 'off')
axis square
maps = {'Edge', 'CubicL', 'CubicYF', 'IsoL'};
set(gcf,'Color',[1 1 1]);

colormap(pmkmp(250, maps{2}))
%colorbar('Ytick',[minVal;maxVal],'YtickLabel',num2str([10^minVal;10^maxVal]));

if (MAKEICON)
    figure, hold on,
    set(gcf, 'position', [256 5 96 96]);
    hand = surf(X,Y,Z,C);
    set(hand, 'edgealpha', 0);
    plot3([-1.4, 1.4], [-1.4 1.4], [0,0], '-', 'linewidth', 1, 'color', 'k') % the L+M direction
    plot3([0, 0], [0 0], [-1.4 1.4], '-', 'linewidth', 1, 'color', 'b')         % S-iso color direction
    plot3([1.2, -1.2], [-1.2 1.2], [0 0], '-', 'linewidth', 1, 'color', 'r')            % L-M color direction
    set(gca, 'view', [172.2 10.161],...
        'visible', 'off')
    axis vis3d
    maps = {'Edge', 'CubicL', 'CubicYF', 'IsoL'};
    set(gcf,'Color',[1 1 1]);
    xlim([-1.3, 1.3]); ylim([-1.3, 1.3]); zlim([-1.3, 1.3]);
    
    views = linspace(0,2*pi,80)*180/pi;
    images = zeros(96,96,length(views));
    [imind,cm] = rgb2ind(im,256);
    outfile = 'HassAnguerya.gif';
    for i = 1:length(views)
        set(gca,'View',[views(i) 10]);
        
        drawnow;
        frame = getframe(gcf);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);       
        % On the first loop, create the file. In subsequent loops, append.
        if i==1
            imwrite(imind,cm,outfile,'gif','DelayTime',0,'loopcount',inf);
        else
            imwrite(imind,cm,outfile,'gif','DelayTime',0,'writemode','append');
        end
    end
end

%%
% Section 10)
% Cone mosaic stuff
% Just pooling signals from a particular cone type within the spatial
% extent of our Gaussian envelope (from the Gabor stimulus).

% cones per mm^2
ecc = sqrt((-5)^2 + (-3.5)^2); % in deg
mmperdeg = .22;
conespermm2 = 150.9*10^3*exp(-1.2*ecc)+35.9*10^3*exp(-.16*ecc)+9.9*10^3*exp(-.03*ecc); % actually cones/mm^2 ~27,000
Sconespermm2 = 2.5*10^3 *exp(-.2*ecc)+1.8*10^3*exp(-.05*ecc);

gaborsddeg = 0.15; % Should be 0.4
gaborsdmm = gaborsddeg*mmperdeg;

% Ok, we're going to make a cone mosaic that extends "n" sds in both
% directions
mm2 = .6; % size of cone array in mm2
nstds_array = sqrt(mm2)/gaborsdmm; % Size of cone array re: Gabor SD
totconesinfield = round(mm2*conespermm2);
Sconesinfield = round(mm2*Sconespermm2);
Mconesinfield = round((totconesinfield-Sconesinfield)/2);

imsize = 2000; ratio = .8; jitter = 0.1; % paramaters taken from Curcio and Sloan paper
[Clist,nCones,nL,nM,nS] = MakeSloanClist(totconesinfield,Mconesinfield,Sconesinfield,imsize,ratio,jitter);
figure; axes; hold on;
L = Clist(:,3) == 1;
plot(Clist(L,1),Clist(L,2),'ro','MarkerFaceColor','red','MarkerSize',2); % 1.5 is good for sd = 0.4°
L = Clist(:,3) == 2;
plot(Clist(L,1),Clist(L,2),'go','MarkerFaceColor','green','MarkerSize',2);
L = Clist(:,3) == 3;
plot(Clist(L,1),Clist(L,2),'bo','MarkerFaceColor','blue','MarkerSize',2);
axis square; axis equal;

% Drawing the gabor on the cone mosaic
% plot(imsize*[1/4 3/4 3/4 1/4 1/4],imsize*[1/4 1/4 3/4 3/4 1/4],'k-');

% Creating the Gaussian envelope of the Gabor
% The Gaussian function is half as wide and half as tall at the imsize x
% insize cone mosaic because we want to be able to convolve a bit without
% falling off the edge of the cone mosaic.
nstds_Gabor = 6; % Size of Gabor in SDs; should be smaller than nstds_array
if (nstds_Gabor > nstds_array)
    error('nstds_Gabor > nstds_array');
end

x = linspace(0,gaborsdmm*nstds_Gabor,imsize*(nstds_Gabor/nstds_array)); % in mm
gauss1d = normpdf(x,mean(x),gaborsdmm); % double check
gauss = gauss1d'*gauss1d;
figure;
imagesc(gauss); colormap(flipud(gray));
axis square; axis equal;
set(gca,'xlim',[0 imsize-1],'ylim',[0 imsize-1])

% Looking at S-cones
Sconematrix = zeros(imsize);
L = Clist(:,3) == 3;
Sconematrix((Clist(L,1)-1)*imsize+Clist(L,2)) = 1;
result = conv2(Sconematrix,gauss,'valid');

% Looking at M-cones
Mconematrix = zeros(imsize);
L = Clist(:,3) == 1;
Mconematrix((Clist(L,1)-1)*imsize+Clist(L,2)) = 1;
result = conv2(Mconematrix,gauss,'valid');

figure;
h = surface(result); set(gca,'View',[-23 -20]);
set(h,'EdgeAlpha',0);
set(gca,'Clim',[max(result(:))*.8 max(result(:))]);
set(gca,'Zlim',[0 max(result(:))*1.1])
[std(result(:)) mean(result(:))]
std(result(:))./mean(result(:))


%%
% Section 11)
% How many cones are within 1 or 2 standard deviations of the Gabor?

gaborcenter = [-6 0]; % deg
mmperdeg = .22;
gaborsddeg = 0.4; % Should be 0.4
gaborsdmm = gaborsddeg*mmperdeg;

% Making a grid of retinal locations
[x,y] = meshgrid([-3:.1:3]*gaborsdmm+gaborcenter(1),[-3:.1:3]*gaborsdmm+gaborcenter(2));
ecc = sqrt(x.^2+y.^2);
binwidth = y(2)-y(1); % mm/bin side
binarea = binwidth.^2;% mm^2/bin

conespermm2 = 150.9*10^3*exp(-1.2*ecc)+35.9*10^3*exp(-.16*ecc)+9.9*10^3*exp(-.03*ecc); % actually cones/mm^2 ~27,000
Sconespermm2 = 2.5*10^3 *exp(-.2*ecc)+1.8*10^3*exp(-.05*ecc);
conesperbin = conespermm2.*binarea;
Sconesperbin = Sconespermm2.*binarea;

dist2fromgaborcenter = ((x-gaborcenter(1)).^2+(y-gaborcenter(2)).^2);
mask = dist2fromgaborcenter < 2*2*gaborsdmm.^2; % first 2 is for mean of chi2

sum(sum(conesperbin.*mask))
sum(sum(Sconesperbin.*mask))

% Trying it the wrong way as a sanity check
mask = abs(x-gaborcenter(1))/gaborsdmm < 2 & abs(y-gaborcenter(2))/gaborsdmm < 2;

%%
% Section 12)
% Opponent vs non-opponent signals.
% cones per mm^2. Similar to the above cell.
ecc = sqrt((-5)^2 + (-3.5)^2); % in deg
mmperdeg = .22;
conespermm2 = 150.9*10^3*exp(-1.2*ecc)+35.9*10^3*exp(-.16*ecc)+9.9*10^3*exp(-.03*ecc); % actually cones/mm^2 ~27,000
Sconespermm2 = 2.5*10^3 *exp(-.2*ecc)+1.8*10^3*exp(-.05*ecc);

centerdiams = [0.005 0.01 0.05];
relativestrengths = linspace(0,1,20);
data = [];
Clist = [];
for i = 1:length(relativestrengths)
    for j = 1:length(centerdiams)
        centersddeg = centerdiams(j); % From Croner and Kaplan 0.05
        centersdmm = centersddeg*mmperdeg;
        surroundsddeg = 0.40; % From Croner and Kaplan 0.43
        surroundsdmm = surroundsddeg*mmperdeg;
        Kc = 1;  % 114 This is the peak height.
        Ks = Kc*(centersddeg/surroundsddeg)^2;  % Equal integral for center and surround
        % Squaring is because convolution kernel is 2D
        Ks = Ks*relativestrengths(i);
        
        % Ok, we're going to make a cone mosaic that extends "n" surround sds
        nstds_array = 8; % Size of cone array re: Gabor SD
        mm2 = surroundsdmm*nstds_array;
        totconesinfield = round(mm2*conespermm2);
        Sconesinfield = round(mm2*Sconespermm2);
        Mconesinfield = round((totconesinfield-Sconesinfield)/2);
        
        imsize = 500; ratio = .8; jitter = 0.1; % paramaters taken from Curcio and Sloan paper
        if (isempty(Clist))
            [Clist,nCones,nL,nM,nS] = MakeSloanClist(totconesinfield,Mconesinfield,Sconesinfield,imsize,ratio,jitter);
        end
        figure; axes; hold on;
        L = Clist(:,3) == 1;
        plot(Clist(L,1),Clist(L,2),'ro','MarkerFaceColor','red','MarkerSize',2); % 1.5 is good for sd = 0.4°
        L = Clist(:,3) == 2;
        plot(Clist(L,1),Clist(L,2),'go','MarkerFaceColor','green','MarkerSize',2);
        L = Clist(:,3) == 3;
        plot(Clist(L,1),Clist(L,2),'bo','MarkerFaceColor','blue','MarkerSize',2);
        axis square; axis equal;
        
        % Now doing the convolution
        nstds_RF = 6; % Size of Gabor in SDs; should be smaller than nstds_array
        if (nstds_RF > nstds_array)
            error('nstds_RF > nstds_array');
        end
        
        x = linspace(0,surroundsdmm*nstds_RF,imsize*(nstds_RF/nstds_array)); % in mm
        gauss1d = normpdf(x,mean(x),surroundsdmm); % double check
        surroundgauss = gauss1d'*gauss1d;
        surroundgauss = Ks*surroundgauss./max(surroundgauss(:));
        gauss1d = normpdf(x,mean(x),centersdmm); % double check
        centergauss = gauss1d'*gauss1d;
        centergauss = Kc*centergauss./max(centergauss(:));
        gauss = centergauss-surroundgauss;
        figure;
        plot(gauss(:));
        
        % Looking at luminance
        tmpmatrix = zeros(imsize);
        Ll = Clist(:,3) == 1;
        Lm = Clist(:,3) == 2;
        tmpmatrix((Clist(Ll | Lm,1)-1)*imsize+Clist(Ll | Lm,2)) = 1;
        LUMresult = conv2(tmpmatrix,gauss,'valid');
        tmpmatrix((Clist(Ll,1)-1)*imsize+Clist(Ll,2)) = 1;
        tmpmatrix((Clist(Lm,1)-1)*imsize+Clist(Lm,2)) = -1;
        RGresult = conv2(tmpmatrix,gauss,'valid');
        
        for k= 1:2
            if (k == 1);
                stat = LUMresult;
            else
                stat = RGresult;
            end
            figure; subplot(2,1,1);
            
            h = surface(stat); set(gca,'View',[-23 -20]);
            set(h,'EdgeAlpha',0);
            subplot(2,1,2);
            hist(stat(:));
            title(['mean: ',num2str(mean(abs(stat(:))))]);
        end
        data = [data; relativestrengths(i) centerdiams(j) mean(abs(LUMresult(:))) mean(abs(RGresult(:)))]
        close all
    end
end
colors = [1 0 0; 0 .5 0; 0 0 0];
figure; axes; hold on;
for j = 1:length(centerdiams)
    L = data(:,2) == centerdiams(j);
    plot(relativestrengths,data(L,3)./data(L,4),'o-','Linewidth',2,'Color',colors(j,:),'MarkerSize',4)
end
set(gca,'Yscale','log')
xlabel('surround/center');
ylabel('LUM sens/RG sens.');
% A factor of 3 difference in our data
plot([0 1],[.5 .5],'k--')

%%
% Section 13)
% Code for making 3D rendered behavioral detectability surfaces and heat
% maps. This is for a new Figure 4 for the cone noise model manuscript.
%
% This stuff was largely taken from "coneNoiseAnalysis.m" line 255, the
% section titled "ISO-DETECTION SURFACES FOR MONKEY AND RETINA"
% cd /Users/greghorwitz/Documents/Manuscripts/Charlie''''s' model'/For' Greg'/

observers = {'kali'};         % Kali_DTNT_022815.mat or Sedna_DTNT_022815.mat
load modIsoSurf;
npsoiters = 1
viewSetting = [172.2 10.161];
FIGURE_5 = 0;  % Set to '1' is you want to calculate things needed for Figure 5;
fig5data = [];
for observeridx = 1
    observer = observers{observeridx};
    load ([observer,'_DTNT_022815.mat']);
    for sfidx = 1
        clear fpar_monkey
        % divide the behavioral data by 100 so that % is b/w 0&1.
        monkeyAlphas = cellfun(@(x,y) x./y, dtnt.alpha',  repmat({100}, numel(dtnt.alpha),1), 'uniformoutput', 0);
        
        % which SF condition should be analyzed? The behvioral data has many SFs,
        % but the model data is only at one SF.
        sfs = unique(cat(1,dtnt.sfs{:}))
        
        % make a logical vector to weed out conditions (if need be)
        l_questBeta = dtnt.l_beta > 0;
        l_bkgndrgb = dtnt.l_bkgndrgb > 0;
        
        % pull out the spatial frequncy condition
        l_sfs = cat(1,dtnt.sfs{:}) == sfs(sfidx);
        
        % grab the appropriate files
        l_valid = [l_sfs & l_questBeta & l_bkgndrgb];
        colors = cat(1, dtnt.colorDirs{l_valid});
        alphas = cat(1, monkeyAlphas{l_valid});
        
        % Not sure if this should be commented out or not
        %l_goodQuestThresh = cat(1, dtnt.validThresh{l_valid});
        %colors = colors(l_goodQuestThresh,:);
        %alphas = alphas(l_goodQuestThresh);
        
        % default init params on the monkey
        [fpar_monkey_default, fval_monkey_default] = fitDetectionSurface(colors, alphas, 'ellipsoid');
        
        % init params from kali at 0.5 cpd (on the monkey)
        initparams_kali = [1.15; -15.82; 38.29; -14.36; -17.29; 6.62; -8.74; -81.82; 68.40; 9.36];  % optional input argument for threshSurfPlot for 0.5cpd
        [fpar_monkey_kali, fval_monkey_kali] = fitDetectionSurface(colors, alphas, initparams_kali);
        
        % particle swarm on the monkey
     %   [fpar_monkey_pso, fval_monkey_pso] = fitDetectionSurface(colors, alphas, 'pso', npsoiters);
        fpar_monkey_pso = fpar_monkey_kali'; % hack
        fval_monkey_pso = 100; % hack
        % find the best method
        method = {'default'; 'kali'; 'pso'};
        fvals = [fval_monkey_default; fval_monkey_kali; fval_monkey_pso]
        fpars = [fpar_monkey_default'; fpar_monkey_kali'; fpar_monkey_pso];
        [~, idx] = min(fvals);
        bestMethod = method{idx}
        fpar_monkey = fpars(idx,:);
        
        % fit the cone noise model's data
        fpar_cones = fitDetectionSurface(modIsoSurf.colorDirs, modIsoSurf.thresholds, 'ellipsoid');
        
        if (sfidx == 4)
            viewSetting = [172.2 45];
        else
            viewSetting = [172.2 10.161];
        end
    
        
        % plot the behavioral data
        for i = 1:2
            figure; hold on;
            title(sprintf('Spatial Frequency = %.2f cpd', sfs(sfidx)))
            coordinates = bsxfun(@times, colors, alphas);
            set(gca, 'linewidth', 1, 'fontsize', 12)
            
            % Plotting the fit
            plotlims = max(abs(coordinates));
            [x,y,z] = meshgrid(linspace(-plotlims(1),plotlims(1),150),linspace(-plotlims(2),plotlims(2),150),linspace(-plotlims(3),plotlims(3),50));
            tmp = [x(:) y(:) z(:)];
            v = sum(abs(tmp * reshape(fpar_monkey(2:end),3,3)).^fpar_monkey(1),2);
            fv = isosurface(x,y,z,reshape(v,size(x,1),size(x,2),size(x,3)),1);
            p = patch(fv);
            set(p,'EdgeColor', 'none', 'FaceAlpha',.8,'FaceColor',[0 .5 0],'Edgealpha',0,'SpecularExponent',.6);
            set(p,'edgecolor','none','facecolor',[0 .6 0]);
            lighting gouraud;
            axis equal;
            xlabel('L'); ylabel('M'); zlabel('S');
            if (i == 1)
                set(gca,'View',viewSetting);
            else
                set(gca,'View',viewSetting+[-90 0]);
            end
            h_light = camlight(20,20);
            plot3(coordinates(:,1),coordinates(:,2),coordinates(:,3),'ko','MarkerFaceColor',[0 .4 0],'MarkerEdgeColor','none','markersize', 4);
            plot3(-coordinates(:,1),-coordinates(:,2),-coordinates(:,3),'ko','MarkerFaceColor',[0 .4 0],'MarkerEdgeColor','none','markersize', 4);
            set(gca,'Zlim',[-.12 .12])
            set(gcf,'renderer','painters'); eval(['print -depsc junk',observer,num2str(sfidx),num2str(i)]);
            set(gca,'Visible','off');% for surface and points
            set(gcf,'Color',[1 1 1]);

            eval(['export_fig ',['junk',observer,num2str(sfidx),num2str(i)],' -tiff -m3;']);
        end
        % make a map of the ratio b/w Monkey and Retina. Where the ratio is large,
        % the monkey does worse than the retina.
        
        [X,Y,Z] = sphere(400);
        spherecolors = [X(:), Y(:), Z(:)];
        spherecolors = bsxfun(@rdivide, spherecolors, sqrt(sum(spherecolors.^2, 2)));
        monkeyThresh = coleThresh(reshape(fpar_monkey(2:end),3,3), fpar_monkey(1), spherecolors);
        conesThresh = coleThresh(reshape(fpar_cones(2:end),3,3), fpar_cones(1), spherecolors);
        C = monkeyThresh ./ conesThresh;
        C = log10(C);
        C = reshape(C, size(X,1), size(X,2)); % for plotting
        minVal = min(C(:))
        minVec = [X(C(:) == minVal), Y(C(:) == minVal), Z(C(:) == minVal)] .* 1.65;
        if size(minVec,1) ==1; minVec = [minVec; -minVec]; end
        maxVal = max(C(:))
        maxVec = [X(C(:) == maxVal), Y(C(:) == maxVal), Z(C(:) == maxVal)] .* 1.65;
        if size(maxVec,1) ==1; maxVec = [maxVec; -maxVec]; end
        
        figure; hold on;
        hand = surf(X,Y,Z,C);
        set(hand, 'edgealpha', 0,'edgecolor','none');
        plot3([-1.4, 1.4], [-1.4 1.4], [0,0], '-', 'linewidth', 5, 'color', 'k') % the L+M direction
        % plot3([-1, 1], [-1 1], [-1 1], '-', 'linewidth', 5, 'color', 'k')           % the achromatic direction
        plot3([0, 0], [0 0], [-1.4 1.4], '-', 'linewidth', 5, 'color', 'b')         % S-iso color direction
        plot3([1.2, -1.2], [-1.2 1.2], [0 0], '-', 'linewidth', 5, 'color', 'r')            % L-M color direction
        %plot3(minVec(:,1), minVec(:,2), minVec(:,3), '--', 'linewidth', 5, 'color', 'm')
        plot3(maxVec(:,1), maxVec(:,2), maxVec(:,3), '--', 'linewidth', 5, 'color', [.3 .3 .3])
        xlim([-1.3, 1.3]); ylim([-1.3, 1.3]); zlim([-1.3, 1.3]);
        set(gca, 'view', viewSetting,...
            'visible', 'off')
        axis square
        maps = {'Edge', 'CubicL', 'CubicYF', 'IsoL'};
        title(sprintf('Spatial Frequency = %.2f cpd', sfs(sfidx)))
        colormap(pmkmp(250, maps{2}))
        colorbar('Ytick',[minVal maxVal],'Yticklabel',num2str([10^minVal 10^maxVal]));
        set(gcf,'renderer','zbuffer'); set(gca,'Visible','off'); eval(['print -dtiff heatmap',observer,num2str(sfidx)])
        drawnow;
        
        % Saving ideal observer and monkey thresholds for the 5 special
        % color directions
        
        if (FIGURE_5)
            close all;
            testcolors = [0 0 0.6364;
                .1273 .1273 .1273;
                0.0636   -0.0636    0.4500;
                 0.0636   -0.0636    -0.4500;
                .09 -.09 0;
                ];
            testcolors = testcolors./repmat(sqrt(sum(testcolors.^2,2)),1,3)
            monkeyThresh = coleThresh(reshape(fpar_monkey(2:end),3,3), fpar_monkey(1), testcolors);
            conesThresh = coleThresh(reshape(fpar_cones(2:end),3,3), fpar_cones(1), testcolors);
            fig5data = [fig5data; observeridx sfidx monkeyThresh' conesThresh'];
        end
        
    end
end

if (FIGURE_5)
    plotcolors = [0 0 1; 0 0 0; 1 0 1; 1 .5 0; 1 0 0];
    figure; axes; hold on;
    h = plot([.5 1 2 4], fig5data(1:4,3:7),'o-','LineWidth',2);
    set(gca,'Yscale','log','XScale','log')
    for i = 1:length(h)
        set(h(i),'Color',plotcolors(i,:));
        set(h(i),'MarkerFaceColor',plotcolors(i,:));
        set(h(i),'MarkerEdgeColor',plotcolors(i,:));
    end
    h = plot([.5 .5 .5 .5], fig5data(1:4,8:12),'o'); % Model preds
    for i = 1:length(h)
        set(h(i),'MarkerEdgeColor',plotcolors(i,:));
    end
    set(gca,'YLim',[.002 .5]);
end


%%
% Section 14)
% Seeing how SMJ cone-isolating stimuli change when we present them to the
% model
% load('/Users/greghorwitz/Documents/Manuscripts/Charlie''s model/For Greg/kali_DTNT_022815.mat');
% load('/Users/greghorwitz/Documents/Manuscripts/Charlie''s model/For Greg/fundamentals.mat');
load ('den_lens_ss');
lenstransmittance = 1./(10.^(den_lens_ss));
lenstransmittance = 1./10.^(den_lens_ss./den_lens_ss(5));
synthfundamentals = fundamentals.*repmat(lenstransmittance',3,1);
synthfundamentals = synthfundamentals./repmat(max(synthfundamentals')',1,81);

lmsCC = [1 -1 0];
%lmsCC = [.14 -.14 .98];

spd = reshape(dtnt.monSpect,81,3);
bkgndrgb = dtnt.bkgndrgb;
load T_cones_smj10
%load T_cones_synthgh2
M1 = T_cones_smj10*spd;
bkgndlms = M1*bkgndrgb';
rgb = inv(M1)*((1+lmsCC').*bkgndlms);

% Going back the other way
(M1*rgb-bkgndlms)./bkgndlms


% New funds
M2 = synthfundamentals*spd;
bkgndlms = M2*bkgndrgb';
a = (M2*rgb-bkgndlms)./bkgndlms
a./norm(a)

%%
% Getting deviations from the intended cone contrasts to those used in the
% model. Now just using the .mat files that Charlie sent me since I feel
% pretty confident now that the numbers are reasonably accurate.
load('/Users/greghorwitz/Library/Containers/com.apple.mail/Data/Library/Mail Downloads/0E348B5F-165F-4A31-A9F9-29BF25467435/smj2.mat')
load('/Users/greghorwitz/Library/Containers/com.apple.mail/Data/Library/Mail Downloads/601B1EC4-49A9-40B1-A02D-FF810C28FC6F/smj10.mat')
data = [smj2; smj10];
colordirs = [0 0 1;
    0.7071 -0.7071 0;
    0.1386 -0.1386 -0.9806;
    0.1386 -0.1386 0.9806];
L = zeros(size(data,1),1);
out = [];
for i = 1:size(colordirs,1)
    L1 = softEq(data(:,1), colordirs(i,1),2) &...
        softEq(data(:,2), colordirs(i,2),2) &...
        softEq(data(:,3), colordirs(i,3),2);
    L2 = softEq(data(:,1), -colordirs(i,1),2) &...
        softEq(data(:,2), -colordirs(i,2),2) &...
        softEq(data(:,3), -colordirs(i,3),2);
    if (sum(L1) > sum(L2))
        L = L1;
    else
        L = L2;        
    end
    
    out = [out; sum(L) mean(data(L,[4 5 6]))];  
end

%%
% Section 15)
% Looking at heat maps of threshold ratios using L:M ratios that are not
% 1:1. And the ideal observer weighting function stuff that he sent on
% 8/16/15.
% Need to load data_new_LtoM_1.mat, etc. Let's navigate to the appropriate
% directory.
cd ('/Users/greghorwitz/Documents/Manuscripts/Completed/Charlie''s model/For Greg');
% Here's where you run "wtfxn_spectrum.m", save a postscript file and
% compare it to the previously submitted figure.

% Now looking at the heat maps of the threshold ratios

whichdatafile = 1;
if (whichdatafile == 1)
    load('data_new_LtoM_1.mat');
elseif (whichdatafile == 2)
    load('data_new_LtoM_1pt5.mat');
elseif (whichdatafile == 3)
    load('data_new_MtoL_1pt5.mat'); 
elseif (whichdatafile == 4)
    load('data_old_LtoM_1.mat');
end
ncones = [sum(cones.num_L(:)) sum(cones.num_M(:)) sum(cones.num_S(:))];
threshold_pts = gab.colorDirs .* repmat(cones.alpha_analytic,1,size(gab.colorDirs,2));
figure; axes; hold on;
tmp = [threshold_pts;-threshold_pts];
plot3(tmp(:,1),tmp(:,2),tmp(:,3),'k.');
% getting radii
D = [tmp(:,1) .* tmp(:,1),...
    tmp(:,2) .* tmp(:,2),...
    tmp(:,3) .* tmp(:,3),...
    2*tmp(:,1) .* tmp(:,2),...
    2*tmp(:,1) .* tmp(:,3),...
    2*tmp(:,2) .* tmp(:,3)];
v = (D' * D) \(D' * ones(size(tmp,1),1));
A = [v(1) v(4) v(5);...
    v(4) v(2) v(6);...
    v(5) v(6) v(3)];
[evecs, evals] = eig(A);
evals = diag(evals);
if ncones(2) > ncones(1)
    evals = evals([2 3 1]);  % back in LMS order. This one works with "data_new_MtoL_1pt5" - why?
else
    evals = evals([3 2 1]);  % back in LMS order. Fragile
end
radii = sqrt(1./evals);
[lcc,mcc,scc] = meshgrid(linspace(-radii(1),radii(1),40),...
    linspace(-radii(2),radii(2),40),...
    linspace(-radii(3),radii(3),40));
fn = (lcc./radii(1)).^2+(mcc./radii(2)).^2+(scc./radii(3)).^2;
p = patch(isosurface(lcc,mcc,scc,fn,1));
set(p,'edgecolor','none','facecolor',[.0 .6 .0])
set(p,'SpecularExponent',1)
lighting gouraud;
axis equal;
radii
% Here are the radii for different L:M cone ratios:
% 1:1 .0022 .0023 .0079
% 1.5:1 .0025 .0021 .0079
% 1:1.5 .0020 .0026 .0079

fpar_cones = fitDetectionSurface(gab.colorDirs, cones.alpha_analytic, 'ellipsoid');
% Now go up to line ~909
% load in a "dtnt" file (e.g. "sedna_DTNT_022815.mat")
% set INHERIT_fpar_cones to 1 
% to compute the heat maps showing threshold ratios
% maxvect (1:1)  = 0.7441    0.7441    1.0690

for whichdatafile = 1:4
    if (whichdatafile == 1)
        load('data_new_LtoM_1.mat');
    elseif (whichdatafile == 2)
        load('data_new_LtoM_1pt5.mat');
    elseif (whichdatafile == 3)
        load('data_new_MtoL_1pt5.mat');
    elseif (whichdatafile == 4)
        load('data_old_LtoM_1.mat');
    end
    fpar_cones_all(:,whichdatafile) = fitDetectionSurface(gab.colorDirs, cones.alpha_analytic, 'ellipsoid');
end
%fpar_cones = fpar_cones_all(:,2)
% Changing from L:M = 1:1 to L:M = 1:1.5 changes the ideal observers'
% threshold by as much as ~12% in the L and M-cone isolating directions.
% Number of L- and M-cones stays constant so a 1:1.5 ratio means
% # L-cones goes down by a factor of .8 and # M-cones goes up by a factor
% of 1.2. 
% Adding 20% more sensors. signal goes up  by 1.2, noise goes up by
% sqrt(1.2). 1.22/sqrt(1.2) = 1.11. Check.

%%
% Section 16)
% Getting some psychophysical/monitor calibration data for Xiaomao Ding,
% David Brainard's undergraduate programmer. He is going to put our
% intermediate color stimuli (L-M-S and L-M+S) through his chromatic
% aberration model to see if the difference is detection thresholds can be
% explained by chromatic aberration.

EQUATECONECONTRASTS = 0; % Make the stimuli identical in cone contrast (using the 
% lower of the two cone contrasts) 

filename = 'S030909003';
start_path = '/Volumes/NO BACKUP/NexFiles/Charlie/Sedna/2009';
stro = nex2stro(findfile(filename,start_path));
spd = reshape(stro.sum.exptParams.mon_spect,81,3);
[thresholds, colordirs, sfs] = DTquestUnpack(stro, 'mode');
M = reshape(stro.sum.exptParams.m_mtx,3,3);
load Dell4BitsCal;
for i = 1:length(cals)
    cal = cals{i};
    calspd = cal.P_device;
    % size(calspd) % always 101x3 or 201x3 
    [calspd] = SplineSpd([380:4:780]', calspd, [380:5:780]');
    if all(calspd == spd)
        break;
    end
    %figure; axes; hold on; 
    %plot(spd);
    %plot(calspd)
end
% "cal" is now the correct calibration structure.
gammaTable = reshape(stro.sum.exptParams.gamma_table,256,3);
RGB = stro.trial(:,[14:16]); % Confirmed these are DAC values
pixperdeg = stro.sum.exptParams.pixperdeg;
sf = pixperdeg./stro.trial(:,21);
colordir = stro.trial(:,26);
bkgndrgb = [stro.sum.exptParams.bkgnd_r;stro.sum.exptParams.bkgnd_g;stro.sum.exptParams.bkgnd_b];
bkgndlms = M*bkgndrgb;
rgb = [interp1(linspace(0,2^16-1,256),gammaTable(:,1),RGB(:,1),'linear'),...
    interp1(linspace(0,2^16-1,256),gammaTable(:,2),RGB(:,2),'linear'),...
    interp1(linspace(0,2^16-1,256),gammaTable(:,3),RGB(:,3),'linear')];
im = {};
for cdir = 1:2
    L = colordir == cdir & sf > 3;
    %plot(RGB(L,:));
    %plot3(cc(1,:),cc(2,:),cc(3,:),'.'); % Checks out
    lasttrialidx = find(L,1,'last');
    stimsizeinsd = stro.sum.exptParams.flash_size;
    theta = stro.trial(lasttrialidx,22);
    lambda = stro.trial(lasttrialidx,21)/10;  % lambda param is wavelength in 10ths of deg
    sigma = stro.trial(lasttrialidx,20)/10; % sigma param is sd of gabor in 10ths of deg (go out 2*sigma)
    gamma = stro.trial(lasttrialidx,23);
    phi = 0;
    if (EQUATECONECONTRASTS & cdir == 2)
        lasttrialidx = find(colordir == 1 & sf > 3,1,'last'); % getting the cc on the last trial of the *other* cdir
        lms = M*rgb(lasttrialidx,:)'; % Getting the cone contrast for the cdir = 1 direction
        cc = (lms(:,end)-repmat(bkgndlms,1,size(lms(:,end),2)))./bkgndlms;
        cc(3) = -cc(3); % flipping the S-cone contrast
        lasttrialidx = find(colordir == 2 & sf > 3,1,'last');
        rgb(lasttrialidx,:) =inv(M)*(bkgndlms.*(cc+1))
    end
    im{cdir} = DrawGaborEdge(bkgndrgb, -bkgndrgb+rgb(lasttrialidx,:)', [0 0 0], theta, lambda, sigma, gamma, phi, 0, 0, 0, 0, [normcdf(stimsizeinsd,0,1) 1-normcdf(stimsizeinsd,0,1)], pixperdeg);
    
    figure; axes; hold on;
    image(im{cdir});
    % Converting to lms
    lms = M*rgb(lasttrialidx,:)';
    cc = (lms-repmat(bkgndlms,1,size(lms,2)))./repmat(bkgndlms,1,size(lms,2));
 
    title(num2str(cc(:,end)));
    axis image;
    
    % Converting back to RGBs (had to convert to rgbs earlier to make the image)
    invGamma = InvertGammaTable(cal.gammaInput, cal.gammaTable, 2^16);
    imgSize = size(im{cdir});
    flashRGB = ones(imgSize); %preallocate so that indexing works
    for gun = 1:3
        tmp = round(im{cdir}(:,:,gun)*2^16); % gun intensities
        flashRGB(:,:,gun) = round(reshape(invGamma(tmp-1, gun), imgSize(1), imgSize(2))*(2^16-1));
    end        
    squeeze(flashRGB(62,62,:)) % This should be the peak of the Gabor
    im{cdir} = flashRGB;
end
    
    
   
   % Extracting thresholds in cone contrasts units from rgbs and the output of
   % DTquestUnpack. Making sure they agree. And they do.
   %thresholds(cdir,1)/100*mkbasis(colordirs(cdir,:)')
   %cc(:,end)

% Trying to figure out what the bkgndRGB is
bkgndrgb_alt = [gammaTable(round(255*cal.bgColor(1))+1,1),...
    gammaTable(round(255*cal.bgColor(2))+1,2),...
    gammaTable(round(255*cal.bgColor(3))+1,3)];
% Good, this is it.

save('ForXiaoMao','cal', 'im');

