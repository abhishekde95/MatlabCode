% Testing the sensitivity of the Benardete and Kaplan LGN P-cell model to
% my stimuli (assuming Poisson variabilty)
% Assuming a baseline of, what, 40 sp/s?

filename = 'A061517002.nex'; % To get baseline firing rate and contrasts
whichparameters = 'Figure 6'; % 'Figure 6' or 'Table'. 'Table' is parvo cell, figure 6 is example midget RGC
i = sqrt(-1); % In case i is already defined
ntrials = 5000;

if strcmp(whichparameters,'Figure 6')
    % Taking parameters from Figure 6 legend
    % centers
    center.A.mean = 184.20; % spikes/sec/unit contrast
    center.Hs.mean = 0.69;% s
    center.Ts.mean = 18.61/1000; % s
    center.Nl.mean = 38;
    center.Tl.mean = 1.23/1000;  % s
    center.D.mean = 4.0/1000; % s
    
    % surrounds
    surround.A.mean = 125.33; % spikes/sec/unit contrast
    surround.Hs.mean = .56; % s
    surround.Ts.mean = 33.28/1000; % s
    surround.Nl.mean = 124;
    surround.Tl.mean = .42/1000;  % s
    surround.D.mean = 4.0/1000; % s
elseif strcmp(whichparameters,'Table')
    %
    % ---------- From the table (population) ----------
    % centers (from Table 5)
    center.A.mean = 34.93; % spikes/sec/unit contrast
    center.Hs.mean = .72;
    center.Ts.mean = 46.19/1000; % s
    center.Nl.mean = 27.13;
    center.Tl.mean = 52.48/center.Nl.mean/1000;  % parameterized this way because it's time to peak of STA
    center.D.mean = 3.4/1000; % s (Based on Table 5)
    
    % surrounds (from Table 6)
    surround.A.mean = 37.41; % spikes/sec/unit contrast
    surround.A.mean = 17.41; % actually median

    surround.Hs.mean = .26;
    surround.Ts.mean = 57.31/1000; % s
    surround.Nl.mean = 61.88;
    surround.Tl.mean = 59.97/surround.Nl.mean/1000;  % parameterized this way because it's time to peak of STA
    surround.D.mean = 3.5/1000; % s (Based on Table 3 and Table 6)
else
    error('Unknown parameter set');
end

TF = [0:1:300];
n = length(TF);
omega = TF*2*pi; % rad/s
freqresponse_c = center.A.mean*exp(-i.*omega.*center.D.mean).*...
    (1-(center.Hs.mean./(1+i.*omega.*center.Ts.mean))).*...
    (1./(1+i.*omega.*center.Tl.mean)).^center.Nl.mean;

freqresponse_s = surround.A.mean*exp(-i.*omega.*surround.D.mean).*...
    (1-(surround.Hs.mean./(1+i.*omega.*surround.Ts.mean))).*...
    (1./(1+i.*omega.*surround.Tl.mean)).^surround.Nl.mean;
% Below from http://www.mathworks.com/matlabcentral/newsreader/view_thread/5658
% Getting rid of the highest half of the frequency spectrum
if (rem(n,2))
     freqresponse_c((n+3)/2:n)=[]; centerkernel=ifft([freqresponse_c conj(freqresponse_c((n+1)/2:-1:2))]);
     freqresponse_s((n+3)/2:n)=[]; surroundkernel=ifft([freqresponse_s conj(freqresponse_s((n+1)/2:-1:2))]);
else
    freqresponse_c((n+4)/2:n)=[]; centerkernel=ifft([freqresponse_c conj(freqresponse_c(n/2:-1:2))]);
    freqresponse_s((n+4)/2:n)=[]; surroundkernel=ifft([freqresponse_s conj(freqresponse_s(n/2:-1:2))]);
end
% Note: integral of kernel doesn't depend on delta-T (max TF)
x = [0:1/(TF(end)):1/(TF(end))*(length(centerkernel)-1)]; 
% Due to the bit of Mathworks code above, the highest frequency we're considering is max(TF/2)
% That's why sampling is is steps of 1/(TF(end)) instead of 1/(2*TF(end)).

% Getting rid of "white space"
centerkernel(x>0.5) = [];
surroundkernel(x>0.5) = [];
x(x>0.5) = [];

figure;
subplot(2,2,1);
plot(x,centerkernel,'.-'); set(gca,'Xlim',[0 .5]);
subplot(2,2,2);
plot(x,surroundkernel,'.-'); set(gca,'Xlim',[0 .5]);
subplot(2,2,3); hold on;
plot(x, centerkernel+surroundkernel,'r-','LineWidth',2);
plot(x, centerkernel-surroundkernel,'k-','LineWidth',2); set(gca,'Xlim',[0 .5]);
deltaT = x(2)-x(1);

% Loading an IsoSamp file so I can get the contrasts and basseline firing
% rate
stro = nex2stro(findfile(filename));
Lcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_l'));
Mcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_m'));
TF = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'tf'));
stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimon_t'));
stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimoff_t'));
dur = mean(stimoff_t-stimon_t);
spikeidx = strcmp(stro.sum.rasterCells(1,:),getSpikenum(stro));

uniquestim = sortrows(unique([Lcc Mcc TF],'rows'),3); % sorting by TF
uniquestim(all(uniquestim == 0,2),:) = []; % Getting rid of blank
Lblank = Lcc == 0 & Mcc == 0 & TF == 0;

baselinespikerates = [];
for idx = find(Lblank)'
    spiketimes = stro.ras{idx,spikeidx};
    nspikes = sum(spiketimes > stimon_t(idx) & spiketimes < stimon_t(idx)+dur);
    baselinespikerates = [baselinespikerates; nspikes./(stimoff_t(idx)-stimon_t(idx))];
end
baseline = mean(baselinespikerates); % sp/s

nbins = round(.666/deltaT); % time bins for stimulus
temporalenvelope = ones(1,nbins);
temporalenvelope(1:round((1/4)*nbins)) = linspace(0,1,round((1/4)*nbins));
temporalenvelope(end:-1:round((3/4)*nbins)+1) = linspace(0,1,round((1/4)*nbins));
timebins = 0:deltaT:1; % time in s
temporalenvelope = [zeros(1,length(timebins)-length(temporalenvelope)),...
    temporalenvelope];
signal = [];
exampletrial = stro.trial(1,:);
stro.trial = []; stro.ras = [];
for i = 1:size(uniquestim,1)
    tf = uniquestim(i,3);
    centerstim = uniquestim(i,1).*temporalenvelope.*sin(2*pi*tf.*timebins);
    surroundstim = uniquestim(i,2).*temporalenvelope.*sin(2*pi*tf.*timebins);
    centerresp = conv(centerstim,centerkernel,'same');
    surroundresp = conv(surroundstim,surroundkernel,'same');
    figure; subplot(2,1,1); hold on;
    plot(timebins, centerresp)
    plot(timebins, surroundresp)
    plot(timebins, centerresp-surroundresp)
    % Quick and dirty energy calculation
    basis1 = exp(-2*pi*sqrt(-1)*tf*[0:length(timebins)-1]/length(timebins))';
    F1 = 2*abs((centerresp-surroundresp)*basis1);
    signal(i) = F1;
    
    % getting a firing rate for a Monte Carlo Poisson model
    % (Monte Carlo allows me to use IsoSampgetDPrime.m)
    fr = max(baseline+(centerresp-surroundresp),0);
    start_counting_t = .1;
    stop_counting_t = start_counting_t+.666;
    nspikes = zeros(ntrials,size(timebins,2));
    exampletrial(strcmp(stro.sum.trialFields(1,:),'stim_l')) = uniquestim(i,1);
    exampletrial(strcmp(stro.sum.trialFields(1,:),'stim_m')) = uniquestim(i,2);
    exampletrial(strcmp(stro.sum.trialFields(1,:),'tf')) = uniquestim(i,3);
    exampletrial(strcmp(stro.sum.trialFields(1,:),'stimon_t')) = start_counting_t; % Empirical, this should be more principled!
    exampletrial(strcmp(stro.sum.trialFields(1,:),'stimoff_t')) = stop_counting_t;
    for j = 1:ntrials
        nspikes(j,:) = poissrnd(fr.*deltaT);
        stro.trial = [stro.trial; exampletrial];
        spiketimes = [];
        for k = 1:max(nspikes(j,:))
            spiketimes = [spiketimes timebins(nspikes(j,:)>=k)+deltaT.*unifrnd(-.5,.5,1,sum(nspikes(j,:)>=k))];
        end
        stro.ras{size(stro.ras,1)+1,1} = spiketimes;
        stro.ras{size(stro.ras,1),2} = nan;
    end
    subplot(2,1,2);
    imagesc(nspikes)
    %plot(timebins,sum(nspikes));
end
% adding blank trials
exampletrial(strcmp(stro.sum.trialFields(1,:),'stim_l')) = 0;
exampletrial(strcmp(stro.sum.trialFields(1,:),'stim_m')) = 0;
exampletrial(strcmp(stro.sum.trialFields(1,:),'tf')) = 0;
for j = 1:ntrials
    stro.trial = [stro.trial; exampletrial];
    stro.ras{size(stro.ras,1)+1,1} = unifrnd(start_counting_t,stop_counting_t,1,poissrnd(baseline*.666));
    stro.ras{size(stro.ras,1),2} = nan;
end
[uniquestim, dprime, signal, noise] = IsoSampGetDPrime(stro,1);


% ploting signal strength as a function of TF for color and luminance
% separately.
% Bottom line: Signal strength is in good qualitative agreement with
% d-primes from parvocells (red-green > luminance; high frequency > low
% frequency)

figure; axes; hold on;
for stimtype = 1:2 % non-opponent, opponent
    if (stimtype == 1)
        Lstimtype = sign(uniquestim(:,1)) == sign(uniquestim(:,2));
    else
        Lstimtype = sign(uniquestim(:,1)) ~= sign(uniquestim(:,2));
    end
    h = plot(uniquestim(Lstimtype,3),dprime(Lstimtype),'.-');
    if (stimtype == 1)
        set(h,'Color','black');
    else
        set(h,'Color','red');
    end
end
set(gca,'Xscale','log');

% Now I need a noise model
% Is this as simple as signal/sqrt(signal)?

% % Figuring out how to convert Lcc&Mcc into luminance contrast
% % (It's not straight forward and depends on the device primaries)
% % Bottom line: It's approximately true that .6*Lcc+.4Mcc = luminance contrast
% % So I can use Lcc (or Mcc) as a good proxy for luminance contrast
% % If Lcc = Mcc = Scc = X, then luminance contrast = X.

% load ProPixx
% P_device = SplineSpd([380:4:780]', cals{end}.P_device, [380:5:780]');
% load T_cones_smj
% M = T_cones_smj*P_device;
% bkgndrgb = [.5 .5 .5];
% bkgndlms = M*bkgndrgb';
% cc = [1 1 0]';
% stimlms = bkgndlms.*(1+cc);
% stimrgb = inv(M)*stimlms
% spectrum_bkgnd = P_device*bkgndrgb';
% spectrum_stim = P_device*stimrgb;
% 
% load T_xyz1964.mat
% vlambda = T_xyz1964(2,:);
% L_fund = T_cones_smj(1,:);
% M_fund = T_cones_smj(2,:);
% 
% % pure radiance increment
% Lcc = (L_fund*spectrum_bkgnd-L_fund*spectrum_stim)/(L_fund*spectrum_bkgnd)
% Mcc = (M_fund*spectrum_bkgnd-M_fund*spectrum_stim)/(M_fund*spectrum_bkgnd)
% LUMcc = (vlambda*spectrum_bkgnd-vlambda*spectrum_stim)/(vlambda*spectrum_bkgnd)

