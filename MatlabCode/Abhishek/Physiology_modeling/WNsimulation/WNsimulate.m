% I am writing this script so as to simulate some of my unanswered
% questions: 1) Does the non-linear spatial integration by a filter
% manifests as multiple filters? 
% 2) How do the nonlinearities affect what filters are recovered?
% I will start with these two questions and see where I reach. This is
% gonna be a massive project, might take me days to implement.
% Author - Abhishek De, 11/17

% First generate a battery of stimulus, then get the response from a
% filter. Implement this thing first and then think about how to estimate
% STA and STC from the spiking responses. You might have to take help from
% some of the functions from code_stc_JP - This part is done

% Generating stimulus frames 
% Each frame is 10x10x3 dimensions
close all; clearvars;
plot_counter = 1;
tottrials = 120;
triallenghts = randi([500 1000],tottrials,1);
LMS1 = randn(1,3);
LMS2 = randn(1,3);
% LMS1 = [-0.25 0.25 0];
% LMS2 = -LMS1;
filter = zeros(10,10,3); % 3 is for L,M,S channels
subunit1 = filter;
subunit2 = filter;
subunit1(3:8,2:5,1) = LMS1(1); subunit1(3:8,2:5,2) = LMS1(2); subunit1(3:8,2:5,3) = LMS1(3);
subunit2(3:8,6:9,1) = LMS2(1); subunit2(3:8,6:9,2) = LMS2(2); subunit2(3:8,6:9,3) = LMS2(3);
filter = subunit1 + subunit2;
RefreshRate = 75;  % Stim refresh rate (Hz)
dtbin = .1; % binsize for Poisson spike generation

for mode = [1 6]
    respmat = [];
    stim = [];
    for ii = 1:tottrials
        L = triallenghts(ii);
        input = randn(300,L);
        r1 = sum(repmat(subunit1(:),1,L).*input,1);
        r2 = sum(repmat(subunit2(:),1,L).*input,1);
        r = (max(0,getcumdrive(r1,r2,mode))).^2; % spiking rectified quadratic nonlinearity
        rbig = repmat(r/RefreshRate*dtbin,1./dtbin,1); % make Poisson spike train
        sp = sum(rand(size(rbig))<rbig)';
        respmat = [respmat sp'];
        stim = [stim input];
    end
    
    % calculating STS and STA
    spikeind = find(respmat>0);
    nspikes = sum(respmat);
    STS = sum(stim(:,spikeind).*repmat(respmat(spikeind),[300,1]),2); % spike triggered ensemble
    STA = STS./nspikes;
    STA = reshape(STA,[10 10 3]);
    STCross = zeros(300,300); % calculating STCross
    for ii = 1:numel(spikeind)
        STCross = STCross + (respmat(spikeind(ii)).^2) * (stim(:,spikeind(ii)) * stim(:,spikeind(ii))');
    end
    
    % calculating STC so that I can extract the principal components - need to calculate STCross
    tmp = STS(:)*STS(:)';
    STCs = (nspikes.*STCross-tmp)/(nspikes*(nspikes-1));
    P = eye(size(STCs)) - STA(:)*inv(STA(:)'*STA(:))*STA(:)'; % WHAT DOES THIS LINE MEAN
    STCs = P*STCs*P';
    [tmp,d] = eig(STCs);
    eig_PC = sort(diag(d)); % storing all the eigenvalues
    v = real(tmp);
    [~, idxs] = sort(diag(d));
    v = v(:,idxs);
    suppresive_PC = 2;
    PC1 = v(:,end); % first principal component
    PC2 = v(:,end-1); % second principal component
    PCend = v(:,2); % last principal component
    PC1 = reshape(PC1,[10 10 3]);
    PC2 = reshape(PC2,[10 10 3]);
    PCend = reshape(PCend,[10 10 3]);
    
    % Projecting onto STA and PC1
    projSTSontoSTA = STA(:)'* stim(:,spikeind);
    projSTSontoPC1 = PC1(:)'* stim(:,spikeind);
    projSTSontoPC2 = PC2(:)'* stim(:,spikeind);
    projSTSontoPCend = PCend(:)'* stim(:,spikeind);
    projSTSontosub1 = subunit1(:)'*stim(:,spikeind);
    projSTSontosub2 = subunit2(:)'*stim(:,spikeind);
    projallontoSTA = STA(:)'* stim;
    projallontoPC1 = PC1(:)'* stim;
    projallontoPC2 = PC2(:)'* stim;
    projallontoPCend = PCend(:)'* stim;
    projallontosub1 = subunit1(:)'* stim;
    projallontosub2 = subunit2(:)'* stim;
    projST = [projSTSontoSTA' projSTSontoPC1' projSTSontoPC2' projSTSontoPCend'];
    projall = [projallontoSTA' projallontoPC1' projallontoPC2' projallontoPCend'];
    min_val = floor(min(projall(:))*100)/100;
    max_val = ceil(max(projall(:))*100)/100;
    num_bins = 15;
    bin_interval = (max_val-min_val)/num_bins;
    nbins1 = min_val:bin_interval:max_val;
    nbins = linspace(prctile(projall(:),5), prctile(projall(:),95),numel(nbins1)+2);
    new_nbins = [nbins(1)-mean(diff(nbins)), nbins, mean(diff(nbins)) + nbins(end)];
    non_linSTAPC1 = getFRmap([projST(:,1) projST(:,2)],[projall(:,1) projall(:,2)] ,new_nbins);
    non_linPC1PC2 = getFRmap([projST(:,3) projST(:,2)],[projall(:,3) projall(:,2)] ,new_nbins);
    non_linPC1PCend = getFRmap([projST(:,4) projST(:,2)],[projall(:,4) projall(:,2)] ,new_nbins);
    xmin = min(nbins); xmax = max(nbins);
    
    % plotting the figures
    figure(plot_counter),subplot(2,5,1),image((0.5*filter./(max(abs(filter(:)))+0.01)) + 0.5); set(gca,'XTick',[],'YTick',[]); axis square; title('Actual filter');
    subplot(2,5,2),image((0.5*STA./(max(abs(STA(:)))+0.01)) + 0.5); set(gca,'XTick',[],'YTick',[]); axis square; title('STA');
    subplot(2,5,3),image((0.5*PC1./(max(abs(PC1(:)))+0.01)) + 0.5); set(gca,'XTick',[],'YTick',[]); axis square; title('PC1');
    subplot(2,5,4),image((0.5*PC2./(max(abs(PC2(:)))+0.01)) + 0.5); set(gca,'XTick',[],'YTick',[]); axis square; title('PC2');
    subplot(2,5,5),image((0.5*PCend./(max(abs(PCend(:)))+0.01)) + 0.5); set(gca,'XTick',[],'YTick',[]); axis square; title('PCend');
    subplot(2,5,7),imagesc([xmin xmax],[xmin xmax],non_linSTAPC1); xlabel('projPC1'), ylabel('projSTA'); title('FR map'); axis square;
    subplot(2,5,8),imagesc([xmin xmax],[xmin xmax],non_linPC1PC2); xlabel('projPC1'), ylabel('projPC2'); title('FR map'); axis square;
    subplot(2,5,9),imagesc([xmin xmax],[xmin xmax],non_linPC1PCend); xlabel('projPC1'), ylabel('projPCend'); title('FR map'); axis square;
    
    % Need to check the FR map when the STS is projected onto individual subunits
    projSTsubunits = [projSTSontosub1' projSTSontosub2'];
    projallsubunits = [projallontosub1' projallontosub2'];
    num_bins = 30;
    min_val = floor(min(projallsubunits(:))*100)/100;
    max_val = ceil(max(projallsubunits(:))*100)/100;
    bin_interval = (max_val-min_val)/num_bins;
    nbins1 = min_val:bin_interval:max_val;
    nbins = linspace(prctile(projallsubunits(:),5), prctile(projallsubunits(:),95),numel(nbins1)+2);
    new_nbins = [nbins(1)-mean(diff(nbins)), nbins, mean(diff(nbins)) + nbins(end)];
    non_linsub1sub2 = getFRmap(projSTsubunits, projallsubunits ,new_nbins);
    xmin = min(nbins); xmax = max(nbins);
    figure(plot_counter),subplot(2,5,6),imagesc([xmin xmax],[xmin xmax],non_linsub1sub2); xlabel('projsub1'), ylabel('projsub2'); title('FR map'); axis square;
    plot_counter = plot_counter + 1;
end

% The conclusion of this exercise is fruitful. So whenever there is
% nonlinear spatial integration, I do observe multiple filters emerge as
% STA and PCs, which is what is mentioned in the paper by Pillow (JOV) and
% Rust (NEURON). What I don't know from this exercise is how to predict
% the nature of nonlinear spatial integration from multiple filters. It seems that when I project the Raw spike triggered stimulus onto the subunits,
% the nonlinear spatial summation comes up. In our experimental setup, the
% contrast in the WhiteNoise paradigm is low such that the FRmap is not
% much informative. Plus I think, I need to hold the neuron for a very long
% time to get a very good estimate of FRmap.
% Now would be a good time to look at the STAPCandspatialsummation/mainanalysis.m 