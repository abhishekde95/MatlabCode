% Author - Abhishek De, 11/17
% This a spinoff of the WNsimulate based on Greg's suggestion. Now I want
% to create a cell having 2 filters and see what does the STA and STC
% analysis reveal .
close all; clearvars;
plot_counter = 1;
tottrials = 300;
triallenghts = randi([500 1000],tottrials,1);
LMS1 = [0.25 -0.25 0]; % chrom opponent
LMS2 = [0 0.25 0.25]; % chrom non-opponent
filter = zeros(10,10,3); % 3 is for L,M,S channels
filter1 = filter; 
filter2 = filter; 
filter1(3:8,2:5,1) = LMS1(1); filter1(3:8,2:5,2) = LMS1(2); filter1(3:8,2:5,3) = LMS1(3);
filter1(3:8,6:9,1) = -LMS1(1); filter1(3:8,6:9,2) = -LMS1(2); filter1(3:8,6:9,3) = -LMS1(3);
filter2(3:8,2:5,1) = LMS2(1); filter2(3:8,2:5,2) = LMS2(2); filter2(3:8,2:5,3) = LMS2(3);
filter2(3:8,6:9,1) = -LMS2(1); filter2(3:8,6:9,2) = -LMS2(2); filter2(3:8,6:9,3) = -LMS2(3);
RefreshRate = 75;  % Stim refresh rate (Hz)
dtbin = .1; % binsize for Poisson spike generation

for mode = 1:4
    respmat = [];
    stim = [];
    for ii = 1:tottrials
        L = triallenghts(ii);
        input = randn(300,L);
        r1 = sum(repmat(filter1(:),1,L).*input,1);
        r2 = sum(repmat(filter2(:),1,L).*input,1);
        r = getcumdrive2(r1,r2,mode); % cumulative drive
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
    projSTSontofilt1 = filter1(:)'*stim(:,spikeind);
    projSTSontofilt2 = filter2(:)'*stim(:,spikeind);
    projallontoSTA = STA(:)'* stim;
    projallontoPC1 = PC1(:)'* stim;
    projallontoPC2 = PC2(:)'* stim;
    projallontoPCend = PCend(:)'* stim;
    projallontofilt1 = filter1(:)'* stim;
    projallontofilt2 = filter2(:)'* stim;
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
    figure(plot_counter),subplot(3,3,1),image((0.5*filter1./(max(abs(filter1(:)))+0.01)) + 0.5); set(gca,'XTick',[],'YTick',[]); axis square; title('filter1');
    subplot(3,3,2),image((0.5*filter2./(max(abs(filter2(:)))+0.01)) + 0.5); set(gca,'XTick',[],'YTick',[]); axis square; title('filter2');
    subplot(3,3,3),image((0.5*STA./(max(abs(STA(:)))+0.01)) + 0.5); set(gca,'XTick',[],'YTick',[]); axis square; title('STA');
    subplot(3,3,4),image((0.5*PC1./(max(abs(PC1(:)))+0.01)) + 0.5); set(gca,'XTick',[],'YTick',[]); axis square; title('PC1');
    subplot(3,3,5),image((0.5*PC2./(max(abs(PC2(:)))+0.01)) + 0.5); set(gca,'XTick',[],'YTick',[]); axis square; title('PC2');
    subplot(3,3,7),imagesc([xmin xmax],[xmin xmax],non_linSTAPC1); xlabel('projPC1'), ylabel('projSTA'); title('FR map'); axis square;
    subplot(3,3,8),imagesc([xmin xmax],[xmin xmax],non_linPC1PC2); xlabel('projPC1'), ylabel('projPC2'); title('FR map'); axis square;
    subplot(3,3,9),plot(flipud(eig_PC),'*'), ylabel('Eigenvalues'),xlabel('index'),title('Principle components');
    
    % Need to check the FR map when the STS is projected onto individual subunits
    projSTfilters = [projSTSontofilt1' projSTSontofilt2'];
    projallfilters = [projallontofilt1' projallontofilt2'];
    num_bins = 30;
    min_val = floor(min(projallfilters(:))*100)/100;
    max_val = ceil(max(projallfilters(:))*100)/100;
    bin_interval = (max_val-min_val)/num_bins;
    nbins1 = min_val:bin_interval:max_val;
    nbins = linspace(prctile(projallfilters(:),5), prctile(projallfilters(:),95),numel(nbins1)+2);
    new_nbins = [nbins(1)-mean(diff(nbins)), nbins, mean(diff(nbins)) + nbins(end)];
    non_linfilt12 = getFRmap(projSTfilters, projallfilters ,new_nbins);
    xmin = min(nbins); xmax = max(nbins);
    figure(plot_counter),subplot(336),imagesc([xmin xmax],[xmin xmax],non_linfilt12); xlabel('projfilt1'), ylabel('projfilt2'); title('FR map'); axis square;
    plot_counter = plot_counter + 1;
end

