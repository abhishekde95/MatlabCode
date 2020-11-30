% A main script for checking how to regularize the STA
% Author - Abhishek De, 10/2020

close all; clearvars;
plot_counter = 1;

% 1--> Let's do this exercise in 1-D

% There parameters could be tunable
maxT = 15;
gamma = 0.3;
sigma = 0.5;
temporal_profile = exp(-1*(1:1:maxT)*gamma).*cos((0:1:maxT-1)*sigma);

% Creating a WN stim
T = 1000; % in sec
sampling_rate = 0.001; % 1 ms
stim = randn(T/sampling_rate,1);

% Projecting the stim into the temporal kernel and simulating a poisson train
slope = 1.0;
drive = max(0,slope*conv(stim,temporal_profile,'same'));
output = poissrnd(floor(drive));

figure(plot_counter); 
subplot(221); plot(1:maxT, temporal_profile,'Linewidth',2); axis square; set(gca,'Tickdir','out'); 
xlabel('time'); ylabel('Temporal kernel');
subplot(222); plot(drive,'k'); axis square; set(gca,'Tickdir','out'); xlabel('time'); ylabel('drive');
subplot(223); plot(output,'k'); axis square; set(gca,'Tickdir','out'); xlabel('time'); ylabel('drive');

% Now checking whether the temporal kernel can be extracted using STCOVmex 
STCOVmex('init', {1 1 maxT});
STCOVmex(stim,output);
out = STCOVmex('return'); 
STS = out{1};  
nspikes = out{3}; 
clear STCOVmex; clear out;
recovered_filter = STS/nspikes;

figure(plot_counter);
subplot(224); plot(-1*recovered_filter,'g','Linewidth',2); axis square; set(gca,'Tickdir','out'); xlabel('time'); ylabel('Magnitude')
title('Recovered filter');
plot_counter = plot_counter + 1;

% STA is able to recover the temporal filter

%% 2-->Trying a different method of projecting a stimulus onto the recovered filter

projs = conv(stim,-1*recovered_filter,'same');

figure(plot_counter);
subplot(121); plot(projs); axis square; set(gca,'Tickdir','out'); 
xlabel('Time'); ylabel('Signal');
subplot(122); histogram(projs(output==0)); hold on; histogram(projs(output>0)); 
axis square; xlabel('Proj'); set(gca,'Tickdir','out'); 
plot_counter = plot_counter + 1;

roc_val = max([1-rocN(projs(output==0),projs(output>0)) rocN(projs(output==0),projs(output>0))]);
%% 3--> Next task is to extract a regularized version of the recovered filter

% Testing raised cosine vectors 
nkt = 15; % number of ms in stim filter
neye = 0; % number of "identity" basis vectors near time of spike;
ncos = 4; % number of raised-cosine vectors to use
kpeaks = [.1 round(nkt/1.2)];  % position of first and last bump (relative to identity bumps)
b = 10;
kdt = 2;  % spacing of x axis must be in units of 1

nlin = @(x)log(x+1e-20);
invnl = @(x)exp(x)-1e-20; % inverse nonlinearity

yrnge = nlin(kpeaks+b);  
db = diff(yrnge)/(ncos-1);  % spacing between raised cosine peaks
ctrs = yrnge(1):db:yrnge(2);  % centers for basis vectors
mxt = invnl(yrnge(2)+2*db)-b; % maximum time bin
kt0 = [0:kdt:mxt]';
nt = length(kt0);        % number of points in iht
ff = @(x,c,dc)(cos(max(-pi,min(pi,(x-c)*pi/dc/2)))+1)/2; % raised cosine basis vector
kbasis0 = ff(repmat(nlin(kt0+b), 1, ncos), repmat(ctrs, nt, 1), db);

% Concatenate identity-vectors
nkt0 = size(kt0,1);
kbasis = [[eye(neye); zeros(nkt0,neye)] [zeros(neye, ncos); kbasis0]];
% kbasis = flipud(kbasis);  % flip so fine timescales are at the end.
nkt0 = size(kbasis,1);
kbasis = normalizecols(kbasis);

figure(plot_counter);
plot(kbasis, 'Linewidth', 2), axis square; 
set(gca,'Tickdir','out');
plot_counter = plot_counter + 1;

% Fitting the function
filter_reg = regularize_filter(stim, output, -1*recovered_filter, [], 1); 
filter_reg_L1 = regularize_filter(stim, output, -1*recovered_filter, [], 2); 
filter_reg_L2 = regularize_filter(stim, output, -1*recovered_filter, [], 3);
filter_reg_GH = regularize_filter(stim, output, -1*recovered_filter, kbasis, 4);

figure(plot_counter);
subplot(311); plot(temporal_profile,'k','Linewidth',2); axis square;
set(gca,'Tickdir','out'); xlabel('Time'); ylabel('Original filter');

subplot(312); plot(-1*recovered_filter,'g','Linewidth',2); axis square;
set(gca,'Tickdir','out'); xlabel('Time'); ylabel('Recovered filter-STA');

subplot(313); plot(filter_reg,'r','Linewidth',2); hold on; 
plot(filter_reg_L1, 'b', 'Linewidth',2);
plot(filter_reg_L2, 'm', 'Linewidth',2);
plot(filter_reg_GH, 'k','Linewidth',2);
axis square; set(gca,'Tickdir','out'); xlabel('Time'); title('Regularized filter')
plot_counter = plot_counter + 1;

% Some stats on the recovered and regularized filter
[r1,p1] = corr(-1*recovered_filter',temporal_profile');
[r2,p2] = corr(filter_reg,temporal_profile');
[r3,p3] = corr(filter_reg_L1,temporal_profile');
[r4,p4] = corr(filter_reg_L2,temporal_profile');
[r5,p5] = corr(filter_reg_GH,temporal_profile');

% Verdict --> Filter with cosine bumps is working and gives the best/smoothed response.

