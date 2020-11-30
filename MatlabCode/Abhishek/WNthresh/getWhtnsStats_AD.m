function out = getWhtnsStats_AD(stro,nframesback,module, initargs, nrandnums_perchannel, spikename, M, allstim)

% Similar implimentation as Greg's but with capability of dealing with subunits
% out = getWhtnsStats(stro,maxT,module,initargs, <spikename>, <M>, <allstim>)
%
%	 Based on getWhtnsStats3, this function will iteratively call a "module"
% function with vectors of RBGs and spikes (once per trial).  "M" is an
% optional 3x3 matrix that transforms the color space prior to the
% projections being computed. GDLH 1/29/17.
%
%    Initargs should be a cell array containing the initialization
%    arguments appropriate for the module called.
%
% GDLH 3/5/08


gammaTable = stro.sum.exptParams.gamma_table; %get gamma_table
gammaTable = reshape(gammaTable,length(gammaTable)/3,3); %reshapse gamma_table into three columns
invgammaTable = InvertGamma(gammaTable,1); %invert gamma_table
ngammasteps = size(invgammaTable,1); %get number of rows of gamma_table (65536)
seedidx = strcmp(stro.sum.trialFields(1,:),'seed'); %get seed index from trialFields
noisetypeidx = strcmp(stro.sum.trialFields(1,:),'noise_type'); %get noise_type index from trialFields
nframesidx = strcmp(stro.sum.trialFields(1,:),'num_frames'); %get nframes index from trialFields
stimonidx = strcmp(stro.sum.trialFields(1,:),'stim_on'); %get stimon index from trialFields
muidxs = [find(strcmp(stro.sum.trialFields(1,:),'mu1')),... %get mu indices into vector from trialFields
         find(strcmp(stro.sum.trialFields(1,:),'mu2')),...
         find(strcmp(stro.sum.trialFields(1,:),'mu3'))];
sigmaidxs = [find(strcmp(stro.sum.trialFields(1,:),'sigma1')),... %get sigma indices into vector from trialFields
            find(strcmp(stro.sum.trialFields(1,:),'sigma2')),...
            find(strcmp(stro.sum.trialFields(1,:),'sigma3'))];
% bkgndRGBidxs = [find(strcmp(stro.sum.trialFields(1,:),'bkgnd_r')),...
%                find(strcmp(stro.sum.trialFields(1,:),'bkgnd_g')),...
%                find(strcmp(stro.sum.trialFields(1,:),'bkgnd_b'))];
if (nargin < 6)
    spikeidx = find(strcmp(stro.sum.rasterCells(1,:),'sig001a')); %get spike index from rasterCells
else
    spikeidx = find(strcmp(stro.sum.rasterCells(1,:),spikename)); %get spike index from rasterCells
end
if (nargin < 7 | isempty(M))
    bigM = [];
elseif (size(M,1) ~= 3 || size(M,2) ~= 3)
   error('M should be 3x3'); 
else
    bigM = kron(M, eye(nstixperside^2));
end
if (nargin < 8)
   allstim = 0; 
end
msperframe = 1000/stro.sum.exptParams.framerate; %calculate msec per frame
if isfield(stro.sum.exptParams,'nrepframes')
    if isnan(stro.sum.exptParams.nrepframes)
        nvblsperstimupdate = 1;
            else
        nvblsperstimupdate = stro.sum.exptParams.nrepframes;
    end
else
    nvblsperstimupdate = 1;
end
ntrials = size(stro.trial,1); %get number of trials

% compute once for gun noise
xx = linspace(stro.sum.exptParams.gauss_locut/1000, stro.sum.exptParams.gauss_hicut/1000, ...
    ngammasteps); %dividing gauss by 1000 and making equal intervals so that there are 65536 values
yy = norminv(xx');

feval(module,'init',initargs);

hwait = waitbar(0,'Finding stimuli...'); %display progress bar
for i = 1:ntrials %get values and insert into given column
    nframes = stro.trial(i,nframesidx);
    if nframes == 0, continue; end
    seed = stro.trial(i,seedidx);
    noisetype = stro.trial(i,noisetypeidx);
    mu = stro.trial(i,muidxs)/1000;
    sigma = stro.trial(i,sigmaidxs)/1000;
    
    if (noisetype == 1)  % Gaussian gun
        invnormcdf = bsxfun(@plus, bsxfun(@times, yy, sigma), mu);
        randnums = getEJrandnums(3*nrandnums_perchannel*nframes, seed);
        randnums = reshape(randnums, [nrandnums_perchannel*3, nframes]);
        for j = 1:3
            idxs = (1:nrandnums_perchannel)+nrandnums_perchannel*(j-1);
            randnums(idxs,:) = reshape(invnormcdf(randnums(idxs,:)+1,j),[length(idxs),nframes]); % order in drawing on each frame frame is: red (all pixels), green (all pixels), blue (all pixels).
        end
        if (~isempty(bigM))
            randnums = bigM*randnums; % at this point, randnums has 3*n rows (where n is nstixperframe*nstixperframe) and nframes columns
        end
    elseif (noisetype == 2)  % Binary cone
        colordirlms = sign(fullfact([2,2,2])-1.5);
        randnums = getEJrandnums(nrandnums_perchannel*nframes, seed);
        randnums = reshape(randnums, [nrandnums_perchannel, nframes]);
        randnums = mod(randnums, size(colordirlms,1))+1;
        tmp = colordirlms(randnums,:);
        tmp = reshape(tmp, [nrandnums_perchannel nframes 3]);
        tmp = permute(tmp,[1 3 2]); % space, cone, time
        randnums = reshape(tmp,[nrandnums_perchannel*3 nframes]);
        if (~isempty(bigM))
            randnums = bigM*randnums;
        end
    else
        continue;
    end
    t_stimon = stro.trial(i, stimonidx);
    spiketimes = (stro.ras{i,spikeidx}-t_stimon)*1000;  % converting to ms
    frametimes = linspace(0, nframes*msperframe*nvblsperstimupdate, nframes)+(msperframe*nvblsperstimupdate/2)';
    spiketimes(spiketimes < nframesback*msperframe) = [];
    spiketimes(spiketimes > frametimes(end)) = [];
    n = hist(spiketimes, frametimes);
    if allstim
        n = ones(size(n)); % Extracting all the stimuli 
    end
    feval(module,randnums(:),n);
    if (ishandle(hwait))
        waitbar(i/ntrials,hwait);
    else
        break;
    end
end

if (ishandle(hwait))
   close(hwait);
   out = feval(module,'return');
end
eval(['clear ',module]);

