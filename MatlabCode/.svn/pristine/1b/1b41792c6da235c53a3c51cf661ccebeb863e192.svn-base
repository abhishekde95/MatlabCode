% Figures for the Annual Reviews of Vision Science paper 2020
% Section 1: Rasters and response surfaces for blue-yellow neurons with
% luminance PC1s

%%
% Section 1

filename = 'K041008001.nex';
%filename = 'K052708001.nex';

if strcmp(filename,'K041008001.nex')
    spikename = 'sig001a'; % OK for both example cells
    image_idx_lims = [418 544];
elseif strcmp(filename,'K052708001.nex')
    spikename = 'sig001a'; % OK for both example cells
    image_idx_lims = [250 419];
end
if ~exist('WN') % For working at home
    WN=nex2stro(findfile(filename));
end
framerate = WN.sum.exptParams.framerate;
nstixperside = WN.sum.exptParams.nstixperside;
ntrials = length(WN.sum.absTrialNum);
stimonidx = find(strcmp(WN.sum.trialFields(1,:),'stim_on'));
stimoffidx = find(strcmp(WN.sum.trialFields(1,:),'all_off'));
nframesidx = find(strcmp(WN.sum.trialFields(1,:),'num_frames'));
noisetypeidx = find(strcmp(WN.sum.trialFields(1,:),'noise_type'));
sigmaidxs = strmatch('sigma',WN.sum.trialFields(1,:));

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

spikeidx = find(strcmp(WN.sum.rasterCells(1,:),spikename));
maxT = 9;
Lconenoise = WN.trial(:,noisetypeidx) == 2;
if any(Lconenoise)
    error
    tmpstro.ras(Lconenoise,:) = [];
    tmpstro.trial(Lconenoise,:) = [];
end
out = getWhtnsStats(WN,maxT,'STCOVmex', {nstixperside^2, 3, maxT}, spikename);
STAs = out{1};
STCs = out{2};
nspikes = out{3};

PCs = [];
for i = 1:size(STCs,2)
    STC = reshape(STCs(:,i), 3*nstixperside^2, 3*nstixperside^2);
    subSTC = STC;
    subSTA = STAs(:,i);
    P = eye(size(subSTC))-subSTA*inv(subSTA'*subSTA)*subSTA';
    subSTC = P*subSTC*P';
    [tmp,d] = eig(subSTC);
    %v = repmat(double(Lmask),[1 size(tmp,2)]);
    %v(Lmask,:) = tmp;
    [newd, idxs] = sort(diag(d));
    v = tmp(:,idxs);
    v = v(:,[end end-1]); 
    PCs = cat(3, PCs, v);
end

% Normalizing images
template = reshape([1:nstixperside^2],nstixperside,nstixperside);
edgepixels = [template(:,1); template(1,[2:end-1])'; template(end,[2:end-1])'; template(:,end)];
edgepixelidxs = [edgepixels; edgepixels+nstixperside^2; edgepixels+2*(nstixperside^2)];
PCelements = PCs(edgepixelidxs,:,:);
PCsds = std(PCelements);    % One std calculated per PC
PCs = PCs.*repmat(std(STAs(:,1))./PCsds,[size(PCs,1),1,1]);

rowidxs = reshape([1:3*nstixperside^2],[nstixperside^2 3]);
maxes = []; mins = [];
imagevectors = [STAs, reshape(PCs,[size(STAs,1) size(PCs,2)*size(PCs,3)])];
for i = 1:3
    maxes = [maxes; max(max(imagevectors(rowidxs(:,i),:)))];
    mins = [mins; min(min(imagevectors(rowidxs(:,i),:)))];
end
potentialnormfactors = [(1-[.5; .5; .5]-eps)./(maxes-[.5; .5; .5]); (-[.5; .5; .5]+eps)./(mins-[.5; .5; .5])];
% 'eps' in above line is a kludge that is required for avoiding
% out of bounds errors.
potentialnormfactors(potentialnormfactors < 0) = []; % if min > mu or max < mu
normfactor = min(potentialnormfactors);
%muvect = reshape(repmat(bkgndrgb',nstixperside^2,1),nstixperside^2*3,1);
muvect = reshape(repmat([.5 .5 .5],nstixperside^2,1),nstixperside^2*3,1);

% Plotting
MARGIN = .3;
AXWIDTH = 1.5;
AXFIRSTROW = 20;
AXLEFTEDGE = 3;
t = [1:size(STAs,2)]*1/WN.sum.exptParams.framerate/2;
figprefs;
for i = 1:size(STAs,2)
    STA = normfactor*(STAs(:,i)-muvect)+muvect;
    STA = reshape(STA,[nstixperside nstixperside 3]);
    axes('position',[AXLEFTEDGE+(i-1)*(AXWIDTH+MARGIN), AXFIRSTROW, AXWIDTH, AXWIDTH]);
    image(STA); axis tight;
    set(gca,'XTick',[],'YTick',[]); axis square;
    title(num2str(round(t(i)*1000)));
    if (i == 1)
        ylabel('STA');
    end
    for j = 1:size(v,2)
        PC = normfactor*(PCs(:,j,i)-muvect)+muvect;
        PC = reshape(PC,[nstixperside nstixperside 3]);
        axes('position',[AXLEFTEDGE+(i-1)*(AXWIDTH+MARGIN), AXFIRSTROW-j*(AXWIDTH+1.5*MARGIN), AXWIDTH, AXWIDTH]);
        image(PC);
        set(gca,'XTick',[],'YTick',[]); axis square;
        if (i == 1)
            ylabel(['PC',num2str(j)]);
        end
        axis tight;
    end
end

% Now plotting firing rates and surfaces
% First getting the indexs of the relevant replay trials

imageidx = find(strcmp(WN.sum.rasterCells(1,:),'synth_image'));
allimages = WN.ras(:,imageidx);
imageidxs = [];
for i = image_idx_lims(1):image_idx_lims(2)
   if ~isnan(allimages{i})
       imageidxs = [imageidxs; i];
   end
end

% Pulling out the images and putting them through the gamma functions
% to get them in normalized intensity units.
nelements = length(allimages{imageidxs(1)});
npix = nelements/3;

npixperside = sqrt(npix);
images = nan(nelements,length(imageidxs));
for i = 1:length(imageidxs)
    im = allimages{imageidxs(i)};
    for j = 1:3
        idxs = [1:npix]+npix*(j-1);
        images(idxs,i) = gammaTable1(im(idxs),j);
    end
end

% Subtracting the background and pulling out the unique images
bkgndrgbmat = repmat(bkgndrgb', npix, 1);
bkgndrgbmat = repmat(bkgndrgbmat(:), 1, length(imageidxs));
[uniqueimages, i, whichimage] = unique((images-bkgndrgbmat)','rows');
nuniqueims = max(whichimage);

% Ordering the images by projections onto two basis vectors
% "rankweights" specifies the X and Y 'coordinates' of each unique image

[u,s,v] = svd(uniqueimages');
dists = sum(v(:,1).^2+v(:,2).^2,2);
corneridx = find(dists == max(dists),1);
thetas = [atan(v(corneridx,2)./v(corneridx,1))+pi/4;...
    atan(v(corneridx,2)./v(corneridx,1))-pi/4];
theta = thetas(abs(thetas) == min(abs(thetas)));
rotmat = [cos(theta) sin(theta); -sin(theta) cos(theta)];
basis = u(:,[1 2])*rotmat' * rotmat*s([1 2], [1 2])*rotmat';
basis = pinv(basis)';  % im = b'*w so w = im*pinv(b)
norms = sqrt(sum(basis.^2));
basis = basis./repmat(norms, size(basis,1),1);
weights = uniqueimages*basis;
normweights = [(weights(:,1)-min(weights(:,1)))/(max(weights(:,1))-min(weights(:,1))),...
    (weights(:,2)-min(weights(:,2)))/(max(weights(:,2))-min(weights(:,2)))];
rankweights = round(normweights*(sqrt(nuniqueims)-1))+1;
nimages = max(rankweights);

% -------
% Images
% -------
AXLEFTEDGE = 3;
AXTOPEDGE = 13;
AXWIDTH = 1.25;
MARGIN = .1;
OFFSET = 0;
for i = 1:nimages(1)
    for j = 1:nimages(2)
        axes('position',[AXLEFTEDGE+(i-1)*(AXWIDTH+MARGIN) AXTOPEDGE-(j-1)*(AXWIDTH+MARGIN) AXWIDTH AXWIDTH])
        scalefactor = .5/max(abs(uniqueimages(:)));
        L = rankweights(:,1) == j & rankweights(:,2) == i;
        image(reshape(scalefactor*uniqueimages(L,:)',[nstixperside, nstixperside, 3])+.5);
        axis tight
        set(gca,'Xtick',[],'Ytick',[],'Box','off');  
    end
end

%-----------------
% Rasters
%-----------------
AXLEFTEDGE = 11.5;

ax_handles = nan*ones(max(rankweights));
meanspikerate = nan*ones(max(rankweights));
for i = 1:nimages(1)
    for j = 1:nimages(2)
        ax_handles(i,j) = axes('position',[AXLEFTEDGE+(j-1)*(AXWIDTH+MARGIN) AXTOPEDGE-(i-1)*(AXWIDTH+MARGIN) AXWIDTH AXWIDTH]);
        idxs = find(whichimage == find(rankweights(:,1) == i & rankweights(:,2) == j));
        counter = 0;
        tmp = [];
        for idx = idxs' % looping over trials
            stimon_t = WN.trial(imageidxs(idx),stimonidx);
            numframes = WN.trial(imageidxs(idx),nframesidx);
            stimoff_t = stimon_t+numframes/framerate;
            spikes = WN.ras{imageidxs(idx),spikeidx};
            spikes(spikes < stimon_t-OFFSET | spikes > stimoff_t+OFFSET) = [];
            nsp = length(spikes);
            plot([spikes spikes]'-stimon_t,[zeros(nsp,1) .5*ones(nsp,1)]'+counter,'k-');
            counter = counter + 1;
            tmp = [tmp; nsp./(stimoff_t-stimon_t)];
        end
        if OFFSET > 0
           % plot([0 0],[0 counter],'-','LineWidth',.5,'color','white')
           % plot([1 1]*numframes/framerate,[0 counter],'-','LineWidth',.5,'color','white')
        end
        set(gca,'Xtick',[],'Ytick',[],'Box','off','Ylim',[0 counter],'Xlim',[0 numframes/framerate]+[-OFFSET OFFSET]);
        meanspikerate(i,j) = mean(tmp);
    end
end

% Coloring the backgrounds of the plots
normspikerate = (meanspikerate-min(meanspikerate(:)))./(max(meanspikerate(:))-min(meanspikerate(:)));
for i = 1:size(ax_handles,1)
    for j = 1:size(ax_handles,2)
        tmp = (normspikerate(i,j)-.5)/2;
        set(ax_handles(i,j),'color',[.5 .5 .5]+[tmp 0 -tmp])
    end
end

% Color bar
ncolors = 20;
cmap = [linspace(0,1,ncolors)' .5*ones(ncolors,1), linspace(1,0,ncolors)'];
axes('position',[AXLEFTEDGE+nimages(1)*(AXWIDTH+MARGIN) AXTOPEDGE-(nimages(2)-1)*(AXWIDTH+MARGIN) .5 nimages(2)*(AXWIDTH+MARGIN)]);
image([1:size(cmap,1)]');
colormap(cmap)
axis tight; set(gca,'Xtick',[],'YAxisLocation','right');
frates = [min(meanspikerate(:))/10:2:max(meanspikerate(:))/10]*10;
b = [min(meanspikerate(:)) max(meanspikerate(:))]'\[1 ncolors; 1 1]'; % regression min:max maps to 1:20
set(gca,'ytick',b*[frates; ones(1,length(frates))],'yticklabel',frates);
ylabel('sp/sec');


% Getting cone weights
energy = sum(STAs.^2);
tmp = reshape(STAs(:,energy == max(energy)),nstixperside^2, 3);
[u,s,v] = svd(tmp)
rgbweights = v(:,1);
coneweights = inv(M')*rgbweights


    