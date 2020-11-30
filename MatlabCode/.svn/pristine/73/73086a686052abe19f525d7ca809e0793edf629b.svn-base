% iSETbio figures
% Figures for talk at Stanford 8/20/15
% Also for a lab meeting on 8/27/15
% Contents
%
% 1) STAs that start off green and become yellow
% 2) STA for an Isosamp(LMTF) file
% 3) Rasters from an IsoSamp(LMTF) data file 
% 4) Isodetection contours in isoluminant plane
% 5) example L-M-S and L-M+S gabors
% 6) White noise stimulus (taken from Abhishek's code)
% 6.1) Showing stimulus distributions in RGB and LMS space
% 7) Detection thresholds in the isoluminant plane
% 8) Trying to denoise Abhishek's STA

%% Section 1 STAs that go from green to yellow
% V1 neuron (from Abhishek) N070815002 
% LGN neuron A052215001
% Use WNAnalysis. Save to postscript. Load into Illustrator, transform
% (scale), export to png, import into Keynote.
%% Section 2 STA for an isosamp LGN file
% A061215001 is a nice STA that goes with A061215002, which is the IsoSamp
% file (interleaved luminance and color runs - cell is more sensitive to
% high frequency color than anything else).
%% Section 3 Single trials from IsoSamp/LMF experiment
% Franerate of Quicktime movie is 25 fps?
filename = 'A061215002';
stro = nex2stro(findfile(filename));
whichtrial = 7; % L-M: 98, 125. L+M: 18,34, 60
MAKEMOVIE = 0;
scalefactor = 40;

l = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_l'));
m = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_m'));
tfs = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'tf'));

% Finding high frequency L-M
L = sign(l)~=sign(m);
maxtf = max(tfs(L));
find(L&tfs == maxtf)
% Finding L+M at the highest frequency for L-M
L = sign(l)~=sign(m);
maxtf = max(tfs(L));
find(~L&tfs == maxtf)


stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimon_t'));
stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimoff_t'));
spikeidx = strcmp(stro.sum.rasterCells(1,:),getSpikenum(stro));
spikes = stro.ras(:,spikeidx);
secsperframe = 1./stro.sum.exptParams.framerate;
stimduration = stimoff_t(whichtrial)-stimon_t(whichtrial);
nframes = round(stimduration./secsperframe);
ncycles = tfs(whichtrial)*stimduration;
nframespercycle = nframes./ncycles;

colordir = [l(whichtrial) m(whichtrial) 0];
M = [0.0761    0.1524    0.0218;
    0.0275    0.1582    0.0321;
    0.0024    0.0118    0.1220];  % Something standard; I think this is Dell 4
bkgndrgb = [.5 .5 .5];
bkgndlms = M*bkgndrgb';
gaborrgb = inv(M)*bkgndlms.*(1+colordir'); 
gaborrgb = gaborrgb-bkgndrgb';% delta

lambda = 5;
phis = linspace(0,2*pi*ncycles,nframes+1);
phis(end) = [];

envelope = [linspace(0,1,round(length(phis)/4)) ones(1,round(length(phis)/2))];
envelope = [envelope, linspace(1,0,length(phis)-length(envelope))];

figure('units','inches','position',[1 1 4 4]);
axes('units','inches','position',[0 0 4 4],'Visible','off','color',[.5 .5 .5])
if (MAKEMOVIE)
    writerObj = VideoWriter(['IsoSampStim',num2str(whichtrial)]);
    open(writerObj);
end
for i = 1:length(phis)
    im = DrawGaborEdge(bkgndrgb, gaborrgb, [0 0 0], pi, lambda, 2, 1, phis(i), 0, 0, 0, 0, .999, 5);
    image((im-.5)*envelope(i)+.5);
    set(gcf,'Color',bkgndrgb);
    set(gca,'XTick',[],'YTick',[],'visible','off');
    axis equal;
    drawnow;
    if (MAKEMOVIE)
        frame = getframe(gca);
        writeVideo(writerObj,frame);
    end
end
if (MAKEMOVIE)
    close(writerObj);
end

timestep = 0.0001;
timevect = -.20:timestep:1;
spiketimes = stro.ras{whichtrial,1}-stimon_t(whichtrial);
spiketimes(spiketimes<timevect(1)) = [];
audiovector = double(hist(spiketimes, timevect) >= 1);
p=audioplayer(audiovector,1./timestep);
play(p);
audiowrite([filename,'_sp_',num2str(whichtrial),'.wav'], audiovector,1./(2*timestep));

% Plotting spikes
figure; axes('position',[0 .75 1 .2]); hold on;
plot(timevect,audiovector,'k-');
set(gca,'Ytick',[],'Xtick',[0 stimduration],'Xticklabel',[]);
set(gca,'Xlim',[timevect(1) timevect(end)]);
print('-dpdf',['sp',num2str(whichtrial)]);

% Plotting stimulus time course
figure; axes('position',[0 .75 1 .2]); hold on;
phis = linspace(0,2*pi*ncycles,nframes+1);
phis(end) = [];
plot(linspace(0, stimduration,length(envelope)),envelope.*sin(phis),'k-');
plot([timevect(1) 0],[0 0],'k-');
plot([stimduration timevect(end)],[0 0],'k-');
set(gca,'Ytick',[],'Xtick',[0 stimduration],'Xticklabel',[]);
set(gca,'Xlim',[timevect(1) timevect(end)]);
print('-dpdf',['st',num2str(whichtrial)]);


%% Section 4 Isodetection contours in isoluminant plane
observer = nexfilepath('Charlie', 'Sedna', 'text files', 'quest.txt'); %sedna
%observer = nexfilepath('Charlie', 'Kali', 'text files', 'questCSFdata.txt'); %kali
%observer = nexfilepath('Charlie', 'CharliePsychophysics', 'Text files', 'charliePsychNotes.txt'); %Charlie
%observer = nexfilepath('nexfilelists', 'Greg', 'DTgreg.txt'); %Greg

nTrials = 20;
perfRange = []; %don't filter on the baisis of performance
[colors, sfs, data] = questBatchProcess(observer, 'mode', nTrials, perfRange);
% data is n x m x k
% n is number of color directions
% m is number of spatial frequencies
% k is number of files
% detectionContours(colors, sfs, data);

Lint1 = (colors*[0.1386   -0.1386   0.9806]' > .999); % L-M+S
Lint2 = (colors*[0.1386   -0.1386   -0.9806]' > .999); % L-M-S
Llsf = sfs == min(sfs);
Lhsf = softEq(sfs,2.2426,1);
for i = 1:2
    if (i == 1)
        Lsf = Llsf;
    else
        Lsf = Lhsf;
    end
    thresholds1 = squeeze(data(Lint1, Lsf,:));
    thresholds2 = squeeze(data(Lint2, Lsf,:));
    
    figure; axes; hold on;
    h = bar([1 2],[nanmean(thresholds1) nanmean(thresholds2)]);
    set(h,'Facecolor','white','Edgecolor','black','linewidth',3);
    plot([1 1],nanmean(thresholds1)+[nanstd(thresholds1) -nanstd(thresholds1)],'k-','linewidth',3);
    plot([2 2],nanmean(thresholds2)+[nanstd(thresholds2) -nanstd(thresholds2)],'k-','linewidth',3)
    set(gca,'Xlim',[0 3]);
end


%% Section 5 L-M-S and L-M+S gabor movies

MAKEMOVIE = 1;
M = [0.0761    0.1524    0.0218;
    0.0275    0.1582    0.0321;
    0.0024    0.0118    0.1220];  % Something standard; I think this is Dell 4
lambda = 10; % 5 is good
nframespercycle = 50;
ncycles = 1.5;
bkgndrgb = [.7 .7 .7];
for i = 1:4
    if (i == 1)
        gaborlms = [ 0.1386   -0.1386    0.9806]/20;  % L-M+S deltas
        moviefilename = 'gaborL-M+S';
    elseif (i == 2)
        gaborlms = [ 0.1386   -0.1386    -0.9806]/20;  % L-M-S deltas
        moviefilename = 'gaborL-M-S';
    elseif (i == 3)
        gaborlms = [0.7071   -0.7071   0]/50;  % L-M deltas
        moviefilename = 'gaborL-M';                
    elseif (i == 4)
        gaborlms = [0.7071   0.7071   0]/15;  % L+M deltas
        moviefilename = 'gaborL+M';        
    else
        gaborlms = [0 0 1]/50;  % S-isolating deltas
        moviefilename = 'gaborS';        
    end
    gaborrgb = inv(M)*gaborlms';
    
    phis = linspace(0,2*pi*ncycles,nframespercycle*ncycles+1);
    phis(end) = [];
    figure('units','inches','position',[1 1 4 3]);
    axes('units','inches','position',[0 0 4 3],'Visible','off','color',[.5 .5 .5])
    
    envelope = [linspace(0,1,round(length(phis)/4)) ones(1,round(length(phis)/2))];
    envelope = [envelope, linspace(1,0,length(phis)-length(envelope))];
    
    
    if (MAKEMOVIE)
        writerObj = VideoWriter(moviefilename);
        open(writerObj);
    end
    for j = 1:length(phis)
        im = DrawGaborEdge(bkgndrgb, gaborrgb, [0 0 0], pi, lambda, 2, 1, phis(j), 0, 0, 0, 0, .999, 5);
        bkgndrgbmat = repmat(permute(bkgndrgb,[3,1,2]),[size(im,1),size(im,2),1]);
        image((im-bkgndrgbmat)*envelope(j)+bkgndrgbmat);
        set(gcf,'Color',bkgndrgb);
        set(gca,'XTick',[],'YTick',[],'visible','off');
        axis equal;
        drawnow;
        if (MAKEMOVIE)
            frame = getframe(gca);
            writeVideo(writerObj,frame);
        end
    end
    if (MAKEMOVIE)
        close(writerObj);
    end
end

%% Section 6 white noise stimulus (from Abhishek)
% Upgraded to show binary cone noise too.

USESUBUNITMASK = 1;
NOISETYPE = 1; % 1 = Gaussian RGB, 2 = Binary cone noise, 3 = LM noise

M =[0.0608    0.1219    0.0175;...
    0.0220    0.1266    0.0257;...
    0.0019    0.0095    0.0976]; % From Dell4BitsCal(5)

if (NOISETYPE == 1)
    writerObj = VideoWriter('gunnoise.mp4');
elseif (NOISETYPE == 2)
    CCs = [.05 .05 .5];
    bkgndlms = M*[.5 .5 .5]';
    writerObj = VideoWriter('conenoise.mp4');
elseif (NOISETYPE == 3)
    CCs = [.5 .05 0]; % L+M and L-M, NOT L and M
    bkgndlms = M*[.5 .5 .5]';
    writerObj = VideoWriter('LMnoise.mp4');
end

if (USESUBUNITMASK)
    if ~exist('stro')
        stro = nex2stro;
    end
    mask = stro.ras{end,end};
    nstixperside = sqrt(length(mask));
    mask = reshape(mask, nstixperside, nstixperside);
else
    nstixperside = 11;
end
close all;
nframes = 200;
figure;
axes;
set(gcf,'Color',[.5 .5 .5]);
set(gca,'XTick',[],'YTick',[],'visible','off');
axis square;

open(writerObj);

for i = 1:nframes
    if (NOISETYPE == 1)
        a = uint8(normrnd(128,30,nstixperside,nstixperside,3));
    elseif (NOISETYPE == 2)
        a = 2*unidrnd(2,nstixperside,nstixperside,3)-3; % -1 and 1
        for i = 1:3
            a(:,:,i) = bkgndlms(i)*(a(:,:,i)*CCs(i)+1); % cone excitations now
        end
        b = inv(M)*reshape(a,[nstixperside*nstixperside, 3])';
        a = uint8(reshape(b',size(a))*255);
    elseif (NOISETYPE == 3)
        tmp = unidrnd(4,nstixperside,nstixperside);
        a = zeros(size(tmp,1),size(tmp,2),3);
        a(:,:,1) = (tmp == 1)*CCs(1);
        a(:,:,1) = a(:,:,1)+(tmp == 2)*-CCs(1); 
        a(:,:,1) = a(:,:,1)+(tmp == 3)*CCs(2);
        a(:,:,1) = a(:,:,1)+(tmp == 4)*-CCs(2);
        a(:,:,2) = (tmp == 1)*CCs(1);
        a(:,:,2) = a(:,:,2)+(tmp == 2)*-CCs(1); 
        a(:,:,2) = a(:,:,2)+(tmp == 3)*-CCs(2);
        a(:,:,2) = a(:,:,2)+(tmp == 4)*CCs(2);
        
        for i = 1:3
            a(:,:,i) = bkgndlms(i)*(a(:,:,i)+1);
        end
        b = inv(M)*reshape(a,[nstixperside*nstixperside, 3])';
        a = uint8(reshape(b',size(a))*255);
    end
    if (USESUBUNITMASK)
        mask(mask == 0) = max(mask(:))+1;
        colors = squeeze(a([1:max(mask)],1,:));
        colors(end,:) = [128 128 128];
        a = uint8(zeros(size(a)));
        for i = 1:size(colors,1)
           tmp = mask == i;
           a(:,:,1) = a(:,:,1)+uint8(tmp)*colors(i,1);
           a(:,:,2) = a(:,:,2)+uint8(tmp)*colors(i,2);
           a(:,:,3) = a(:,:,3)+uint8(tmp)*colors(i,3);
        end
    end
    image(a);
    colormap(gray)
    axis square;
    set(gcf,'Position',[566 375 500 500]);
    set(gca,'XTick',[],'YTick',[],'Visible','off');
    axis image;
    frame = getframe(gca);
    writeVideo(writerObj,frame);
end

close(writerObj);

%% Section 6.1 White noise stimuli in RGB and LMS space
%
MAKEMOVIE = 1;
PLOTCONENOISE = 1;
PLOTGUNNOISE = 1;
M =[0.0608    0.1219    0.0175;...
    0.0220    0.1266    0.0257;...
    0.0019    0.0095    0.0976]; % From Dell4BitsCal(5)
CCs = [.05 .05 .5];
bkgndlms = M*[.5; .5; .5];
lms_binary = [];
for i = 1:3
   lms_binary(i) = bkgndlms(i)*(CCs(i)+1); % cone excitations now
end
lms_wrt_bkgnd = lms_binary-bkgndlms';
corners = 2*fullfact([2 2 2])-3;

nstim = 1000;
rgb_gauss = normrnd(.5,.2,nstim,3);
lms_gauss = rgb_gauss*M';

figure(1); axes; hold on; grid on;
set(gcf,'Color',[.5 .5 .5]);
if (PLOTGUNNOISE)
    plot3(rgb_gauss(:,1)-.5,rgb_gauss(:,2)-.5,rgb_gauss(:,3)-.5,'k.');
end
if (PLOTCONENOISE)
    for i = 1:size(corners,1)
        rgb_binary = inv(M)*(corners(i,:).*lms_wrt_bkgnd)';
        h = plot3(rgb_binary(1), rgb_binary(2), rgb_binary(3),'rs');
        set(h,'MarkerFaceColor','red');
    end
end
set(gca,'Xlim',[-.5 .5],'Ylim',[-.5 .5],'Zlim',[-.5 .5]);
set(gca,'View',[45, 30]); axis vis3d;
xlabel('R','FontSize',20); ylabel('G','FontSize',20); zlabel('B','FontSize',20);

figure(2); axes; hold on; grid on;
set(gcf,'Color',[.5 .5 .5]);
if (PLOTGUNNOISE)
    plot3(lms_gauss(:,1)-bkgndlms(1),lms_gauss(:,2)-bkgndlms(2),lms_gauss(:,3)-bkgndlms(3),'k.');
end
if (PLOTCONENOISE)
    for i = 1:size(corners,1)
        h = plot3(corners(i,1)*lms_wrt_bkgnd(1), corners(i,2)*lms_wrt_bkgnd(2), corners(i,3)*lms_wrt_bkgnd(3),'rs');
        set(h,'MarkerFaceColor','red');
    end
end
set(gca,'View',[45, 30]); axis vis3d;
xlabel('L','FontSize',20); ylabel('M','FontSize',20); zlabel('S','FontSize',20);

views = 45:3:360+45;
views = [views; 30*ones(1,size(views,2))]';

if (MAKEMOVIE)
    for j = 1:2
        figure(j)
        writerObj = VideoWriter(['cloud',num2str(j)','.mp4']);
        open(writerObj);
        for i = 1:size(views,1)
            view(views(i,1),views(i,2))
            frame = getframe(gca);
            writeVideo(writerObj,frame);
        end
        close(writerObj);
    end
end


%% Section 7 detection thresholds in the isoluminant plane

observer = nexfilepath('Charlie', 'Sedna', 'text files', 'quest.txt'); %sedna
nTrials = 20;
perfRange = []; %don't filter on the baisis of performance
[colors, sfs, data] = questBatchProcess(observer, 'mode', nTrials, perfRange);
detectionContours(colors, sfs, data);

%% Section 8
% Trying to denoise STAs
load('/Users/greghorwitz/Library/Containers/com.apple.mail/Data/Library/Mail Downloads/4FBDA399-D1E3-4F75-B07F-F9A3DA650CB1/STAs.mat');

normfact = 1/(2*max(abs(STAs(:))));
for i = 1:size(STAs,2)
    subplot(4,4,i);
    image(normfact*reshape(STAs(:,i),10,10,3)+.5);
    axis square;
    set(gca,'Xtick',[],'Ytick',[]);
end

avgSTA = mean(STAs(:,[4:7]),2);
normfact =  1/(2*max(abs(avgSTA)));
image(normfact*reshape(avgSTA,10,10,3)+.5);