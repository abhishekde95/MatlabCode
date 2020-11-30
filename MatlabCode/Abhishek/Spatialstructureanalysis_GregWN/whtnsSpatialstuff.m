%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Stuff for testing out color3.stt
%  the paradigm with the spatially
%  varying stimulus.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

name = '/home/horwitz/Cortex/Data/JUNK.1';
[time_arr,event_arr,eog_arr,header,trialcount]  = get_data(name);
[foo, index, raster, genparams] = readcort2 (name);
indices_whtnsSpatial;

% The first thing I have to do is make sure the random number seeds and
% random numbers are working correctly.  Much of this stuff was copied
% from above (when I was testing out the uniform color patch paradigm)
% and is modified here to work with the spatially-varying paradigm.

nframes = sum(event_arr == 119);
seeds = [];
% RGBs = nan*ones(sum(nframes),3);
times = [];
rowcounter = 1;
for i = 1:trialcount
   i
   % Getting the random number seed
   for j = 0:3   
      a = find(event_arr(:,i) == 120+j);
      L = event_arr(:,i) ~= 102;
      tmp = min(find(find(L)>a));
      seed(j+1) = event_arr(tmp,i) - 1000;
   end
   seed = seed(1)+seed(2)*2^8+seed(3)*2^16+seed(4)*2^24;
   seeds = [seeds, seed];

   stimon = find(event_arr(:,i) == 23);
   stimoff = find(event_arr(:,i) == 24);
   a = find(event_arr(:,i) == 119);
   times = [times; diff(time_arr(a,i))];
%   for j = 1:length(a)
%      rgb = [];
%      k = 1;
%      while(length(rgb) < 3)
%	 currcode = event_arr(a(j)+k,i)
%         if (currcode > 500)
%	    rgb = [rgb, currcode - 1000];
%         end
%         k=k+1;
%      end
%      RGBs(rowcounter,:) = rgb;
%      rowcounter = rowcounter+1;
%   end
end

str = ['[monID,invGammaTable] = ReadGam(''/home/horwitz/PR650/SONYF500_2706097/gamma'');']
eval(str);

BKGNDRGB = [genparams(GEN_BKGNDR),
            genparams(GEN_BKGNDG),
            genparams(GEN_BKGNDB)];
for i = 1:3
  lower = min(find(invGammaTable(:,i) == BKGNDRGB(i)));
  upper = max(find(invGammaTable(:,i) == BKGNDRGB(i)));
  RGBmid(i) = floor((lower+upper)/2);  
end

% Making the master lookup table
invgamres = genparams(GEN_NGAMELS)*10;
 maxgamidx = invgamres-1;
tmp = linspace(1,invgamres,invgamres)';
tmp = tmp/(maxgamidx+2);
tmp = EJnorminv(tmp);
RADIUS = genparams(GEN_STDEV)/100;
tmp = round(tmp*RADIUS*maxgamidx);

bigtable =[];
for gun = 1:3
   bigtable(:,gun) = tmp+RGBmid(gun)-1;
   bigtable(:,gun) = invGammaTable(bigtable(:,gun)+1,gun);
end
plot(bigtable)

% There , we've made bigtable....
% Format of newRGBs is (frame,pixel,gun)
newRGBs = nan*ones(sum(nframes),genparams(GEN_STIMNELEMX)*genparams(GEN_STIMNELEMY),3);
MAXRAND = 65535;
for trial = 1:trialcount
   trial
   if (nframes(trial) >= 1)
      % Getting the random numbers for this trial
      nrands = 3*nframes(trial)*genparams(GEN_STIMNELEMX)*genparams(GEN_STIMNELEMY);
      str = ['!/home/horwitz/src/random2 ',num2str(nrands),' ',num2str(seeds(trial)),'  > randout']; 
      eval(str);
      eval('load randout');
      eval('!rm randout');

      % Getting the dropped random number
      tmp = zeros(2,1);
      a = find(event_arr(:,trial) == 124);
      tmp(1) = event_arr(a+1,trial)-1000;
      a = find(event_arr(:,trial) == 125);
      tmp(2) = event_arr(a+1,trial)-1000;
      randomnum = tmp(1)+tmp(2)*2^8
      if (randomnum ~= randout(1))
         error('Mismatched random numbers!');
      end
      counter = 1;
      for i = 1:nframes(trial)
         for j = 1:genparams(GEN_STIMNELEMX)*genparams(GEN_STIMNELEMY)
            r = randout(counter)*maxgamidx/MAXRAND;
            g = randout(counter+1)*maxgamidx/MAXRAND;
            b = randout(counter+2)*maxgamidx/MAXRAND;

            rgb = round([r g b]);
   
            R = bigtable(rgb(1)+1,1);
            G = bigtable(rgb(2)+1,2);
            B = bigtable(rgb(3)+1,3);

            newRGBs(i,j,:) = [R G B];
            counter = counter+3;
         end
      end
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Just trying some basic
%   things like looking at
%   the STA, etc.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	

name = '/home/horwitz/Cortex/Data/MAX/M001.1'
paranum = get_paradigm(name,1000);
[foo, index, raster, genparams,e1h,e1v,e1t] = readdata (name,1);

[out] = getWhtnsStats(foo,index,raster,genparams,20,0,'spike',...
                      paranum,'STAmod');

STA = out{1};
STV = out{2};
minval = min(STA(:));
maxval = max(STA(:));
minval = min([minval -maxval]);
maxval = max([maxval -minval]);
figure;
for i = 1:size(STA,3)
   im = squeeze(STA(:,:,i));
   im = reshape(im,sqrt(size(im,1)),sqrt(size(im,1)),size(im,2));
   im = (im-minval)/(maxval-minval);
   subplot(ceil(sqrt(size(STA,3))),ceil(sqrt(size(STA,3))),i)
   image(im);
   colormap(hot);
   axis('square')
   set(gca,'XTickLabel',[],'YTickLabel',[]);
end
title(max(max(max(abs(STA)))));
[x,y] = ginput;
x = round(x); y = round(y);

if (paranum == 60)
   indices_whtnsSpatial2;
else
   indices_whtnsSpatial;
end

idxs = y+genparams(GEN_STIMNELEMX)*(x-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Getting Allstim for a 
%  ROI so we can do the 
%  analysis of non-linearity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maxT = 20;
nframes = ceil(sum(index(:,SIN_STIMOFF)-index(:,SIN_STIMON))/10);
AllStimDim =  [nframes 1 3 maxT];
AllStimmod('init',{idxs,AllStimDim});

[out] = getWhtnsStats(foo,index,raster,genparams,20,0,'all',...
                      paranum,'AllStimmod');

AllStim = out{1};
Lspike = out{2};

%  Now we can go use the code in EeroStuff
%  or in RFsearch


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Trying eye position correction
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
name = '/home/horwitz/Cortex/Data/CORN/C146.5'
paranum = get_paradigm(name,1000);
[foo, index, raster, genparams,e1h,e1v,e1t] = readcort2 (name,103);

[out] = getWhtnsStats2(foo,index,raster,genparams,e1h,e1v,e1t,13,0,'spike',...
                      paranum,'STAmod');

STA = out{1};
minval = min(STA(:));
maxval = max(STA(:));
minval = min([minval -maxval]);
maxval = max([maxval -minval]);
figure;
for i = 1:size(STA,3)
   im = squeeze(STA(:,:,i));
   im = reshape(im,sqrt(size(im,1)),sqrt(size(im,1)),size(im,2));
   im = (im-minval)/(maxval-minval);
   subplot(ceil(sqrt(size(STA,3))),ceil(sqrt(size(STA,3))),i)
   image(im);
   colormap(hot);
   axis('square')
   set(gca,'XTickLabel',[],'YTickLabel',[]);
end
title(max(max(max(abs(STA)))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%  Trying out ways of computing covariance
%  matrices for each pixel in the stimulus
%  array.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

name = '/home/horwitz/Cortex/Data/CORN/C127.2c3'
paradigmcd = get_paradigm(name,1000);
[foo, index, raster, genparams] = readcort2 (name,102);

maxT = 11;
frameskip = 0;


opt = 'spike';
module = 'STCOVmod4';
foo = foo([1:50],:);
raster = raster([1:50],:);
index = index([1:50],:);
getWhtnsSetup;

feval(module,'init',[genparams(GEN_STIMNELEMX)*genparams(GEN_STIMNELEMY) 3 6 6]);
%out = feval(module,'return');
[out] = getWhtnsStats3(foo,index,raster,genparams,maxT,0,'spike',paradigmcd,module);

STCOV = out{2};
OPTIONS.issym = 1;
OPTIONS.disp = 1;

tic
[v,d,flag] = eigs(STCOV,2,OPTIONS)
toc

v = [out{1} v];
figure;
for j = 1:size(v,2)
   vector = reshape(v(:,j),[3,64,maxT]);
   vector = permute(vector,[2,1,3]);
   minval = min(vector(:));
   maxval = max(vector(:));
   minval = min([minval -maxval]);
   maxval = max([maxval -minval]);
   for i = 1:size(vector,3)
      im = squeeze(vector(:,:,i));
      im = reshape(im,sqrt(size(im,1)),sqrt(size(im,1)),size(im,2));
      im = (im-minval)/(maxval-minval);
      axes('Units','normalized','position',[i/(maxT+2) j/(size(v,2)+2) .05 .05]);
      axis square;
      image(im);
      colormap(hot);
      axis('square')
      set(gca,'XTickLabel',[],'YTickLabel',[]);
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Eye position correction stuff
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

name = '/home/horwitz/Cortex/Data/CORN/C153.5'
paradigmcd = get_paradigm(name,1000);
[foo, index, raster, genparams, e1h, e1v, e1t] = readcort2 (name,102);

maxT = 13;
[SpikeStim,TrialVect,TimeVect] =...
              getWhtnsStim(foo,index,raster,genparams,maxT,0,'spike',paradigmcd);

indices_whtnsSpatial;
arrayHeight = genparams(GEN_STIMH)/10;
arrayWidth = genparams(GEN_STIMW)/10;
nx = genparams(GEN_STIMNELEMX)
ny = genparams(GEN_STIMNELEMY)
framedur = getGlobalParams('framedur');

degperblock = arrayHeight/ny
degperblock = arrayWidth/nx  % These things should give the same answer
fixwinheight = 1;

maxshift = ceil(fixwinheight/degperblock);  % Maximum number of pixels we'll ever have to shift
template = reshape([1:nx*ny],nx,ny);  % Representation of the bitmap
templates = zeros(nx,ny,maxshift,maxshift);
shiftmap = [-maxshift:maxshift];
for i = 1:length(shiftmap)
   for j = 1:length(shiftmap)
      templates(:,:,i,j) = scroll(template,[shiftmap(i),shiftmap(j)]);
   end
end

NewSpikeStim = zeros(size(SpikeStim));
counter = 1;

eyesamptime = getGlobalParams('eyesamplerate');
nframes = size(SpikeStim,4);
integtime = framedur*nframes;
neyesamps = round(integtime/eyesamptime);  % Number of eye position samples per stimulus
for i = 1:length(TrialVect)
   i
   trial = TrialVect(i);
   t_end = TimeVect(i)-e1t(trial)+framedur/2; % End time of last frame relative to start of eye collection (start of fixation)
   samp_end = floor(t_end/eyesamptime);  % Using "floor" here is a bit odd, but rounding can push me outside of the data
   Hs = e1h(trial, [samp_end-neyesamps+1:samp_end]);  % Hs go from begining of first frame to end of last frame
   Vs = e1v(trial, [samp_end-neyesamps+1:samp_end]);
   Hs = Hs-nanmedian(e1h(trial,:));
   Vs = Vs-nanmedian(e1v(trial,:));
   eyesampbins = round(linspace(0,neyesamps,nframes+1));

   Stim = squeeze(SpikeStim(i,:,:,:));
   for j = 1:nframes
      hshift = median(Hs([eyesampbins(j)+1:eyesampbins(j+1)]));
      vshift = median(Vs([eyesampbins(j)+1:eyesampbins(j+1)]));
      hshift = round(hshift/degperblock);
      vshift = round(vshift/degperblock);
      hshift  = max([hshift -maxshift]);
      hshift  = min([hshift maxshift]);
      vshift  = max([vshift -maxshift]);
      vshift  = min([vshift maxshift]);
      tmp = templates(:,:, find(shiftmap == -hshift), find(shiftmap == vshift));  % See scroll.m sign conventions
      NewSpikeStim(counter,:,:,nframes-j+1) = Stim(tmp(:),:,nframes-j+1);
   end
   counter = counter + 1;
end
if (counter < size(NewSpikeStim,1))
   NewSpikeStim(counter:end,:,:,:) = [];
end

STA1 = squeeze(mean(SpikeStim,1));
STA2 = squeeze(mean(NewSpikeStim,1));

n1 = size(SpikeStim,1);
STV1 = (n1*squeeze(sum(SpikeStim.^2,1))-((n1*STA1).^2))/(n1*(n1-1));
n2 = size(NewSpikeStim,1);
STV2 = (n2*squeeze(sum(NewSpikeStim.^2,1))-((n2*STA2).^2))/(n2*(n2-1));

figure
subplot(2,2,1)
hist(STA1(:),100);
title(num2str([min(STA1(:)) max(STA1(:))]));
subplot(2,2,2)
hist(STA2(:),100);
title(num2str([min(STA2(:)) max(STA2(:))]));
subplot(2,2,3)
hist(STV1(:),100);
title(num2str([min(STV1(:)) max(STV1(:))]));
subplot(2,2,4)
hist(STV2(:),100);
title(num2str([min(STV2(:)) max(STV2(:))]));

for i = 1:2;
   if i == 1
      tit = 'Before correction';
      STA = STA1;
   else
      tit = 'After correction';
      STA = STA2;
   end
   minval = min(STA(:));
   maxval = max(STA(:));
   if (minval < 0)
      minval = min([minval -maxval]);
      maxval = max([maxval -minval]);
   end
   figure;
   for i = 1:size(STA,3)
      im = squeeze(STA(:,:,i));
      im = reshape(im,sqrt(size(im,1)),sqrt(size(im,1)),size(im,2));
      im = (im-minval)/(maxval-minval);
      subplot(ceil(sqrt(size(STA,3))),ceil(sqrt(size(STA,3))),i)
      image(im);
      colormap(hot);
      axis('square')
      set(gca,'XTickLabel',[],'YTickLabel',[]);
   end
   title(tit);
end

% Looking at noise in the eye position traces

figure; hold on;
h = nanmean(e1h');
v = nanmean(e1v');
plot(h,'g-')
plot(v,'m-')

a = xcorr(h-mean(h));
b = xcorr(v-mean(v));
figure; hold on;
plot(b,'m-');
plot(a,'g-');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Doing Fourier transforms on STA to estimate
%  preferred spatial frequency.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

name = '/home/horwitz/Cortex/Data/KIMBALL/K081.3'
paradigmcd = get_paradigm(name,1000);
[foo, index, raster, genparams, e1h, e1v, e1t] = readcort2 (name,102);

maxT = 13;
[out] = getWhtnsStats(foo,index,raster,genparams,maxT,0,'spike',...
                      paradigmcd,'STAmod');

STA = out{1};
tmp = reshape(STA,size(STA,1)*size(STA,2),size(STA,3));
tmp = max(tmp);
peakT = find(tmp == max(tmp));find(tmp == max(tmp));
im = squeeze(STA(:,:,peakT));
im = sum(im,2);
im = reshape(im, [sqrt(size(STA,1)) sqrt(size(STA,1))]);
imagesc(im);

fftim = fftshift(abs(fft2(im)));
fftim = fftim([1:5],:)
% Preferred spatial frequencies below 
% Shifted so 0 is DC, 1 is 1 cycle, 2 is 2 cycles, etc.
peakfft = max(max(fftim));
p_omega_h = find(max(fftim) == peakfft)-5
p_omega_v = 5-find(max(fftim' == max(max(fftim'))))

indices_whtnsSpatial;
arrayHeight = genparams(GEN_STIMH)/10;
arrayWidth = genparams(GEN_STIMW)/10;

p_sf = sqrt(p_omega_v^2+p_omega_h^2)/arrayWidth  % In cycles/degree
p_orient = atan2(p_omega_v, p_omega_h)*180/pi+90 % In degrees


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Getting a bunch of STAs and STCOVs
%  This is just so that I can set
%  the script running over night and,
%  in the morning, I should be in a good
%  position to make basis vectors and
%  find stimulus projections.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

path = '/home/horwitz/Cortex/Data';
files = /CORN/c127.2c3'
paradigmcd = get_paradigm(name,1000);
[foo, index, raster, genparams] = readcort2 (name,102);

maxT = 11;
frameskip = 0;

opt = 'spike';
module = 'STCOVmod3';
foo = foo([1:50],:);
raster = raster([1:50],:);
index = index([1:50],:);
getWhtnsSetup;

feval(module,'init',[genparams(GEN_STIMNELEMX)*genparams(GEN_STIMNELEMY) 3 maxT]);
%out = feval(module,'return');
[out] = getWhtnsStats3(foo,index,raster,genparams,maxT,0,'spike',paradigmcd,'STCOVmod3');

STCOV = out{2};
OPTIONS.issym = 1;
OPTIONS.disp = 1;

tic
[v,d,flag] = eigs(STCOV,5,OPTIONS)
toc

v = [out{1} v];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Obsolete stuff below this point.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Taking a look at subunit 
%  interactions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Should have created the STA map by now
% Click on two squares, then hit return
figure(1);
[x,y] = ginput;
x = round(x); y = round(y);
idxs = genparams(GEN_STIMNELEMX)*(y-1)+x;
if (length(idxs) ~= 2)
   error('Should have selected exactly two squares');
end

maxT = 20;
nframes = ceil(sum(index(:,SIN_STIMOFF)-index(:,SIN_STIMON))/10);
AllStimDim =  [nframes 2 3 maxT];
AllStimmod('init',{idxs,AllStimDim});

[out] = getWhtnsStats(foo,index,raster,genparams,maxT,0,'all',...
                      paranum,'AllStimmod');

AllStim = out{1};
Lspike = out{2};

Template1 = squeeze(mean(AllStim(Lspike > 0,1,:,:),1));
Template2 = squeeze(mean(AllStim(Lspike > 0,2,:,:),1));

x = vectorizestim(AllStim(:,1,:,:));
y = vectorizestim(AllStim(:,2,:,:));

dotprods1 = x*Template1(:);
dotprods2 = y*Template2(:);

nbins = 7;
xbins = linspace(min(dotprods1)/2,max(dotprods1)/2,nbins+1)
ybins = linspace(min(dotprods2)/2,max(dotprods2)/2,nbins+1)
% these are bin boundaries.
mat = zeros(nbins, nbins);
for i = 1:nbins
   for j = 1:nbins
      L = dotprods1 > xbins(i) & dotprods1 < xbins(i+1) &...
          dotprods2 > ybins(j) & dotprods2 < ybins(j+1);
      response = Lspike(L);
      response = sum(response > 0);
      mat(j,i) = response/sum(L);
   end
end
figure
imagesc(flipud(mat))

cntrlines = [0 logspace(-2,0,25)];
figure
contourf(mat,cntrlines);
axis('square')
xlabel('Proj. onto template 1');
ylabel('Proj. onto template 2');


figure; set(gcf,'DefaultAxesColorOrder',[1 0 0; 0 1 0; 0 0 1]);
subplot(2,1,1);
plot(Template1');
title('STA at location 1');
subplot(2,1,2);
plot(Template2');
title('STA at location 2');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Automated procedure for finding the 
%  dimensions of the receptive field
%  in both space and time.  The algorithim
%  will be to run through the data doing
%  tests on the mean and variance of
%  each pixel and gun.  Each pixel will
%  then be considered "in the RF" if 
%  at least one gun is significant in at least
%  one of the two tests.  Then, we'll go
%  through and find the largest region of
%  continguous pixels.  This will be the
%  RF.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MAXRAND = 65535;
ntimepts = 20;
name = '/home/horwitz/Cortex/Data/CORN/C127.2c3';
[foo, index, raster, genparams,e1h,e1v,e1t] = readdata(name,102);
paradigmcd = get_paradigm(name,1000);

% Doing the work
[out] = getWhtnsStats(foo,index,raster,genparams,ntimepts,0,'spike',paradigmcd,'STAmod');
STA = out{1}; 
STV = out{2}; 
n = out{3};

x = [0:MAXRAND];
if (paradigmcd == 60)
   indices_whtnsSpatial2;
   ngamels = genparams(GEN_NGAMELS1)+genparams(GEN_NGAMELS2)*2^8;
   x = bitand(x,ngamels-1);
else
   indices_whtnsSpatial;
   ngamels = genparams(GEN_NGAMELS)*10;
   x = round(x*(ngamels-1)/MAXRAND);
end
% Setting up for the statistics

monserial = genparams(GEN_MONSERIAL1)+...
            genparams(GEN_MONSERIAL2)*2^8+...
            genparams(GEN_MONSERIAL3)*2^16;
bigtable = getbigtable(genparams(GEN_BKGNDR),genparams(GEN_BKGNDG),genparams(GEN_BKGNDB),genparams(GEN_STDEV)/100,ngamels,monserial);
cal = getGlobalParams('cal',monserial);
weights = FindModelWeights(cal.P_ambient,cal.P_device)';
gam = cal.gammaTable+repmat(weights,size(cal.gammaTable,1),1);
transformedbigtable = zeros(size(bigtable));
for gun = 1:3
   transformedbigtable(:,gun) = ...
         gam(bigtable(:,gun)+1,gun);
end

LinGunVals = zeros(MAXRAND+1,3);
for gun = 1:3
   LinGunVals(:,gun) = ...
         transformedbigtable(x+1,gun);
end
ExpectedValues = mean(LinGunVals);
Variances = var(LinGunVals);

%%%%%%% End of statistics set up %%%%%%%

s = size(STA);
pmat1 = zeros(s);   % p values from STA
pmat2 = zeros(s);   % p values from STV
tmat1 = zeros(s);   % t values from STA
tmat2 = zeros(s);   % t values from STV
for gun = 1:size(STA,2)
   for timept = 1:size(STA,3)
      % STA
      im = STA(:,gun,timept);
      im = im./(sqrt(Variances(gun)/n));
      tmat1(:,gun,timept) = im; 
      pmat1(:,gun,timept) = (1-normcdf(abs(im')));
      % STV
      im = STV(:,gun,timept);
      im = (n-1)*im/Variances(gun);
      pmat2(:,gun,timept) = chi2cdf(im,n-1);
      [m,v] = chi2stat(n-1);
      tmat2(:,gun,timept) = (im-m)./sqrt(v);
   end
end

alpha = 0.0001;
pmat = zeros(size(pmat1,1),size(pmat1,3));
for i = 1:ntimepts
   tmp = min(pmat1(:,:,i),pmat2(:,:,i));
   tmp = pmat1(:,:,i);
   pmat(:,i) = any(tmp<alpha,2);
end

pmat = reshape(pmat,sqrt(size(pmat,1)),sqrt(size(pmat,1)),size(pmat,2));
pmat = permute(pmat,[2 1 3]);
for time = 1:ntimepts
   subplot(ceil(sqrt(ntimepts)),ceil(sqrt(ntimepts)),time);
   colormap(gray(2));
   image(pmat(:,:,time)+1);
   axis('square');
end

% Number of false positives that we
% expect by chance (assuming everything's
% independent, which of course it's not).

2*prod(size(pmat1))*alpha

% Now running through and eliminating pixels
% that aren't neighboring other significanct
% pixels.
idxs = find(pmat);
clusterids = [1:length(idxs)]';
for i = 1:length(idxs)
   [i1, i2, i3] = ind2sub(size(pmat),idxs(i));
   i1 = i1+[-1 0 1];
   i2 = i2+[-1 0 1];
   i3 = i3+[-1 0 1];
   i1 = i1(i1>0 & i1<=size(pmat,1));
   i2 = i2(i2>0 & i2<=size(pmat,2));
   i3 = i3(i3>0 & i3<=size(pmat,3));
   neighbors = [];
   for j = fullfact([length(i1), length(i2), length(i3)])'
      if (pmat(i1(j(1)),i2(j(2)),i3(j(3))))
         neighbors = [neighbors; sub2ind(size(pmat),i1(j(1)),i2(j(2)),i3(j(3)))];
      end
   end
   L = logical(ismember(idxs,neighbors));
   clusterids(L) = min(clusterids(L));
end
[a,b] = hist(clusterids,unique(clusterids))
whichcluster = b(find(a==max(a)));
pmat(idxs(clusterids ~= whichcluster)) = 0;

figure;
for time = 1:ntimepts
   subplot(ceil(sqrt(ntimepts)),ceil(sqrt(ntimepts)),time);
   colormap(gray(2));
   image(pmat(:,:,time)+1);
   axis('square');
end

[i1, i2, i3] = ind2sub(size(pmat),find(pmat));

[min(i3) max(i3)]    % Latency and duration
RFperimeter = convhull(i2,i1);  % Convhull wants x,y not row,column
if (isempty(RFperimeter))  % Colinear RF
   RFperimeter = unique([i1, i2],'rows');
end 

allpts = fullfact([8 8]);
L = inpolygon(allpts(:,2),allpts(:,1),i2(RFperimeter),i1(RFperimeter));
% i1(RFperimeter) are the "xs"
% i2(RFperimeter) are the "ys"

RFmap = zeros(8,8);
RFmap(L>0) = 1;
idxs = find(RFmap);
tmp(sub2ind(size(RFmap),i1(RFperimeter),i2(RFperimeter))) = 2;  % Checking vertices
figure;
image(RFmap+1);
colormap(gray(3));

figure
image(sum(pmat,3))

tmp = vectorizestim(permute(STA(idxs,:,:),[1 3 2]));
[pc,score,latent] = princomp(tmp);
figure; set(gcf,'DefaultAxesColorOrder',[1 0 0; 0 1 0; 0 0 1]);
subplot(1,3,1)
plot(latent./sum(latent),'k.');
subplot(1,3,2)
plot(reshape(pc(:,1),20,3))
subplot(1,3,3)
plot(reshape(pc(:,2),20,3))

for i = 1:2
   tmp = RFmap;
   tmp(idxs) = score(:,i);
   figure;
   imagesc(tmp);
   colormap(gray);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  How many spikes went into an STA calculation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

name = '/home/horwitz/Cortex/Data/KIMBALL/K076.3c4';
maxT = 20; frameskip = 0;
[foo, index, raster, genparams,e1h,e1v,e1t] = readdata(name,102);
paradigmcd = get_paradigm(name,1000);
indices_whtnsSpatial2;
early = repmat(index(:,SIN_STIMON)+maxT+frameskip,1,size(raster,2));
late = repmat(index(:,SIN_STIMOFF),1,size(raster,2));
sum(sum(raster>early & raster < late))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Breaking down a RF into color and (time/space)
%  Plot the points in 3-D.  Is the surround really 
%  spectrally opposed to the center (do the points
%  lie on a line that passes through the origin?)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

maxT = 20;
alpha = 0.0001;
name = '/home/horwitz/Cortex/Data/CORN/C129.2c6';
[foo, index, raster, genparams,e1h,e1v,e1t] = readdata(name,102);
paranum = get_paradigm(name,1000);
[ExpectedValues, Variances] = StimStats(genparams, paranum);

RFmod('init',{ExpectedValues,Variances,alpha});
[out] = getWhtnsStats(foo,index,raster,genparams,maxT,0,'spike',paranum,'RFmod');

n = out{1};
framestart = out{2};
framestop = out{3};
RFmap = out{4};
STA = out{5};
clear out;

idxs = find(RFmap);

data = [];
for frame = framestart:framestop
   data = [data; STA(idxs,:,frame)];
end

figure
plot3(data(:,1),data(:,2),data(:,3),'k.')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Comparing model fits (additive vs. multiplicative for
%  2-D non-linearity)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gotta run BG startup before coming here
nbins = 6;
scale = 2.4;
proj = input.proj(:,[1 2]);
stds = std(proj);
bins = scale*[-stds(1) -stds(2); stds(1) stds(2)];
bins = [nbins nbins; bins];
[Na,Xa] = hist2(proj,bins);
[Ns,Xs] = hist2(proj(input.Lspike==1,:),bins);
for j = 2:max(input.Lspike)
  [tmp1,Xs] = hist2(proj(input.Lspike==j,:),bins);
  Ns = Ns+j*tmp1;
end

fr = Ns./Na;
marginal1 = sum(fr);
marginal2 = sum(fr')';
tot = sum(sum(fr));
pred1 = marginal2*marginal1
pred1 = pred1 * tot/sum(sum(pred1))
pred2 = repmat(marginal1,length(marginal2),1)+repmat(marginal2,1,length(marginal1));
pred2 = pred2./length(marginal1);
correctionterm = tot/(length(marginal1)*length(marginal2));
pred2 = pred2-correctionterm;

figure
subplot(3,2,1)
imagesc(fr);
subplot(3,2,3)
imagesc(pred1);
title(sum(sum((pred1-fr).^2)))
subplot(3,2,4)
imagesc(pred1-fr);
subplot(3,2,5)
imagesc(pred2);
title(sum(sum((pred2-fr).^2)))
subplot(3,2,6)
imagesc(pred2-fr);
colormap hot

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Testing the additive vs. multiplicative 
%  models for a population of cells.
%  (Spatial domain -- single time slice)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OPTIONS.issym = 1;
OPTIONS.disp = 1;
nbins = 6;
scale = 2;
maxT = 20;
alpha = (0.0001)/6;   % For defining the RF (remember, 6 tests per pixel!)
monks = {'C','CORN';'K','KIMBALL';'T','TIGER'};
path = '/home/horwitz/Cortex/Data';
list = '/home/horwitz/Cortex/lists/bluelist';
[filelist, spikecodes] = readfilelist(list,monks);
for whichcell = 1:length(filelist)
   % First, loading files
   filename = char(filelist{whichcell});
   spikecd = str2num(spikecodes{whichcell});
   disp([filename,': ',num2str(spikecd)]);
   [foo, index, raster, genparams] = readcort2([path,'/',filename],spikecd);
   paranum = get_paradigm([path,'/',filename],1000);
   if (paranum == 60)
      indices_whtnsSpatial2;
   else
      indices_whtnsSpatial;
   end

   [E,V] = StimStats(genparams,paranum);

   RFmod('init',{E,V,alpha});
   [out] = getWhtnsStats(foo,index,raster,genparams,maxT,0,'spike',paranum,'RFmod');
   nspikes = out{1};
   startframe = out{2};
   endframe = out{3};
   nframes = endframe-startframe+1;
   RF = out{4};
   STA = out{5};
   STV = out{6};

   [i,j,STAmaxframe] = ind2sub(size(STA),find(STA == max(max(max(STA)))));
   [i,j,STVmaxframe] = ind2sub(size(STV),find(STV == max(max(max(STV)))));
   whichframe = STVmaxframe;

   feval('STCOVmod4','init',[genparams(GEN_STIMNELEMX)*genparams(GEN_STIMNELEMY) 3 whichframe whichframe]);
   [out] = getWhtnsStats3(foo,index,raster,genparams,whichframe,0,'spike',paranum,'STCOVmod4');
   STCOV = out{2};
   STA = out{1};

   % Getting rid of stuff that's "not in RF"
   inrf = [1:192];
   RFbool = permute(repmat(RF,[1 1 3]),[3 1 2]);
   RFbool = logical(RFbool(:));
   inrf = inrf(RFbool);

   STCOV = STCOV(inrf,inrf);
   STA = STA(inrf);
   STV = diag(STCOV)-mean(diag(STCOV));

   [v1,d1,flag] = eigs(STCOV,1,OPTIONS);

   basis = zeros(192,2);
   basis(inrf,1) = STA;
   basis(inrf,2) = v1;
   basis = mkbasis(basis);

   stimdurations = index(:,SIN_STIMOFF) - index(:,SIN_STIMON);
   STPROJmod('init',basis,whichframe-1,sum(stimdurations)/10,[64 3 1])
   [out] = getWhtnsStats3(foo,index,raster,genparams,whichframe,0,'spike',paranum,'STPROJmod');
   proj = out{1}; 
   Lspike = out{2};

   stds = std(proj);
   bins = scale*[-stds(1) -stds(2); stds(1) stds(2)];
   bins = [nbins nbins; bins];
   [Na,Xa] = hist2(proj,bins);
   [Ns,Xs] = hist2(proj(Lspike==1,:),bins);
   for j = 2:max(Lspike)
      [tmp1,Xs] = hist2(proj(Lspike==j,:),bins);
      Ns = Ns+j*tmp1;
   end

   fr = Ns./Na;
   marginal1 = sum(fr);
   marginal2 = sum(fr')';
   tot = sum(sum(fr));
   pred1 = marginal2*marginal1
   pred1 = pred1 * tot/sum(sum(pred1))
   pred2 = repmat(marginal1,length(marginal2),1)+repmat(marginal2,1,length(marginal1));
   pred2 = pred2./length(marginal1);
   correctionterm = tot/(length(marginal1)*length(marginal2));
   pred2 = pred2-correctionterm;

   figure
   subplot(3,2,1)
   imagesc(fr);
   title(filename);
   subplot(3,2,3)
   imagesc(pred1);
   title(sum(sum((pred1-fr).^2)))
   subplot(3,2,4)
   imagesc(pred1-fr);
   subplot(3,2,5)
   imagesc(pred2);
   title(sum(sum((pred2-fr).^2)))
   subplot(3,2,6)
   imagesc(pred2-fr);
   colormap hot
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Comparing multiplicative and additive 
%  models for the interaction between STA and
%  PC1 (in the temporal domain -- for a single
%  pixel).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nbins = 6;
scale = 2;
maxT = 20;
monks = {'C','CORN';'K','KIMBALL';'T','TIGER'};
path = '/home/horwitz/Cortex/Data';
list = '/home/horwitz/Cortex/lists/bluelist';
[filelist, spikecodes] = readfilelist(list,monks);
for whichcell = 1:length(filelist)
   % First, loading files
   filename = char(filelist{whichcell});
   spikecd = str2num(spikecodes{whichcell});
   disp([filename,': ',num2str(spikecd)]);
   [foo, index, raster, genparams] = readdata([path,'/',filename],spikecd);
   paranum = get_paradigm([path,'/',filename],1000);
   if (paranum == 60)
      indices_whtnsSpatial2;
   else
      indices_whtnsSpatial;
   end

   [E,V] = StimStats(genparams,paranum);

   [out] = getWhtnsStats(foo,index,raster,genparams,maxT,0,'spike',paranum,'STAmod');
   STA = out{1};
   STV = out{2};
   n = out{3};

   [i,j,k] = ind2sub(size(STA),find(STA == max(max(max(STA)))));

   nframes = ceil(sum(index(:,SIN_STIMOFF)-index(:,SIN_STIMON))/10);
   AllStimDim =  [nframes 1 3 maxT];
   AllStimmod('init',{i,AllStimDim});

   [out] = getWhtnsStats(foo,index,raster,genparams,maxT,0,'all',...
                      paranum,'AllStimmod');

   AllStim = out{1};
   Lspike = out{2};
   AllStim = vectorizestim(permute(AllStim,[1 2 4 3]));
   SpikeStim = AllStim(Lspike > 0,:);
   [pc,score,latent] = princomp(SpikeStim);

   STA_t = squeeze(STA(i,:,:))';
   STA_t = STA_t(:);
   basis = mkbasis([STA_t, pc(:,1)]);

   proj = AllStim*basis;  

   stds = std(proj);
   bins = scale*[-stds(1) -stds(2); stds(1) stds(2)];
   bins = [nbins nbins; bins];
   [Na,Xa] = hist2(proj,bins);
   [Ns,Xs] = hist2(proj(Lspike==1,:),bins);
   for j = 2:max(Lspike)
      [tmp1,Xs] = hist2(proj(Lspike==j,:),bins);
      Ns = Ns+j*tmp1;
   end

   fr = Ns./Na;
   marginal1 = sum(fr);
   marginal2 = sum(fr')';
   tot = sum(sum(fr));
   pred1 = marginal2*marginal1;
   pred1 = pred1 * tot/sum(sum(pred1));
   pred2 = repmat(marginal1,length(marginal2),1)+repmat(marginal2,1,length(marginal1));
   pred2 = pred2./length(marginal1);
   correctionterm = tot/(length(marginal1)*length(marginal2));
   pred2 = pred2-correctionterm;

   figure
   subplot(3,2,1)
   imagesc(fr);
   title(filename);
   subplot(3,2,3)
   imagesc(pred1);
   title(sum(sum((pred1-fr).^2)))
   subplot(3,2,4)
   imagesc(pred1-fr);
   subplot(3,2,5)
   imagesc(pred2);
   title(sum(sum((pred2-fr).^2)))
   subplot(3,2,6)
   imagesc(pred2-fr);
   colormap hot
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Looking at the timing difference (latency) between
%  STA and PCs in spatial white noise experiments.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maxT = 20;
monks = {'C','CORN';'K','KIMBALL';'T','TIGER'};
path = '/home/horwitz/Cortex/Data';
list = '/home/horwitz/Cortex/lists/BYlist2';
[filelist, spikecodes] = readfilelist(list,monks);
for whichcell = 1:length(filelist)
   % First, loading files
   filename = char(filelist{whichcell});
   spikecd = str2num(spikecodes{whichcell});
   disp([filename,': ',num2str(spikecd)]);
   [foo, index, raster, genparams] = readdata([path,'/',filename],spikecd);
   paranum = get_paradigm([path,'/',filename],1000);
   if (paranum == 60)
      indices_whtnsSpatial2;
   else
      indices_whtnsSpatial;
   end

   [E,V] = StimStats(genparams,paranum);
   [out] = getWhtnsStats(foo,index,raster,genparams,maxT,0,'spike',...
                      paranum,'STAmod');

   STA = out{1};
   minval = min(STA(:));
   maxval = max(STA(:));
   maxval = max([maxval; -minval]);
   [i,j,k] = ind2sub(size(STA),find(abs(STA)==maxval))
   STA_t = squeeze(STA(i,:,:));
   STA_t = STA_t./(sum(sum(STA_t)))

   nframes = ceil(sum(index(:,SIN_STIMOFF)-index(:,SIN_STIMON))/10);
   AllStimDim =  [nframes 1 3 maxT];
   AllStimmod('init',{i,AllStimDim});

   [out] = getWhtnsStats(foo,index,raster,genparams,maxT,0,'all',...
                      paranum,'AllStimmod');


   AllStim = out{1};
   Lspike = out{2};
   AllStim = vectorizestim(permute(AllStim,[1 2 4 3]));
   SpikeStim = AllStim(Lspike > 0,:);
   [pc,score,latent] = princomp(SpikeStim);

   PC1energy = pc(:,1).^2;
   PC1energy  = PC1energy./sum(PC1energy);
   PC1energy = sum(reshape(PC1energy,maxT,3)')
   PC1tpeak =  find(PC1energy == max(PC1energy));

   PC2energy = pc(:,2).^2;
   PC2energy  = PC2energy./sum(PC2energy);
   PC2energy = sum(reshape(PC2energy,maxT,3)')
   PC2tpeak =  find(PC2energy == max(PC2energy));

   PC12energy = (latent(1)*pc(:,1)).^2+(latent(2)*pc(:,2)).^2;
   PC12energy  = PC12energy./sum(PC12energy);
   PC12energy = sum(reshape(PC12energy,maxT,3)')
   PC12tpeak =  find(PC12energy == max(PC12energy));   

   STAenergy = sum(STA_t.^2./sum(sum(STA_t.^2)));
   STAtpeak = find(STAenergy == max(STAenergy));

   data = [data; k STAtpeak PC1tpeak PC2tpeak PC12tpeak]
end

hist(data(:,3)-data(:,2));
[h,p] = ttest(data(:,3)-data(:,2))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Looking at the spectral sensitivity of a complex cell.
%   Doing PCA on the single frame/pixel where the variance(?)
%   is the highest.  This is to test the hypothesis that 
%   Complex cells are not color opponent, but some simple
%   cells are.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

monks = {'C','CORN';'K','KIMBALL';'T','TIGER'};
datapath = '/home/horwitz/Cortex/Data';
list = '/home/horwitz/Cortex/lists/complex_swn';
[filelist, spikecodes] = readfilelist(list,monks);
data = [];
for whichcell = 1:length(filelist)
   % First, loading files
   filename = [datapath,'/',char(filelist{whichcell})];
   spikecd = str2num(spikecodes{whichcell});
   disp([filename,': ',num2str(spikecd)]);
   [foo, index, raster, genparams,e1h,e1v,e1t] = readdata(filename,spikecd);
   paradigmcd = get_paradigm(filename,1000);

   % Doing the work
   [out] = getWhtnsStats(foo,index,raster,genparams,ntimepts,0,'spike',paradigmcd,'STAmod');
   STA = out{1}; 
   STV = out{2}; 
   n = out{3};
   [whichpix,whichgun,whichtime] = ind2sub(size(STV),find(STV == max(max(max(STV)))));
   nspikes = sum(sum(~isnan(raster)));
   SpikeStimDim =  [nspikes 1 3 whichtime];
   AllStimmod('init',{whichpix,SpikeStimDim});
   [out] = getWhtnsStats(foo,index,raster,genparams,whichtime,0,'spike',...
                    paradigmcd,'AllStimmod');
   SpikeStim = out{1};
   SpikeStim = squeeze(SpikeStim(:,:,:,end));
   [pc,score,latent] = princomp(SpikeStim);
   data = [data; pc(:,1)']
end
data = data./repmat(sum(abs(data'))',1,3);
plot(data(:,2),data(:,3),'k.')
set(gca,'Xlim',[-1 1],'Ylim',[-1 1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%  Comparing various models accounting for the cell's
%  spiking response.  The models that I'm currently
%  considering include:
%     1) Mean spike rate only (as a baseline)
%     2) STA and some estimate of the static NL
%     3) STA/PC1 and some estimate of the 2-D NL
%     4) STA/Random-other-vector and some estimate of 2-D NL
%
%  The goal here is to quantify how big an effect this 
%  second excitatory axis is.
%
%  Make the stuff for finding the basis it's own
%  little routine.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nbins = 120;
scale = 2;
modeldataprop = .5;
iter = 10;
data = [];
proj = input.proj;
Lspike = input.Lspike;

% By this point that you get here you should have a variable called
% "proj" (which is an nx3 vector of stimulus projections) and 
% a vector called "Lspike" which is a nx1 vector of spike counts.

for i = 1:iter
   nstim = size(proj,1);
   idxs = randperm(nstim)';
   L = logical(idxs < round(modeldataprop*nstim));
   proj_mod = proj(L,:);
   proj_test = proj(~L,:);
   Lspike_mod = Lspike(L);
   Lspike_test = Lspike(~L);

%%%%%%%%%%%%%%%%%%
%%% SNL models %%%
%%%%%%%%%%%%%%%%%%

   %%% Fitting the model
   stds = mean(std(proj_mod));
   bininfo = scale*repmat([-stds; stds],1,3);
   bininfo = [repmat(nbins,1,3);bininfo];
   errors = [];
   for j = 1:2   % which basis: STA/PC1 or STA/Randomvector
      [Na,Xa] = hist2(proj_mod(:,[1 j+1]),bininfo(:,[1 j+1]));
      [Ns,Xs] = hist2(proj_mod(Lspike_mod==1,[1 j+1]),bininfo(:,[1 j+1]));
      for k = 2:max(Lspike_mod)
         [tmp1,Xs] = hist2(proj_mod(Lspike_mod==k,[1 j+1]),bininfo(:,[1 j+1]));
         Ns = Ns+k*tmp1;
      end
      bins = Xa{1};
      resolution = bins(2)-bins(1);
      meanmat = [bins*sum(Na/sum(sum(Na)))' bins*sum(Na'/sum(sum(Na)))'];
      X2mat =  [(bins.^2)*sum(Na/sum(sum(Na)))' (bins.^2)*sum(Na'/sum(sum(Na)))'];
      varmat = X2mat-meanmat.^2;
      hmat = 2.214*sqrt(varmat).*(length(proj_mod)^(-1/6));

      xkernel1 = [-hmat(1):resolution:hmat(1)];
      kernel1 = 3/4*(1-(xkernel1/hmat(1)).^2);
      xkernel2 = [-hmat(2):resolution:hmat(2)];
      kernel2 = 3/4*(1-(xkernel2/hmat(2)).^2);
      kernel = kernel2'*kernel1;   % Taking the outerproduct
      kernel = kernel./sum(sum(kernel));
      smoothNa = conv2(Na,kernel,'valid');  % For now using a single kernel
      smoothNs = conv2(Ns,kernel,'valid');  % For now using a single kernel
      bins1 = linspace(Xa{1}(1),Xa{1}(end),size(smoothNa,1))*size(smoothNa,1)/nbins;
      bins2 = linspace(Xa{2}(1),Xa{2}(end),size(smoothNa,2))*size(smoothNa,2)/nbins;
      fr_2D = smoothNs./smoothNa;

      fr_1D = conv(sum(kernel),sum(Ns))./conv(sum(kernel),sum(Na));
      fr_1D([1:length(kernel)-1]) = [];
      fr_1D([end-length(kernel)+2:end]) = [];
  
      binwidth = [bins1(2)-bins1(1) bins2(2)-bins2(1)];
      leftedges = [bins1(1)-binwidth(1)/2 bins2(1)-binwidth(2)/2];
      rightedges = [bins1(end)+binwidth(1)/2 bins2(end)+binwidth(2)/2];

%%% Taking a look at the set of test stimuli

%%% Getting rid of the test stimuli that fall outside the acceptance range
      L = proj_test(:,[1 j+1]) > repmat(leftedges,size(proj_test,1),1);
      L = L & proj_test(:,[1 j+1])  < repmat(rightedges,size(proj_test,1),1);
      L = logical(L(:,1)&L(:,2));
      proj_t = proj_test(L,[1 j+1]);
      Lspike_t = Lspike_test(L,:);
      nstim_t = size(proj_t,1);

% Simplest model
      Lspike_tmp = repmat(mean(Lspike_mod),nstim_t,1);
      err1 = mean((Lspike_tmp-Lspike_t).^2);

% 1-D SNL model
      whichbins = (proj_t(:,1)-leftedges(1))/...
            (rightedges(1)-leftedges(1))*(length(fr_1D));
      whichbins = ceil(whichbins);
      pred_fr = fr_1D(whichbins)';
      err2 = mean((pred_fr-Lspike_t).^2);

% 2-D SNL model (depending on what 'j' is we use different basis vectors)

      whichbins = [(proj_t(:,1)-leftedges(1))/...
            (rightedges(1)-leftedges(1))*(length(fr_1D))...
            (proj_t(:,2)-leftedges(2))/...
            (rightedges(2)-leftedges(2))*(length(fr_1D))];
      whichbins = ceil(whichbins);
% Reversal of whichbins in line below because we're switching
% from using (X,Y) coordinates to (row,column) indices.  Also,
% flipping the rows since numbering for rows starts at the top.
      pred_fr = fr_2D(sub2ind(size(fr_2D),length(bins1)-whichbins(:,2)+1,whichbins(:,1)));
      err3 = mean((pred_fr-Lspike_t).^2);
      if (j == 1)
         errors = [err1 err2 err3];
      else
         errors = [errors err3];
         data = [data; errors];
         bar(errors)
      end
   end  % of "j" loop
end

figure;
subplot(2,1,1);
hist((data(:,2)-data(:,4))./(data(:,1)-data(:,2)))
mean((data(:,2)-data(:,4))./(data(:,1)-data(:,2)))
subplot(2,1,2);
hist((data(:,1)-data(:,3))./(data(:,1)-data(:,2)))
mean((data(:,1)-data(:,3))./(data(:,1)-data(:,2)))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Looking at the Principal Components
%  decomposition of a single pixel.
%  Trying to see which eigenvalues are 
%  significant by a randomization test.
%  (Tried the Tracy-Widom way, but it gave
%  too many false positives.)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

niter = 2000;
maxT = 20;
alpha = 0.05;
monks = {'C','CORN';'K','KIMBALL';'T','TIGER'};
datapath = '/home/horwitz/Cortex/Data';
list = '/home/horwitz/Cortex/lists/BYlist2';
[filelist, spikecodes] = readfilelist(list,monks);
data = {};
for whichcell = 1:length(filelist)
   % Initializing data structure
   data{whichcell}.name = [];
   data{whichcell}.spikecd = [];
   data{whichcell}.pvalues = [];
   data{whichcell}.vectors = [];

   % Loading files
   filename = [datapath,'/',char(filelist{whichcell})];
   spikecd = str2num(spikecodes{whichcell});
   disp([filename,': ',num2str(spikecd)]);
   [foo, index, raster, genparams,e1h,e1v,e1t] = readdata(filename,spikecd);
   paradigmcd = get_paradigm(filename,1000);
   [ExpectedValues, Variances] = StimStats(genparams, paradigmcd);
   [out] = getWhtnsStats(foo,index,raster,genparams,maxT,0,'spike',paradigmcd,'STAmod');
   STA = out{1}; 
   STV = out{2}; 
   n = out{3};
   [whichpix,whichgun,whichtime] = ind2sub(size(STV),find(STV == max(max(max(STV)))));
   indices_whtnsSpatial2;
   nframes = ceil(sum(index(:,SIN_STIMOFF)-index(:,SIN_STIMON))/10);
   AllStimDim =  [nframes 1 3 maxT];
   AllStimmod('init',{whichpix AllStimDim});
   [out,infostruct] = getWhtnsStats(foo,index,raster,genparams,maxT,0,'all',paradigmcd,'AllStimmod');
   AllStim = vectorizestim(squeeze(out{1}));
   Lspike = logical(out{2});
   stdmat = repmat(sqrt(Variances),maxT,1);
   stdmat = stdmat(:)';
   AllStim = AllStim./repmat(stdmat,size(AllStim,1),1);
   tmpStim = AllStim(Lspike>0,:);
   n = sum(Lspike>0);
   p_rand = 0;
   while (p_rand <= alpha)
      xx = tmpStim'*tmpStim;	% Covariance matrix
      [v,d] = eig(xx);
      d = diag(d);
      [d,sortidxs] = sort(d);
      v = v(:,sortidxs);
   
      % Randomization test
      orig_d = max(d);
      ds = [];
      for j = 1:niter
         startidx = unidrnd(length(Lspike)-2)+1;   % 2 to end-1
         L = Lspike([startidx:end, 1:startidx-1]);
         xx = AllStim(L,:)'*AllStim(L,:);		% Covariance matrix
         [v1,d1] = eig(xx);
         ds = [ds; max(diag(d1))];
      end
      p_rand = sum(ds>=orig_d)./niter;

      % Orthogonalizing stimulus vectors
      projections = tmpStim*v(:,end);
      vmat = repmat(v(:,end)',n,1);
      tmpStim = tmpStim-vmat.*repmat(projections,1,size(tmpStim,2));
      
      % Storing the data away
      tmpv = reshape(v(:,end),3,maxT)';  % Reordering the vector elements
      data{whichcell}.name = filename;
      data{whichcell}.spikecd = spikecd;
      data{whichcell}.pvalues = [data{whichcell}.pvalues, p_rand];
      data{whichcell}.vectors = [data{whichcell}.vectors, tmpv(:)];
   end
end
nsig = [];
for i = 1:length(data)
   nsig = [nsig; sum(data{i}.pvalues <= 0.05)];
end
hist(nsig,[0:max(nsig)]);
xlabel('# signficant eigenvectors');
ylabel('count');

for i = 1:length(data)
   nvects = length(data{i}.pvalues);
   figure; 
   set(gcf,'DefaultAxesColorOrder',[1 0 0; 0 1 0; 0 0 1]);
   set(gcf,'DefaultAxesXtick',[]);
   for j = 1:nvects
      subplot(nvects,1,j);
      plot(reshape(data{i}.vectors(:,j),maxT,3))
      title(num2str(data{i}.pvalues(j)));
   end
end

% Spectral profile of significant eigenvectors
rgbs = [];
for i = 1:length(data)
   nvects = sum(data{i}.pvalues <= 0.05);
   for j = 1:nvects
      tmp = reshape(data{i}.vectors(:,j),maxT,3); 
      L = any(abs(tmp) == max(max(abs(tmp))),2);
      rgbs = [rgbs; tmp(L,:)];
   end
end
rgbs = rgbs./repmat(sum(abs(rgbs),2),1,3);
figure;
plot(rgbs(:,2),rgbs(:,3),'k.');
set(gca,'XLim',[-1 1],'Ylim',[-1 1]);
axis square; xlabel('Green'); ylabel('Blue');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Working on an analysis method to ask the
%  question, are subunits of color cells 
%  matched or not?
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maxT = 20;
name = '/home/horwitz/Cortex/Data/DONDERS/D001.1';
[foo, index, raster, genparams,e1h,e1v,e1t] = readdata(name,1);
paranum = get_paradigm(name,1000);
[ExpectedValues, Variances] = StimStats(genparams, paranum);
[out] = getWhtnsStats(foo,index,raster,genparams,maxT,0,'spike',paranum,'STAmod');
STA = out{1}; 
STV = out{2}; 
n = out{3};
[whichpix,whichgun,whichtime] = ind2sub(size(STA),find(STA == max(max(max(STA)))));
peakguns = STA(whichpix,:,whichtime);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tmpSTA = STA(:,:,whichtime);
[u,s,v] = svd(tmpSTA);
pred = u(:,1)*s(1,1)*v(:,1)';
weights = u(:,1);
data = [s(1,1) max(weights) min(weights)];
for i = 1:2000
   tmpSTA = slotMachine(STA(:,:,whichtime)')';
   [u,s,v] = svd(tmpSTA);
   pred = u(:,1)*s(1,1)*v(:,1)';
   weights = s(1,1)*u(:,1);
   data = [data; s(1,1) max(weights) min(weights)];
end
p1 = sum(data([2:end],1)>=data(1,1))./(size(data,1)-1)
p2 = sum(data([2:end],2)>=data(1,2))./(size(data,1)-1)
p3 = sum(data([2:end],3)<=data(1,3))./(size(data,1)-1)

%%%%%%% I don't really like this approach so much.  I can't tell
%%%%%%% the difference between only one subunit and two subunits
%%%%%%% with different color tuning.  Plus it's really convoluted.

% For significance testing
tmat = zeros(size(STA,1),size(STA,2));
for gun = 1:size(STA,2)
   im = STA(:,gun,whichtime);
   im = im./(sqrt(Variances(gun)/n));
   tmat(:,gun) = im;
end
chi2mat = sum(tmat.^2,2)

% Finding the template
[whichpix,whichgun,whichtime] = ind2sub(size(STA),find(STA == max(max(max(STA)))));
tmpSTA = STA(:,:,whichtime);
peakguns = STA(whichpix,:,whichtime);
peakguns = peakguns./norm(peakguns);
proj = peakguns*tmpSTA';
off = find(proj == min(proj));
data = proj(off);

for i = 1:2000
   tmpSTA = slotMachine(STA(:,:,whichtime)')';
   [whichpix,whichgun] = ind2sub(size(tmpSTA),find(tmpSTA == max(max(max(tmpSTA)))));
   peakguns = tmpSTA(whichpix,:);
   peakguns = peakguns./norm(peakguns);
   proj = peakguns*tmpSTA';
   data = [data; min(proj)];
end
p1 = sum(data([2:end],1)<=data(1,1))./(size(data,1)-1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First I need a method for picking out
% cells with subunits.
% Then I can worry about how to compare colors.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Several ways of looking at the RF.
% The idea here is to come up with some
% sort of way of saying that the luminance
% gain control stuff is confined to the RF
% that you get by spike-triggered averaging.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

maxT = 20;
name = '/home/horwitz/Cortex/Data/KIMBALL/K025.4c5';
[foo, index, raster, genparams,e1h,e1v,e1t] = readcort2 (name,102);
paradigmcd = get_paradigm(name,1000);
if (paradigmcd == 60)
   indices_whtnsSpatial2;
else
   indices_whtnsSpatial;
end

[ExpectedValues,Variances] = StimStats(genparams,paradigmcd);
[out] = getWhtnsStats(foo,index,raster,genparams,maxT,0,'spike',paradigmcd,'STAmod');
STA = out{1}; 
STV = out{2}; 
n = out{3};

s = size(STA);
zmat1 = zeros(s);   % STA
zmat2 = zeros(s);   % STV
zmat3 = zeros(s);   % PC1

for gun = 1:size(STA,2)
   for timept = 1:size(STA,3)
      % STA
      im = STA(:,gun,timept);
      im = im./(sqrt(Variances(gun)/n));
      zmat1(:,gun,timept) = im;
      % STV
      im = STV(:,gun,timept);
      im = (n-1)*im/Variances(gun);
      zmat2(:,gun,timept) = (im-(n-1))/sqrt(2*(n-1));
   end
end

% Now doing PC1
% It seems to be the case the whitening everything up front causes
% a systematic decrease in the eigenvalue so that we keep falling
% off the left tail of the Tracy-Widom distribution.  This might be
% because once we do that we're no longer sampling from a huge 
% population, but rather from a relatively modest one with 
% sigma = I (exactly).

ndims = 3*maxT;
stdmat = repmat(sqrt(Variances),maxT,1);
stdmat = stdmat(:)';
zmat3 = zeros(s(1),1);
for whichpix = 1:64
   nframes = ceil(sum(index(:,SIN_STIMOFF)-index(:,SIN_STIMON))/10);
   AllStimDim =  [nframes 1 3 maxT];
   AllStimmod('init',{whichpix AllStimDim});
   [out,infostruct] = getWhtnsStats(foo,index,raster,genparams,maxT,0,'all',paradigmcd,'AllStimmod');
   AllStim = vectorizestim(permute(out{1},[1 2 4 3]));
   SpikeStim = [];
   Lspike = out{2};
   SpikeStim = AllStim(Lspike > 0,:);
%   AllStimMean = mean(AllStim);
%   c = cov(AllStim);
%   [u,s,v] = svd(c);
%   wz = v*sqrt(s);
%   wz = inv(wz);
%   SpikeStim1 = SpikeStim*wz'-repmat(AllStimMean,size(SpikeStim,1),1);
%   SpikeStim1 = SpikeStim*wz';
%   [v,d] = eig(SpikeStim1'*SpikeStim1);
   SpikeStim1 = SpikeStim./repmat(stdmat,size(SpikeStim,1),1);
   [v,d] = eig(SpikeStim1'*SpikeStim1);
   d = max(diag(d));

   ndims = 3*maxT;
   n = size(SpikeStim,1);
   mu_eig = (sqrt(n-1)+sqrt(ndims))^2;
   sigma_eig = (sqrt(n-1)+sqrt(ndims))*((1/sqrt(n-1)+1/sqrt(ndims))^(1/3));
   zeig = (d-mu_eig)/sigma_eig;

   zeig = (zeig+1.21)/1.27;
   zmat3(whichpix) = zeig
end
figure;
subplot(3,1,1);
[whichpixel, whichgun, whichframe] = ind2sub(size(zmat1),...
                       find(zmat1 == max(max(max(zmat1)))));
imagesc(reshape(squeeze(zmat1(:,whichgun,whichframe)>4),8,8))
axis image; colormap(gray);
subplot(3,1,2);
[whichpixel, whichgun, whichframe] = ind2sub(size(zmat2),...
                       find(zmat2 == max(max(max(zmat2)))));
imagesc(reshape(squeeze(zmat2(:,whichgun,whichframe)>4),8,8))
axis image; colormap(gray);
subplot(3,1,3);
imagesc(reshape(zmat3>4,8,8))
axis image; colormap(gray);

tmp1 = squeeze(zmat2(:,whichgun,whichframe)>4);
tmp2 = zmat3(:)>4;
tmp3 = squeeze(zmat2(:,whichgun,whichframe))./norm(squeeze(zmat2(:,whichgun,whichframe)));
tmp4 = zmat3(:)./norm(zmat3(:));
corrcoef(tmp3,tmp4)
sum(tmp1&tmp2)./sum(tmp2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Looking at the single-frame (i.e. spatial)
%  1st PC for only pixels "in the RF" (by
%  some criterion or other).  The idea here
%  is to document the fact that the luminance
%  gain control is spatial contrast too, not 
%  just temporal contrast.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OPTIONS.issym = 1;
OPTIONS.disp = 1;
maxT = 20;
alpha = 0.00001/6;   % For defining the RF (remember, 6 tests per pixel!)
filename = '/home/horwitz/Cortex/Data/KIMBALL/K121.3';
spikecd = 102;
[foo, index, raster, genparams] = readdata(filename,spikecd);
paranum = get_paradigm(filename,1000);
if (paranum == 60)
  indices_whtnsSpatial2;
else
   indices_whtnsSpatial;
end

[E,V] = StimStats(genparams,paranum);

RFmod('init',{E,V,alpha});
[out] = getWhtnsStats(foo,index,raster,genparams,maxT,0,'spike',paranum,'RFmod');
nspikes = out{1};
startframe = out{2};
endframe = out{3};
nframes = endframe-startframe+1;
RF = out{4};
STA = out{5};
STV = out{6};

[i,j,STAmaxframe] = ind2sub(size(STA),find(STA == max(max(max(STA)))));
[i,j,STVmaxframe] = ind2sub(size(STV),find(STV == max(max(max(STV)))));
[STAmaxframe STVmaxframe]
whichframe = STVmaxframe;

feval('STCOVmod4','init',[genparams(GEN_STIMNELEMX)*genparams(GEN_STIMNELEMY) 3 whichframe whichframe]);
[out] = getWhtnsStats3(foo,index,raster,genparams,whichframe,0,'spike',paranum,'STCOVmod4');
STCOV = out{2};
STA = out{1};

% Getting rid of stuff that's "not in RF"
inrf = [1:192];
RFbool = permute(repmat(RF,[1 1 3]),[3 1 2]);
RFbool = logical(RFbool(:));
inrf = inrf(RFbool);

STCOV = STCOV(inrf,inrf);
STA = STA(inrf);
STV = diag(STCOV)-mean(diag(STCOV));

[v1,d1,flag] = eigs(STCOV,1,OPTIONS);
STAmat = zeros(64,3);
STAmat(find(RF),:) = reshape(STA,3,sum(sum(RF)))';
fact = max(max(abs(STAmat)));
STAmat = (STAmat+fact)./(2*fact)

PC1mat = zeros(64,3);
PC1mat(find(RF),:) = reshape(v1,3,sum(sum(RF)))';
fact = max(max(abs(PC1mat)));
PC1mat = (PC1mat+fact)./(2*fact)

figure;
subplot(2,1,1); 
image(reshape(STAmat,8,8,3))
axis image;
subplot(2,1,2);
image(reshape(PC1mat,8,8,3))
axis image;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Looking at a single frame PC1, all 
%  all pixel, but each RGB projected 
%  onto a luminance axis before PCA.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LUMrgb = [6.8 25.1 9.9]; 
LUMrgb = LUMrgb./norm(LUMrgb');
OPTIONS.issym = 1;
OPTIONS.disp = 1;
maxT = 20;
filename = '/home/horwitz/Cortex/Data/KIMBALL/K121.3';
spikecd = 102;
[foo, index, raster, genparams] = readdata(filename,spikecd);
paranum = get_paradigm(filename,1000);
if (paranum == 60)
  indices_whtnsSpatial2;
else
   indices_whtnsSpatial;
end

[out] = getWhtnsStats(foo,index,raster,genparams,20,0,'spike',...
                      paranum,'STAmod');

STA = out{1};
STV = out{2};

[i,j,STAmaxframe] = ind2sub(size(STA),find(STA == max(max(max(STA)))));
[i,j,STVmaxframe] = ind2sub(size(STV),find(STV == max(max(max(STV)))));
[STAmaxframe STVmaxframe]
whichframe = STVmaxframe;

feval('STCOVmod4','init',[genparams(GEN_STIMNELEMX)*genparams(GEN_STIMNELEMY) 3 whichframe whichframe]);
[out] = getWhtnsStats3(foo,index,raster,genparams,whichframe,0,'spike',paranum,'STCOVmod4');
STCOV = out{2};
c = [LUMrgb(1); zeros(191,1)];
r = [LUMrgb, zeros(1,189)];
w = toeplitz(c,r);
w = w([1:3:192],:);
STCOV = w*STCOV*w';

[v1,d1,flag] = eigs(STCOV,1,OPTIONS);
imagesc(reshape(v1,8,8))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Trying to improve our estimates
%  of time/color principal components.
%  For each pixel in the RF I'm going
%  calculate a covariance matrix, then
%  I'm going to average these covariance
%  matrices together and look at the
%  eigenvectors.  This actually worked
%  pretty well.  Might want to make it a 
%  standard practice.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maxT = 20;
alpha = 0.000001/6;
filepath = '/home/horwitz/Cortex/Data';
monks = {'C','CORN';'K','KIMBALL';'T','TIGER'};
list = '/home/horwitz/Cortex/lists/BYlist2';
[filelist, spikecodes] = readfilelist(list,monks);
for whichcell = 1:length(filelist)
   filename = char(filelist{whichcell});
   spikecd = str2num(spikecodes{whichcell});
   disp([filename,': ',num2str(spikecd)]);
   [foo, index, raster, genparams] = readdata([filepath,'/',filename],spikecd);
   paranum = get_paradigm([filepath,'/',filename],1000);
   if (paranum == 60)
     indices_whtnsSpatial2;
   else
     indices_whtnsSpatial;
   end

   [E,V] = StimStats(genparams,paranum);

   RFmod('init',{E,V,alpha});
   [out] = getWhtnsStats(foo,index,raster,genparams,maxT,0,'spike',paranum,'RFmod');
   nspikes = out{1};
   RF = out{4};
   
   MasterCovMat = [];
   pix = find(RF)
   for i = pix';
      i
      nspike = sum(sum(~isnan(raster)));
      SpikeStimDim =  [nspike 1 3 maxT];
      AllStimmod('init',{i,SpikeStimDim});

      [out] = getWhtnsStats(foo,index,raster,genparams,20,0,'spike',...
                      paranum,'AllStimmod');

      SpikeStim = vectorizestim(permute(out{1},[1 2 4 3]));
      MasterCovMat(:,:,find(i==pix)) = cov(SpikeStim);
   end
   L = logical(ones(length(pix),1));  % For selecting a subset of pixels
   [v,d] = eig(mean(MasterCovMat(:,:,L),3));
   d = diag(d);
   [y,i] = sort(d);
   v = v(:,i);
   maxd = max(d);
   maxds = maxd;
   grandv = v(:,end);

   figure;
   for i = 1:length(pix)
      side = ceil(sqrt(length(pix)+1));
      subplot(side, side, i);
      [v,d] = eig(MasterCovMat(:,:,i));
      d = diag(d);
      [y,j] = sort(d);
      v = v(:,j);
      maxd = max(d);
      plot(reshape(v(:,end),20,3));
      title(pix(i));
      maxds = [maxds; max(d)];
   end
   subplot(side, side, i+1);
   plot(reshape(grandv,20,3));
   title('Grand Mean');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Trying out a new definition of 
%  color-opponency: if some confidence
%  ellipsoid around the STA in LMS space
%  doesn't contain any points where all
%  cone weights have the same sign, then
%  the cell is color-opponent.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maxT = 20;
monks = {'C','CORN';'K','KIMBALL';'T','TIGER'};
filepath = '/home/horwitz/Cortex/Data';
list = '/home/horwitz/Cortex/lists/simple_swn';
[filelist, spikecodes] = readfilelist(list,monks);
data = {};
for whichcell = 1:length(filelist)
   filename = char(filelist{whichcell});
   spikecd = str2num(spikecodes{whichcell});
   disp([filename,': ',num2str(spikecd)]);
   [foo, index, raster, genparams] = readdata([filepath,'/',filename],spikecd);
   paranum = get_paradigm([filepath,'/',filename],1000);
   if (paranum == 60)
      indices_whtnsSpatial2;
   else
      indices_whtnsSpatial;
   end
   monserial = genparams(GEN_MONSERIAL1)+...
               genparams(GEN_MONSERIAL2)*2^8+...
               genparams(GEN_MONSERIAL3)*2^16;
   cal = getGlobalParams('cal',monserial);
   fundamentals = getGlobalParams('fundamentals');
   M = fundamentals*cal.P_device;
   [out] = getWhtnsStats(foo,index,raster,genparams,20,0,'spike',...
                         paranum,'STAmod');

   STA = out{1};
   STV = out{2};
   n = out{3};
   [whichpix,whichgun,whichtime] = pickapix(STA,STV,n,1);
   STArgb = squeeze(STA(whichpix,:,whichtime))';
   STVrgb = squeeze(STV(whichpix,:,whichtime))';
   COVrgb = diag(STVrgb)/n;
   STAlms = inv(M)'*STArgb;
   COVlms = inv(M)'*COVrgb*inv(M);
   points = [];
   for i = 1:iter
       vect = mkbasis(normrnd(0,1,3,1));  
       point = STAlms+1.96*sqrtm(COVlms)*vect;
       points = [points; point'];
   end

   data = [data; sum(sign(points))./iter]
end

%% Definition of cone-opponency: At least one cone 
%% has to be going consistently up and another has
%% to be going consistently down.

L1 = any(transpose(data == 1));
L2 = any(transpose(data == -1));
acceptlist = L1&L2;

%% Here's the breakdown of cone opponency...
%% There are 12 categories of cells: 6 for cells for which
%% the signs of all of the cones is consistent plus an additional
%% 6 for cells that have consistency of only two cone types.

tmpdata = data(acceptlist,:);
tmpdata(abs(tmpdata) < 1) = 0
design = fullfact([2 2 2]);
design(design == 2) = -1;
design([1,end],:) = []
design = [design; scroll([1 -1 0; -1 1 0],[0 0])]
design = [design; scroll([1 -1 0; -1 1 0],[1 0])]
design = [design; scroll([1 -1 0; -1 1 0],[2 0])]
counts = zeros(size(design,1),1);
for i = 1:size(counts,1)
   mask = repmat(design(i,:),size(tmpdata,1),1);
   counts(i) = sum(all(transpose(mask == tmpdata)));
end
bar(counts)
%% Lots of L vs. M.
%% M goes with S often (just like Conway 2001).  Gain control artifact? 
%% Doesn't seem to be since I get the same result when I just look
%% at the simple cells.  Cone fundamental related error?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Trying a method for seeing whether two 
%  neighboring pixels might have the same
%  color/time STA.  Basically doing a some 
%  sort of multivariate t-test on the two STAs.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

maxT = 20;
filename = '/home/horwitz/Cortex/Data/CORN/C127.2c3';
spikecd = 102;
[foo, index, raster, genparams] = readdata(filename,spikecd);
paranum = get_paradigm(filename,1000);
if (paranum == 60)
  indices_whtnsSpatial2;
else
   indices_whtnsSpatial;
end

[E,V] = StimStats(genparams,paranum);

[out] = getWhtnsStats(foo,index,raster,genparams,maxT,0,'spike',paranum,'STAmod');
STA = out{1};
STV = out{2};
nspikes = out{3};

[whichpix,whichgun,whichtime] = pickapix(STA,STV,nspikes,1,E,V);
SpikeStimDim =  [nspikes 1 3 maxT];
AllStimmod('init',{whichpix,SpikeStimDim});
[out] = getWhtnsStats(foo,index,raster,genparams,maxT,0,'spike',...
                    paranum,'AllStimmod');
SpikeStim = out{1};
SpikeStim = vectorizestim(permute(SpikeStim,[1 2 4 3]));
STAcenter = mean(SpikeStim);
STCOVcenter = cov(SpikeStim);

% Now we have the STA and STCOV for the hottest pixel.  Now we're
% going to look at a bunch of the surrounding pixels to see whether
% their STAs might be the same as the one in the middle.

[Icenter,Jcenter] = ind2sub([8 8],whichpix);

data = [Icenter Jcenter 0];
STAs = STAcenter;
for i = Icenter+[-1 0 1]
   for j = Jcenter+[-1 0 1]
      if ((i > 0) & (i < 9) & (j > 0) & (j < 9) & ~(i == Icenter & j ==Jcenter))
         whichpix = sub2ind([8 8],i,j);
         AllStimmod('init',{whichpix,SpikeStimDim});
         [out] = getWhtnsStats(foo,index,raster,genparams,maxT,0,'spike',...
                            paranum,'AllStimmod');
         SpikeStim = out{1};
         SpikeStim = vectorizestim(permute(SpikeStim,[1 2 4 3]));
         STAalt = mean(SpikeStim);
         STCOValt = cov(SpikeStim);

         Sp = (STCOVcenter+STCOValt)./2;
         T2 = (STAcenter-STAalt) * Sp * (STAcenter-STAalt)';
         n = 2*size(SpikeStim,1);
         p = size(SpikeStim,2);
         F = T2 * (n-p-1)/p*(n-2);
         pval = fcdf(F,p,n-p-1)
         STAs = [STAs; STAalt];
         data = [data; i j pval];
      end
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  As above, but this time we're just
%  looking at the angle between the 
%  two vectors (the correlation coefficient).
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maxT = 20;
filename = '/home/horwitz/Cortex/Data/CORN/C111.4c5';
spikecd = 102;
[foo, index, raster, genparams] = readdata(filename,spikecd);
paranum = get_paradigm(filename,1000);
if (paranum == 60)
  indices_whtnsSpatial2;
else
   indices_whtnsSpatial;
end

[E,V] = StimStats(genparams,paranum);

[out] = getWhtnsStats(foo,index,raster,genparams,maxT,0,'spike',paranum,'STAmod');
STA = out{1};
STV = out{2};
nspikes = out{3};

[whichpix,whichgun,whichtime] = pickapix(STA,STV,nspikes,1,E,V);
%[Icenter,Jcenter] = ind2sub([8 8],whichpix);
STAcenter = squeeze(STA(whichpix,:,:))';
STAcenter = STAcenter(:);

STAmat = vectorizestim(permute(STA,[1 3,2]))';
r = corrcoef(STAmat);
r = r(whichpix,:);

t = r.*sqrt(size(STAmat,1)-2)./sqrt(1-r.^2)
p = 1-tcdf(t,size(STAmat,1)-2)
goodpixels = find(p<0.001);
figure;
subplot(3,1,1);
plot(STAmat(:,goodpixels))
subplot(3,1,2);
plot(mean(STAmat(:,goodpixels)'))
subplot(3,1,3);
plot(fliplr(reshape(mean(STAmat(:,goodpixels)'),20,3)))

[out, goodpixels] = extractColorTimeSTA(STA, 0.0001);
figure; axes; set(gcf,'DefaultAxesColorOrder',[1 0 0; 0 1 0; 0 0 1]);
plot(reshape(squeeze(mean(STA(goodpixels,:,:)))',20,3))

%  The assumptions behind the correlation hypothesis test
%  is totally busted.  This is a flawed analysis.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Testing for color/time separability.
%
%   The idea is that we're going to do an SVD
%   on the color/time STA and see whether the 
%   separable prediction is within a confidence
%   ellipsoid around the STA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

maxT = 20;
filename = '/home/horwitz/Cortex/Data/CORN/C147.3';
spikecd = 102;
[foo, index, raster, genparams] = readdata(filename,spikecd);
paranum = get_paradigm(filename,1000);
if (paranum == 60)
  indices_whtnsSpatial2;
else
   indices_whtnsSpatial;
end

[out] = getWhtnsStats(foo,index,raster,genparams,maxT,0,'spike',paranum,'STAmod');
STA = out{1};
STV = out{2};
nspikes = out{3};
[whichpix,whichgun,whichframe] = pickapix(STA,STV,nspikes,1);

meanmat = squeeze(STA(whichpix,:,:));
meanvect = reshape(meanmat',1,3*maxT);
[u,s,v] = svd(meanmat);
predmat = u(:,1)*s(1,1)*v(:,1)';
predvect = reshape(predmat',1,3*maxT);
residvect = meanvect-predvect;
residmat = reshape(residvect,maxT,3)';

nspikes = sum(sum(~isnan(raster)));
SpikeStimDim =  [nspikes 1 3 maxT];
AllStimmod('init',{whichpix,SpikeStimDim});
[out] = getWhtnsStats(foo,index,raster,genparams,maxT,0,'spike',...
                 paranum,'AllStimmod');
SpikeStim = out{1};
SpikeStim = squeeze(vectorizestim(permute(SpikeStim,[1 2 4 3])));
nstim = size(SpikeStim,1);
ndims = size(SpikeStim,2);
S = cov(SpikeStim)./nstim;   % Covariance of the mean
dist = (meanvect-predvect)*inv(S)*(meanvect-predvect)'
stat = dist*(nstim-ndims)/((nstim-1)*ndims);
p = 1-fcdf(stat,ndims,nstim-ndims);

%%% Plotting stuff
figure; set(gcf,'DefaultAxesColorOrder',[1 0 0; 0 1 0; 0 0 1]);
subplot(3,1,1);
plot(fliplr(meanmat)');
subplot(3,1,2);
plot(fliplr(predmat)');
subplot(3,1,3);
plot(fliplr(residmat)');

%%% Same thing in spatial domain

feval('STCOVmod4','init',[genparams(GEN_STIMNELEMX)*genparams(GEN_STIMNELEMY) 3 whichframe whichframe]);
[out] = getWhtnsStats3(foo,index,raster,genparams,whichframe,0,'spike',paranum,'STCOVmod4');

% Next three lines only if we've done STCOVmod
meanvect = out{1}';
STCOV = out{2}/nstim;
meanmat = reshape(meanvect,3,genparams(GEN_STIMNELEMX)*genparams(GEN_STIMNELEMY));

% Come here if we haven't done STCOVmod
meanmat = squeeze(STA(:,:,whichframe))';
meanvect = meanmat(:)';

ndims = 3*genparams(GEN_STIMNELEMX)*genparams(GEN_STIMNELEMY);
[u,s,v] = svd(meanmat);
predmat = u(:,1)*s(1,1)*v(:,1)';
predvect = reshape(predmat,1,ndims);
residvect = meanvect-predvect;
residmat = reshape(residvect,3,genparams(GEN_STIMNELEMX)*genparams(GEN_STIMNELEMY));

dist = (meanvect-predvect)*inv(STCOV)*(meanvect-predvect)'
stat = dist*(nstim-ndims)/((nstim-1)*ndims);
p = 1-fcdf(stat,ndims,nstim-ndims);

%%%% Plotting stuff

normresidvect = residvect*0.5/max(abs(residvect))+0.5;
normresidmat = reshape(normresidvect,3,genparams(GEN_STIMNELEMX),genparams(GEN_STIMNELEMY));

normpredvect = predvect*0.5/max(abs(predvect))+0.5;
normpredmat = reshape(normpredvect,3,genparams(GEN_STIMNELEMX),genparams(GEN_STIMNELEMY));

normmeanvect = meanvect*0.5/max(abs(meanvect))+0.5;
normmeanmat = reshape(normmeanvect,3,genparams(GEN_STIMNELEMX),genparams(GEN_STIMNELEMY));

figure;
subplot(3,1,1); image(permute(normmeanmat,[2,3,1]))
axis square; set(gca,'Visible','off')
subplot(3,1,2); image(permute(normpredmat,[2,3,1]))
axis square; set(gca,'Visible','off')
subplot(3,1,3); image(permute(normresidmat,[2,3,1]))
axis square; set(gca,'Visible','off')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Loking at orientation tuning and space/time
%  separability.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maxT = 20;
filename = '/home/horwitz/Cortex/Data/CORN/C101.2'
spikecd = 102;
[foo, index, raster, genparams] = readdata(filename,spikecd);
paranum = get_paradigm(filename,1000);
if (paranum == 60)
  indices_whtnsSpatial2;
else
   indices_whtnsSpatial;
end

[out] = getWhtnsStats(foo,index,raster,genparams,maxT,0,'spike',...
                      paranum,'STAmod');

STA = out{1};
STV = out{2};
nspikes = out{3};

[whichpix,whichgun,whichframe] = pickapix(STA,STV,nspikes,1);
STA_s = squeeze(STA(:,:,whichframe));
[u,s,v] = svd(STA_s);
spacepredvect = u(:,1);
rgbpredvect = v(:,1);
predmat = reshape(spacepredvect,genparams(GEN_STIMNELEMX),genparams(GEN_STIMNELEMY));

fftim = fftshift(abs(fft2(predmat./norm(predmat(:)))));
p_tuning = max(max(fftim))./norm(fftim(:));
fftim = fftim([1:5],:)
% Preferred spatial frequencies below 
% Shifted so 0 is DC, 1 is 1 cycle, 2 is 2 cycles, etc.
peakfft = max(max(fftim));
p_omega_h = max(find(max(fftim) == peakfft))-5
p_omega_v = 5-find(max(fftim' == max(max(fftim'))))

indices_whtnsSpatial;
arrayHeight = genparams(GEN_STIMH)/10;
arrayWidth = genparams(GEN_STIMW)/10;

p_sf = sqrt(p_omega_v^2+p_omega_h^2)/arrayWidth  % In cycles/degree
p_orient = atan2(p_omega_v, p_omega_h)*180/pi+90 % In degrees

for j = 1:2   % X or Y
   STA_t = [];
   STA_tmp = reshape(STA,genparams(GEN_STIMNELEMX),genparams(GEN_STIMNELEMY),3,maxT);
   [whichpix_row, whichpix_col] = ind2sub([genparams(GEN_STIMNELEMX),genparams(GEN_STIMNELEMY)],whichpix)
   if (j == 1)
      STA_tmp = squeeze(STA_tmp(whichpix_row,:,:,:))
   else
      STA_tmp = squeeze(STA_tmp(:,whichpix_col,:,:))
   end
   for frame = 1:maxT
      STA_t(:,frame) = STA_tmp(:,:,frame)*rgbpredvect;
   end
   [u,s,v] = svd(STA_t);
   stats(j) = s(1,1)/trace(s);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Testing for double opponency.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha = 0.005;
maxT = 20;
filename = '/home/horwitz/Cortex/Data/CORN/C152.3'
spikecd = 102;
[foo, index, raster, genparams] = readdata(filename,spikecd);
paranum = get_paradigm(filename,1000);
if (paranum == 60)
  indices_whtnsSpatial2;
else
   indices_whtnsSpatial;
end

[ExpectedValues, Variances] = StimStats(genparams, paranum);

[out] = getWhtnsStats(foo,index,raster,genparams,maxT,0,'spike',...
                      paranum,'STAmod');

STA = out{1};
STV = out{2};
nspikes = out{3};
[whichpix,whichgun,whichframe] = pickapix(STA,STV,nspikes,1);
peakrgb = squeeze(STA(whichpix,:,whichframe));
rgbs = STA(:,:,whichframe);
rgbprojections = rgbs*(peakrgb/norm(peakrgb))'
projvar = Variances*(peakrgb.^2/norm(peakrgb.^2))'
se = sqrt(projvar/nspikes);
sigpix = logical(abs(rgbprojections/se) > norminv(1-(alpha/2)));
sigpixsigns = sign(rgbprojections(sigpix));
[i,j] = ind2sub([8 8],find(sigpix));
[i_peak,j_peak] = ind2sub([8 8],whichpix);

pixpos = sqrt((i - i_peak).^2 + (j - j_peak).^2);
excPixPos = pixpos(pixpos > 0 & sigpixsigns == 1);
inhPixPos = pixpos(sigpixsigns == -1);
meanExcDist = mean(excPixPos);
stdExcDist = std(excPixPos);
nExc = length(excPixPos);
meanInhDist = mean(pixpos(sigpixsigns == -1));
stdInhDist = std(pixpos(sigpixsigns == -1));
nInh = length(inhPixPos);
[meanExcDist stdExcDist nExc]
[meanInhDist stdInhDist nInh]

% Prob on observed nInh based on expected based on Type I error rate
p = binopdf(nInh,63,alpha/2)+(1-binocdf(nInh,63,alpha/2))
% Probability we see this many or more

% Making a plot
template = zeros(8,8);
template(sigpix == 1 & sign(rgbprojections) == 1) = 1;
template(sigpix == 1 & sign(rgbprojections) == -1) = -1;
imagesc(template)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Trying Duane Nykamp's method of sorting
%  simple from complex cells.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maxT = 12;
filename = '/home/horwitz/Cortex/Data/CORN/C101.2'; spikecd = 102;

filename = '/home/horwitz/Cortex/Data/KIMBALL/K069.3'; spikecd = 102;

[foo, index, raster, genparams] = readdata(filename,spikecd);
paranum = get_paradigm(filename,1000);
if (paranum == 60)
  indices_whtnsSpatial2;
else
   indices_whtnsSpatial;
end

[ExpectedValues, Variances] = StimStats(genparams, paranum);
RFmod('init',{ExpectedValues,Variances,0.000001/6});
[out] = getWhtnsStats(foo,index,raster,genparams,maxT,0,'spike',...
                      paranum,'RFmod');

nspikes = out{1};
startframe = out{2};
endframe = out{3};
RFmap = out{4};
STA = out{5};
STV = out{6};

%%% Stripping out the stuff that's not in the RF
STA(:,:,[1:startframe-1,endframe+1:end]) = [];
STV(:,:,[1:startframe-1,endframe+1:end]) = [];
STA(find(~RFmap),:,:) = [];
STV(find(~RFmap),:,:) = [];
npix = sum(RFmap(:));
nSTAframes = endframe-startframe+1;

nframes = foo(:,I_NFRAMES1)+foo(:,I_NFRAMES2).*2^8;
nframes = max(nframes - maxT + 1, 0);
nframes = sum(nframes(:));
ER = nspikes/nframes;

STS = STA.*nspikes;
STSnorm = STS./sqrt(repmat(Variances,[npix 1 nSTAframes]));
EX = STSnorm./nframes;

STS2 = (STV*(nspikes*(nspikes-1))+STS.^2)/nspikes;
EX2 = STS2./(nframes.*repmat(Variances,[npix 1 nSTAframes]));

% EX is sum(x)/nframes -- average x where x ~ N(0,1)
% EX2 is sum(x^2)/nframes -- average squared x where x ~ N(0,1)

g = ((EX2./ER) - (EX./ER).^2).^(-1);
f = (EX/ER).*g;

v1 = f;
v2 = (1-g)/sqrt(2);

v1norm = sum(f(:).^2);
v2norm = sum((1-g(:)).^2);
alpha = v2norm/(v1norm+v2norm)

greg = permute(v2,[1 3 2]);
greg = reshape(greg,size(greg,1)*size(greg,2),size(greg,3));
[u,s,v] = svd(greg);
rgb = v(:,1);

%%%%  A little simulation
nframes = 10000;
ndims = 64*3*12;
x = normrnd(0,1,nframes,ndims);
s = binornd(1,.5,nframes,1);
xx = x(logical(s),:);
nspikes = sum(s);

EX2 = sum(xx.^2)/nframes;
EX = sum(xx)/nframes;
ER = nspikes/nframes;

g = ((EX2./ER) - (EX./ER).^2).^(-1);
f = (EX/ER).*g;

v1 = sum(f(:).^2);
v2 = sum((1-g(:)).^2);
alpha = v2/(v1+v2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Trying some stuff to get the PC1 in space/color
%  after projecting the data onto the subspace orthogonal
%  to the STA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maxT = 20;
filename = '/home/horwitz/Cortex/Data/KIMBALL/K069.3'; spikecd = 102;
[foo, index, raster, genparams] = readdata(filename,spikecd);
paranum = get_paradigm(filename,1000);
if (paranum == 60)
  indices_whtnsSpatial2;
else
   indices_whtnsSpatial;
end
[out] = getWhtnsStats(foo,index,raster,genparams,maxT,0,'spike',paranum,'STAmod');
STA = out{1};
STV = out{2};
n = out{3};
[ExpectedValues, Variances] = StimStats(genparams, paranum);
[whichpix,whichgun,whichframe] = pickapix(STA,STV,n,2,ExpectedValues,Variances);
feval('STCOVmod4','init',[genparams(GEN_STIMNELEMX)*genparams(GEN_STIMNELEMY) 3 whichframe whichframe]);
[out] = getWhtnsStats3(foo,index,raster,genparams,whichframe,0,'spike',paranum,'STCOVmod4');
STAv = out{1};
STCOV = out{2};
P = STAv*inv(STAv'*STAv)*STAv';
plot(P*normrnd(0,1,192,1),'m-');   % Cool, P*anything is a scaled version of the STA
K = (eye(size(P))-P);   % K projects orthogonal to STA -- this is trhe matrix we want

OPTIONS.issym = 1;
OPTIONS.disp = 1;
[v1,d1,flag] = eigs(STCOV,1,'LM',OPTIONS);
[v2,d1,flag] = eigs(K*STCOV*K',1,'LM',OPTIONS);
v1'*STAv
v2'*STAv
plotvariables = {'STAv','v1','v2'}
for j = 1:length(plotvariables)
   tmp = eval(char(plotvariables(j)));
   tmp = tmp*0.5/max(max(abs(tmp)))+0.5;
   tmp = permute(reshape(tmp,3,8,8),[2 3 1]);
   % Choosing sign of PC1 to be codirectional with STA
   cormat = corrcoef([tmp(:) STAv(:)]);
   if (cormat(1,2) < 0)
      tmp = -(tmp-.5)+.5;
   end
   subplot(3,1,j);
   image(tmp);
   axis ij;
   axis image;
   set(gca,'Visible','off');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   Looking at the signal to noise ratio
%   of the STA and PC1 and it's relationship
%   to alpha (the Nykamp parameter).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OPTIONS.issym = 1;
OPTIONS.disp = 1;
niter = 20;
maxT = 12;
filename = '/home/horwitz/Cortex/Data/CORN/C101.2'; spikecd = 102;

filename = '/home/horwitz/Cortex/Data/KIMBALL/K069.3'; spikecd = 102;

[foo, index, raster, genparams] = readdata(filename,spikecd);
paranum = get_paradigm(filename,1000);
if (paranum == 60)
  indices_whtnsSpatial2;
else
   indices_whtnsSpatial;
end

[ExpectedValues, Variances] = StimStats(genparams, paranum);
RFmod('init',{ExpectedValues,Variances,0.000001/6});
[out] = getWhtnsStats(foo,index,raster,genparams,maxT,0,'spike',...
                      paranum,'RFmod');

nspikes = out{1};
startframe = out{2};
endframe = out{3};
RFmap = out{4};
STA = out{5};
STV = out{6};

nframes = foo(:,I_NFRAMES1)+foo(:,I_NFRAMES2).*2^8;
nframes = max(nframes - maxT + 1, 0);
nframes = sum(nframes(:));
[whichpix, whichgun, whichframe] = pickapix(STA,STV,nspikes,6,Variances)
out{2} = whichframe;
out{3} = whichframe;
[alpha, v1, v2] = NykampAlpha(out, Variances, nframes)

data = [];
for i = 1:niter+1
   if (i == 1)
      opt = 'spike';
   else
      opt = 'spikerand';
   end
   feval('STCOVmod4','init',[genparams(GEN_STIMNELEMX)*genparams(GEN_STIMNELEMY) 3 whichframe whichframe]);
   [out] = getWhtnsStats3(foo,index,raster,genparams,whichframe,0,opt,paranum,'STCOVmod4');
   STA = out{1};
   STCOV = out{2};
   [v,d,flag] = eigs(STCOV,1,'LM',OPTIONS);
   data = [data; norm(STA), d(1)]
end
(data(1,1)-mean(data([2:end],1)))./std(data([2:end],1))
(data(1,2)-mean(data([2:end],2)))./std(data([2:end],2))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Trying a few nonlinear transformations
%   of the gun values to improve the contrast
%   of figures (STA, PC1, etc.)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

name = '/home/horwitz/Cortex/Data/DONDERS/D018.1'
paranum = get_paradigm(name,1000);
[foo, index, raster, genparams,e1h,e1v,e1t] = readdata (name,1);

[out] = getWhtnsStats(foo,index,raster,genparams,20,0,'spike',...
                      paranum,'STAmod');

STA = out{1};
STV = out{2};
maxval = max(max(max(abs(STA))));
alpha = 1;
normfactor = 0.5/(1/(1+exp(-alpha/2))-0.5);
figure;
for i = 1:size(STA,3)
   im = squeeze(STA(:,:,i));
   im = reshape(im,sqrt(size(im,1)),sqrt(size(im,1)),size(im,2));
   subplot(ceil(sqrt(size(STA,3))),ceil(sqrt(size(STA,3))),i)
   imboost(im, alpha, maxval);
end

figure; axes; hold on;
tmp = [-0.5:.01:0.5];
y = 1./(1+exp(-alpha*(tmp)))-.5;

plot(tmp,y)
plot(tmp,normfactor*y,'m-')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Mixing STA and PCs to maximize color/
%  space separability.
%
%  Doesn't work well because of noise in
%  the PCs
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

niter = 10000;
name = '/home/horwitz/Cortex/Data/TIGER/T006.3'
spikecd = 102;
paranum = get_paradigm(name,1000);
[foo, index, raster, genparams,e1h,e1v,e1t] = readdata (name,spikecd);

[out] = getWhtnsStats(foo,index,raster,genparams,20,0,'spike',...
                      paranum,'STAmod');
if (paranum == 60)
   indices_whtnsSpatial2;
else
   indices_whtnsSpatial;
end

STA = out{1};
STV = out{2};
n = out{3};
[whichpix, whichgun, whichframe] = pickapix(STA, STV, n, 1);
STAmat = squeeze(STA(:,:,whichframe));
feval('STCOVmod4','init',[genparams(GEN_STIMNELEMX)*genparams(GEN_STIMNELEMY) 3 whichframe whichframe]);
[out] = getWhtnsStats3(foo,index,raster,genparams,whichframe,0,'spike',paranum,'STCOVmod4');
STCOV = out{2};
[v,d] = eigs(STCOV,3,'LM')
basis = mkbasis([STAmat(:) v(:,[1 2])]);
mixes = unifrnd(0,1,size(basis,2)-1,niter);
mixes = [mixes; 1-sum(mixes)];
data = zeros(1,niter);
for i = 1:size(mixes,2)
   vect = basis*mixes(:,i);
   vmat = reshape(vect,3,genparams(GEN_STIMNELEMX)*genparams(GEN_STIMNELEMY));
   [u,s,v] = svd(vmat);
   data(i) = s(1)/trace(s);
end
max(data)
idx = find(data == max(data));
mixes(:,idx)
v = basis*mixes(:,idx);
plot(reshape(v,3,64)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Creating some sort of Monte Carlo test in the
%%  spatial domain.  This still takes way too long.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

name = '/home/horwitz/Cortex/Data/DONDERS/D034.2c6'
spikecd = 1;
paranum = get_paradigm(name,1000);
[foo, index, raster, genparams,e1h,e1v,e1t] = readdata (name,spikecd);
if (paranum == 60)
   indices_whtnsSpatial2;
else
   indices_whtnsSpatial;
end

whichframe = 7;
nframes = ceil(sum(index(:,SIN_STIMOFF)-index(:,SIN_STIMON))/10);
AllStimDim =  [nframes 64 3 1];
AllStimmod2('init',{whichframe,AllStimDim});

[out] = getWhtnsStats(foo,index,raster,genparams,20,0,'spike',...
                      paranum,'AllStimmod2');

SpikeStim = out{1};
Lspike = out{2};

trueS = cov(vectorizestim(SpikeStim));
[v,d] = eig(trueS);
[d,i] = sort(diag(d));
orig_d = max(d);

%% Randomization test part
niter = 2000;
ds = [];
s = size(SpikeStim);
for i = 1:niter
   for j = 1:s(3)
      j
      idx = reshape(randperm(s(1)*s(2)),[s(1) s(2)]);
      tmp = squeeze(SpikeStim(:,:,j));
      SpikeStim(:,:,j) = tmp(idx);
   end
   S = cov(vectorizestim(SpikeStim));
   [v,d] = eig(S);
   ds = [ds; max(diag(d))];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Trying it again, but just creating random
% covariance matrices with the appropriate
% mean and variance.  This still takes a pretty
% long time...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

name = '/home/horwitz/Cortex/Data/DONDERS/D034.2c6'
spikecd = 1;
paranum = get_paradigm(name,1000);
[foo, index, raster, genparams,e1h,e1v,e1t] = readdata (name,spikecd);
if (paranum == 60)
   indices_whtnsSpatial2;
else
   indices_whtnsSpatial;
end
[ExpectedValues, Variances] = StimStats(genparams, paranum);

whichframe = 7;
nframes = ceil(sum(index(:,SIN_STIMOFF)-index(:,SIN_STIMON))/10);
AllStimDim =  [nframes 64 3 1];
AllStimmod2('init',{whichframe,AllStimDim});

[out] = getWhtnsStats(foo,index,raster,genparams,20,0,'spike',...
                      paranum,'AllStimmod2');

SpikeStim = out{1};
Lspike = out{2};

trueS = cov(vectorizestim(SpikeStim));
[v,d] = eig(trueS);
[d,i] = sort(diag(d));
orig_d = max(d);

%% Randomization test part
niter = 2000;
ds = [];
s = size(SpikeStim);
for i = 1:niter
   i
   for j = 1:3
      SpikeStim(:,:,j) = normrnd(0,sqrt(Variances(j)),s(1),s(2));
   end
   S = cov(vectorizestim(SpikeStim));
   [v,d] = eig(S);
   ds = [ds; max(diag(d))];
end

p = sum(ds > orig_d)/niter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%  Figuring out the spatial frequency content 
%  of the flickering checkerboard stimulus.
%  Power spectral density of the stimulus.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pixperelem = 15;
nelem = 10;
pixperdeg = pixperelem/0.22;
a = zeros(nelem*pixperelem,nelem*pixperelem);
a([2*pixperelem:1:3*pixperelem],[2*pixperelem:1:3*pixperelem]) = 1;
figure;
imagesc(a);
axis image;
ticks = get(gca,'XTick');
set(gca,'XTickLabel',ticks.*0.22/pixperelem);
set(gca,'YTickLabel',ticks.*0.22/pixperelem);
colormap(gray);

b = xcorr2(a);
figure;
imagesc(b);
axis image;
midpoint = (size(b,1)+1)/2;
colormap(gray);

maxdist = midpoint/pixperdeg;   % dist from 0 to edge in deg
linecoefs = [(size(b,1)-midpoint)/maxdist midpoint]

ticks = [-1:1];
set(gca,'XTick',ticks*linecoefs(1)+linecoefs(2))
set(gca,'XTickLabel',ticks);
set(gca,'YTick',ticks*linecoefs(1)+linecoefs(2))
set(gca,'YTickLabel',ticks);
degperwidth = size(b,1)/pixperdeg;  % size of edge in degrees

c = fftshift(abs(fft2(b)));
figure;
imagesc(c);
axis image;
colormap(gray)

maxfreq = (size(b,1)/2)/degperwidth  % in c/d
%freqs = round(linspace(-maxfreq, maxfreq,9))
freqs = [-30:10:30];
linecoefs = [(size(c,1)-midpoint)/maxfreq midpoint]
set(gca,'XTick',freqs*linecoefs(1)+linecoefs(2))
set(gca,'XTickLabel',freqs);
set(gca,'YTick',freqs*linecoefs(1)+linecoefs(2))
set(gca,'YTickLabel',freqs);

idxs = round([-2 2].*linecoefs(1)+linecoefs(2))
p1 = sum(sum(c([idxs(1):idxs(2)],[idxs(1):idxs(2)])));
ptot = sum(sum(c));

p1/ptot

%% Approximately 60% of the total power of this
%% stimulus is below 2 cycles/degree

idxs = round([-5 5].*linecoefs(1)+linecoefs(2))
p1 = sum(sum(c([idxs(1):idxs(2)],[idxs(1):idxs(2)])));
p1/ptot


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Distribution of spectra in my experiment
%  Going out 2 standard deviations for each gun.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = 380:5:780;
STDEV = 0.09;
cal = getGlobalParams('cal');
weights = FindModelWeights(cal.P_ambient,cal.P_device)';
gam = cal.gammaTable+repmat(weights,size(cal.gammaTable,1),1);
bkgndRGB = getGlobalParams('bkgndrgb');
bkgndrgb = [gam(bkgndRGB(1)+1,1) gam(bkgndRGB(2)+1,2) gam(bkgndRGB(3)+1,3)];
device = cal.P_device*1000;
mean = bkgndrgb*device';
bound1 = ((bkgndrgb+2*STDEV)*device');
bound2 = ((bkgndrgb-2*STDEV)*device');
figure; axes; hold on;
h = patch([x, fliplr(x)],[bound1,fliplr(bound2)],[.5 .5 .5]);
plot(x,mean,'k-','LineWidth',2);
set(gca,'XLim',[380 780]);
plot(x,(bkgndrgb+2*STDEV./sqrt(10))*device','k--','LineWidth',2);
plot(x,(bkgndrgb-2*STDEV./sqrt(10))*device','k--','LineWidth',2);

figure; axes; hold on;
plot(x,(bound1-mean)*1000,'k-','LineWidth',2);
set(gca,'XLim',[380 780], 'YLim',[0 3]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Getting cone weights for a 
%  few example cells.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

name = '/home/horwitz/Cortex/Data/KIMBALL/K092.3'
spikecd = 102;
paranum = get_paradigm(name,1000);
[foo, index, raster, genparams,e1h,e1v,e1t] = readdata (name,spikecd);
if (paranum == 60)
   indices_whtnsSpatial2;
else
   indices_whtnsSpatial;
end

monserial = genparams(GEN_MONSERIAL1)+...
            genparams(GEN_MONSERIAL2)*2^8+...
            genparams(GEN_MONSERIAL3)*2^16;
cal = getGlobalParams('cal',monserial);
fundamentals = getGlobalParams('fundamentals');
M = fundamentals*cal.P_device;

[ExpectedValues, Variances] = StimStats(genparams, paranum);
maxT = 20;
[out] = getWhtnsStats(foo,index,raster,genparams,maxT,0,'spike',...
                      paranum,'STAmod');
STA = out{1};
STV = out{2};
n = out{3};

im = squeeze(STA(:,:,5));
im = reshape(im,sqrt(size(im,1)),sqrt(size(im,1)),size(im,2));
im = (im-min(min(min(im))))/(max(max(max(im)))-min(min(min(im))));
image(im);

whichpix = 36;
% Cell C101.2: 20, 27, 43
% Cell K092.3: 27, 28, 35
STA_t = squeeze(STA(whichpix,:,:));

plot(STA_t')
[u,s,v] = svd(STA_t');
lms = inv(M')*v(:,1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Looking at STAs/PC1s near and far from 
%  the occurance of fixational saccades.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SACWIN = 250;
name = '/home/horwitz/Cortex/Data/KIMBALL/K050.4c6'
spikecd = 102;
paranum = get_paradigm(name,1000);
[foo, index, raster, genparams,e1h,e1v,e1t] = readdata (name,spikecd);
if (paranum == 60)
   indices_whtnsSpatial2;
else
   indices_whtnsSpatial;
end
samplerate = get_msperADsample(name);
sacstats = getSacData(index, e1h, e1v, e1t, samplerate);
raster_sac = zeros(size(raster));
raster_nosac = zeros(size(raster));
for i = 1:size(foo,1)
   Ltrial = sacstats.trialnums == i;
   starttimes = sacstats.starttimes(Ltrial);
   starttimes = starttimes(starttimes > index(i,SIN_STIMON) &...
                           starttimes < index(i,SIN_STIMOFF));
   Lsac = logical(zeros(size(raster,2),1));
   for j = 1:size(starttimes)
      tmp = raster(i,:) - starttimes(j);
      Lsac1 = tmp > 0 & tmp < SACWIN;
      Lsac2 = tmp < 0 & tmp > -SACWIN;
   end
   raster_sac(i,:) = [raster(i,Lsac1), zeros(1,size(raster,2)-sum(Lsac1))];
   raster_nosac(i,:) = [raster(i,Lsac2), zeros(1,size(raster,2)-sum(Lsac2))];
end
[out] = getWhtnsStats(foo,index,raster_sac,genparams,20,0,'spike',...
                      paranum,'STAmod');
STA_sac = out{1};
STV_sac = out{2};
n_sac = out{3};

[out] = getWhtnsStats(foo,index,raster_nosac,genparams,20,0,'spike',...
                      paranum,'STAmod');

STA_nosac = out{1};
STV_nosac = out{2};
n_nosac = out{3};

STA = STV_nosac;
minval = min(STA(:));
maxval = max(STA(:));
minval = min([minval -maxval]);
maxval = max([maxval -minval]);
figure;
for i = 1:size(STA,3)
   im = squeeze(STA(:,:,i));
   im = reshape(im,sqrt(size(im,1)),sqrt(size(im,1)),size(im,2));
   im = (im-minval)/(maxval-minval);
   subplot(ceil(sqrt(size(STA,3))),ceil(sqrt(size(STA,3))),i)
   image(im);
   colormap(hot);
   axis('square')
   set(gca,'XTickLabel',[],'YTickLabel',[]);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Trying to get standard errors
%  of Mahalanobis distances estimates.
%  This should be pretty straightforward
%  for STAs, but more complicated 
%  for PC1s.  Gotta get the covariance
%  matrix of the characteristic [r,g,b]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

niter = 200;
maxT = 20;
name = '/home/horwitz/Cortex/Data/KIMBALL/K050.4c6'
spikecd = 102;
paranum = get_paradigm(name,1000);
[foo, index, raster, genparams,e1h,e1v,e1t] = readdata (name,spikecd);
if (paranum == 60)
   indices_whtnsSpatial2;
else
   indices_whtnsSpatial;
end
[out] = getWhtnsStats(foo,index,raster,genparams,maxT,0,'spike',...
                      paranum,'STAmod');

STA = out{1};
STV = out{2};
n = out{3};
[whichpix, whichgun, whichtime] = pickapix(STA,STV,n,1);

nspikes = sum(sum(~isnan(raster)));
AllStimDim =  [nspikes 1 3 maxT];
AllStimmod('init',{whichpix,AllStimDim});

[out] = getWhtnsStats(foo,index,raster,genparams,maxT,0,'spike',...
                      paranum,'AllStimmod');

SpikeStim = vectorizestim(out{1});
Lspike = out{2};
nstim = size(SpikeStim,1);

STAv = reshape(mean(SpikeStim),3,20)';
[u,s,v] = svd(STAv);
STArgb = v(:,1);

S = cov(SpikeStim);
[pc, latent] = eigs(S,1);
[u,s,v] = svd(reshape(pc,3,20)');
PC1rgb = v(:,1);

data = [];
for i = 1:niter
   i
   L = unidrnd(nstim,nstim,1);
   STAv = reshape(mean(SpikeStim(L,:)),3,20)';
   [u1,s1,v1] = svd(STAv);

   S = cov(SpikeStim(L,:));
   [pc, latent] = eigs(S,1);
   [u2,s2,v2] = svd(reshape(pc,3,20)');

   data = [data; v1(:,1)' v2(:,1)'];
end
S_STA = cov(data(:,[1 2 3]));
S_PC1 = cov(data(:,[4 5 6]));

figure;
plot3(data(:,4),data(:,5),data(:,6),'k.')
set(gca,'XLim',[-1 1],'YLim',[-1 1],'ZLim',[-1 1]);

LUMrgb = [6.8 25.1 9.9]'; 
LUMrgb = LUMrgb./norm(LUMrgb);

seSTArgb = std(data(:,[1 2 3])*LUMrgb);
sePC1rgb = std(data(:,[4 5 6])*LUMrgb);

dist1 = (LUMrgb-STArgb)'*inv(S_STA)*(LUMrgb-STArgb);
dist2 = (LUMrgb-PC1rgb)'*inv(S_PC1)*(LUMrgb-PC1rgb);

out = [LUMrgb'*STArgb, seSTArgb, dist1, LUMrgb'*PC1rgb, sePC1rgb, dist2];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Trying to get error estimates (error volumes)
%  for normalized cone weights.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load /home/horwitz/matlab/PopAnal/Data/LumCorBoot.mat
Lumdata = data;
load /home/horwitz/matlab/PopAnal/Data/Nykamp.mat
Nykampdata = data;
i = 5;
M = Nykampdata{i}.M;
STArgb = Lumdata{i}.STArgb;
S = Lumdata{i}.S_STA;
STAlms = inv(M')*STArgb;
C = inv(M')*S*inv(M);
normfactor = sum(abs(STAlms));
STAlms = STAlms./normfactor;
C = C./(normfactor^2)
[v,d] = eig(C);
[eigenvalues,i] = sort(diag(d));
v = v(:,i);    % reordering eigenvectors
scalefactors = sqrt(eigenvalues);
newd = diag(scalefactors);
newv = v*newd;
points = [];
for i = 1:50;
   vect = mkbasis(normrnd(0,1,3,1));
   points = [points; (newv*vect)'];
end
figure; axes; hold on;
plot3(STAlms(1),STAlms(2),STAlms(3),'k.');
plot3(points,points,points,'y.');

%%% This approach didn't work so well (for one thing,
%%% the covariance matrix in RGB is nearly singular
%%% and is even more singular in LMS.  We need to 
%%% get a bootstrap estimate of the normalized LMS 
%%% cone weight estimate error.

% First, loading files
fundamentals = getGlobalParams('fundamentals');
maxT = 20;
path = '/home/horwitz/Cortex/Data';
filename = 'CORN/C111.4c5';
spikecd = 102;
disp([filename,': ',num2str(spikecd)]);
paranum = get_paradigm([path,'/',filename],1000);
[foo, index, raster, genparams] = readdata([path,'/',filename],spikecd);
if (paranum == 60)
   indices_whtnsSpatial2;
else
   indices_whtnsSpatial;
end

monserial = genparams(GEN_MONSERIAL1)+...
            genparams(GEN_MONSERIAL2)*2^8+...
            genparams(GEN_MONSERIAL3)*2^16;
cal = getGlobalParams('cal',monserial);
M = fundamentals*cal.P_device;

[out] = getWhtnsStats(foo,index,raster,genparams,maxT,0,'spike',...
                      paranum,'STAmod');
STA = out{1};
STV = out{2};
n = out{3};
[whichpix,whichgun,whichtime] = pickapix(STA,STV,n,1);

nspikes = sum(sum(~isnan(raster)));
AllStimDim =  [nspikes 1 3 maxT];
AllStimmod('init',{whichpix,AllStimDim});

[out] = getWhtnsStats(foo,index,raster,genparams,20,0,'spike',...
                      paranum,'AllStimmod');

SpikeStim = out{1};
Lspike = out{2};
SpikeStim = vectorizestim(permute(SpikeStim,[1 2 4 3]));

niter = 200;
lmsdata = [];
for i = 1:niter+1
   if (i == 1)
      samplevect = [1:length(Lspike)];
   else
      samplevect = unidrnd(length(Lspike),length(Lspike),1);
   end
   STAv = reshape(mean(SpikeStim(samplevect,:)),20,3)';
   [u,s,v] = svd(STAv);
   STArgb = u(:,1);
   peak = find(max(sum(abs(STAv))) == sum(abs(STAv)));
   if (STAv(:,peak)'*STArgb < 0)  % Undoing a potential sign flip
      STArgb = -STArgb;
   end
   lms = inv(M')*STArgb;
   lmsdata = [lmsdata; lms'./sum(abs(lms))];
end

figure; axes; hold on;
plot3(lmsdata(:,1),lmsdata(:,2),lmsdata(:,3),'k.');
set(gca,'XLim',[-1 1],'YLim',[-1 1],'ZLim',[-1 1]);
axis square;
xlabel('L'); ylabel('M'); zlabel('S');

%%% Plotting 1 SE ellipsoid
%%% Skip ahead for spline stuff...
lin = linspace(0,2*pi,1000);
sine = [cos(lin);sin(lin)];
S = cov(lmsdata);
muhat = lmsdata(1,:);

for axidx = 1:3
   greg = zeros(3,length(lin));
   if (axidx == 1)
      idxs = [1 2]
   elseif (axidx == 2)
      idxs = [2 3]
   else
      idxs = [1 3]
   end
   [v,d] = eig(S(idxs,idxs));
   [eigenvalues,i] = sort(diag(d));
   v = v(:,i);    % reordering eigenvectors
   newd = diag(sqrt(eigenvalues))
   greg(idxs,:) = v*newd*sine;
   plot3(greg(1,:)+muhat(1),greg(2,:)+muhat(2),greg(3,:)+muhat(3),'r-');
end
plot3(muhat(1),muhat(2),muhat(3),'k*')


%%%% Trying to show confidence splines...
lmsdata = lmsdata-repmat(lmsdata(1,:),size(lmsdata,1),1);
x = lmsdata(:,1);
y = lmsdata(:,2);
binwidth = pi/20;
kernel = normpdf([-3:.2:3],0,1);
bins = [-pi+binwidth/2:binwidth:pi-binwidth/2];
[t,r] = cart2pol(x,y);
plot(t,r,'k.')
v = [];
for i = 1:length(bins)
   v(i) = sum(r(t>bins(i)-binwidth & t<bins(i)+binwidth));
end
v = conv(repmat(v,1,3),kernel);
v = v((length(v)-length(bins))/2:(length(v)+length(bins))/2)
w = linspace(-pi,pi,length(v)-1);
w = [w,w(1)];
[xs,ys] = pol2cart(w,v);
plot(xs,ys)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Searching for mechanism directions by: 
%  1) Taking the STA and PC1 and creating color/time
%     seperable versions.
%  2) Project onto these and assess seperability 
%     by chi-square
%  3) rotate, independently, the STA and PC1 away from 
%     the acute angle and return to step 1.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
name = '/home/horwitz/Cortex/Data/CORN/C175.3c4'
maxT = 20;
spikecd = 102;
paranum = get_paradigm(name,1000);
[foo, index, raster, genparams,e1h,e1v,e1t] = readdata (name,spikecd);
if (paranum == 60)
   indices_whtnsSpatial2;
else
   indices_whtnsSpatial;
end
[ExpectedValues, Variances] = StimStats(genparams, paranum);
[out] = getWhtnsStats(foo,index,raster,genparams,maxT,0,'spike',...
                      paranum,'STAmod');
STA = out{1};
STV = out{2};
n = out{3};
[whichpix,whichgun,whichtime] = pickapix(STA,STV,n,6,Variances)

%%% Getting stimuli

nframes = ceil(sum(index(:,SIN_STIMOFF)-index(:,SIN_STIMON))/10);
AllStimDim =  [nframes 1 3 maxT];
AllStimmod('init',{whichpix,AllStimDim});
[out] = getWhtnsStats(foo,index,raster,genparams,maxT,0,'all',...
                      paranum,'AllStimmod');
AllStim = out{1};
Lspike = out{2};
AllStim = vectorizestim(permute(AllStim,[1 2 4 3]));
SpikeStim = AllStim(Lspike > 0,:);
STA = mean(SpikeStim);
[pc,score,latent] = princomp(SpikeStim);
if (sign(STA*pc(:,1)) > 0)
   basis = mkbasis([STA', pc(:,1)]);
else
   basis = mkbasis([STA', -pc(:,1)]);
end

nbins = 8;
nangles = 20;
maxangle = pi/6;
scale = 1.5;
std = min(sqrt(Variances));
angles = linspace(0,maxangle,nangles);
data = zeros(length(angles),length(angles));
for theta1 = angles
   for theta2 = angles
      w1hat = cos(theta1)*basis(:,1) -...
              sin(theta1)*basis(:,2);
      w2hat = cos(theta2)*basis(:,2) -...
              sin(theta2)*basis(:,1);
      %%% Doing SVD to create separable basis functions %%%
      [u,s,v] = svd(reshape(w1hat,maxT,3));
      tmp = u(:,1)*s(1,1)*v(:,1)';
      w1hat = tmp(:)/norm(tmp(:));

      [u,s,v] = svd(reshape(w2hat,maxT,3));
      tmp = u(:,1)*s(1,1)*v(:,1)';
      w2hat = tmp(:)/norm(tmp(:));

      %%% Fitting the response surfaces

      tmpbasis = mkbasis([w1hat, w2hat]);
      proj = AllStim*tmpbasis;

      bins = scale*[-std -std; std std];
      bins = [nbins nbins; bins];
      [Na,Xa] = hist2(proj,bins);
      [Ns,Xs] = hist2(proj(Lspike==1,:),bins);
      for j = 2:max(Lspike)
         [tmp1,Xs] = hist2(proj(Lspike==j,:),bins);
         Ns = Ns+j*tmp1;
      end
      fr = Ns./Na;
      marginal1 = sum(fr);
      marginal2 = sum(fr')';
      tot = sum(sum(fr));
      pred1 = (marginal2*marginal1)./tot;
      Nshat = pred1.*Na;
      err = sum(sum((Ns-Nshat).^2./Nshat));
      data(find(angles == theta2),find(angles == theta1)) = err
   end
end
%% X dimension (constant row) is w1
%% Y dimension (constant column) is w2

figure; axes;
imagesc(data);
axis square;
set(gca,'XTickLabel',angles,'YTickLabel',angles);
xlabel('W1 angle');
ylabel('W2 angle');

tmp = (data == min(min(data)));
theta1 = angles(logical(sum(tmp)));
theta2 = angles(logical(sum(tmp')));
w1hat = cos(theta1)*basis(:,1) -...
        sin(theta1)*basis(:,2);
w2hat = cos(theta2)*basis(:,2) -...
        sin(theta2)*basis(:,1);
figure; set(gcf,'Defaultaxescolororder',[1 0 0; 0 1 0; 0 0 1]);
subplot(2,2,1);
plot(reshape(basis(:,1),maxT,3));
subplot(2,2,2);
plot(reshape(w1hat,maxT,3));
subplot(2,2,3);
plot(reshape(basis(:,2),maxT,3));
subplot(2,2,4);
plot(reshape(w2hat,maxT,3));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Taking the stimuli, converting them to a multivariate 
%   uniform density, and doing cluster analysis on the
%   spike-triggered stimuli.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

name = '/home/horwitz/Cortex/Data/CORN/C096.2c6'
name = '/home/horwitz/Cortex/Data/CORN/C111.4c5'
maxT = 20;
spikecd = 102;
paranum = get_paradigm(name,1000);
[foo, index, raster, genparams,e1h,e1v,e1t] = readdata (name,spikecd);
if (paranum == 60)
   indices_whtnsSpatial2;
else
   indices_whtnsSpatial;
end

[ExpectedValues, Variances] = StimStats(genparams, paranum);
[out] = getWhtnsStats(foo,index,raster,genparams,maxT,0,'spike',...
                      paranum,'STAmod');
STA = out{1};
STV = out{2};
n = out{3};
[whichpix,whichgun,whichtime] = pickapix(STA,STV,n,6,Variances);
nspikes = sum(sum(~isnan(raster)));
SpikeStimDim =  [nspikes 1 3 maxT];
AllStimmod('init',{whichpix,SpikeStimDim});
[out] = getWhtnsStats4(foo,index,raster,genparams,maxT,0,'spike',...
                      paranum,'AllStimmod');
SpikeStim = out{1};
SpikeStim = vectorizestim(permute(SpikeStim,[1 2 4 3]));

for i = 1:4
   [Classes,Centers,FinalDistance]=dcKMeans(SpikeStim,i);
   figure; set(gcf,'DefaultAxesColorOrder',[1 0 0; 0 1 0; 0 0 1]);
   for j = 1:size(Centers,1)
      subplot(2,2,j)
      plot(reshape(Centers(j,:),maxT,3))
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%  Calculations for the complex
%  cell paper.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load Nykamp.mat;
M = data{5}.M;
rgb = data{5}.STArgb;
lms = inv(M')*rgb;
lms./sum(abs(lms))
alpha = data{5}.alpha

M = data{83}.M
rgb = data{83}.PC1rgb;
lms = inv(M')*rgb;
lms./sum(abs(lms))
alpha = data{83}.alpha

M = data{96}.M
rgb = data{96}.STArgb;
lms = inv(M')*rgb;
lms./sum(abs(lms))
alpha = data{96}.alpha

name = '/home/horwitz/Cortex/Data/KIMBALL/K069.3'
paranum = get_paradigm(name,1000);
[foo, index, raster, genparams,e1h,e1v,e1t] = readdata (name,102);
indices_whtnsSpatial2;
[nspikes,dt] = countspikes(raster, index, SIN_STIMON, 0, SIN_STIMOFF, 0);
mean(nspikes./dt)*1000
[nspikes,dt] = countspikes(raster, index, SIN_FIX, 0, SIN_STIMON, 0);
mean(nspikes./dt)*1000

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Trying hard to find something interesting
%  in the "non-classical" surround.  Basically,
%  I'm going to find the classical RF, pool
%  across those pixels that are not in the 
%  the classical RF, and look at the PCs.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maxT = 20;
alpha = 0.0001;
name = '/home/horwitz/Cortex/Data/CORN/C129.2c6';
[foo, index, raster, genparams,e1h,e1v,e1t] = readdata(name,102);
paranum = get_paradigm(name,1000);
[ExpectedValues, Variances] = StimStats(genparams, paranum);
if (paranum == 60)
   indices_whtnsSpatial2;
else
   indices_whtnsSpatial;
end

RFmod('init',{ExpectedValues,Variances,alpha});
[out] = getWhtnsStats(foo,index,raster,genparams,maxT,0,'spike',paranum,'RFmod');

n = out{1};
framestart = out{2};
framestop = out{3};
RFmap = out{4};
STA = out{5};
nspikes = sum(sum(~isnan(raster)));

STCOV = zeros(3*maxT, 3*maxT);
Varvect = repmat(Variances,maxT,1);
Varvect = Varvect(:)';
for i = find(~RFmap)'
   i
   SpikeStimDim = [nspikes 1 3 maxT];
   AllStimmod('init',{i,SpikeStimDim});
   [out] = getWhtnsStats(foo,index,raster,genparams,maxT,0,'spike',...
                    paranum,'AllStimmod');
   SpikeStim = out{1};
   SpikeStim = vectorizestim(permute(SpikeStim,[1 2 4 3]));
   SpikeStim = SpikeStim./repmat(Varvect,size(SpikeStim,1),1);
   STCOV = STCOV + cov(SpikeStim)
end
[v,d] = eig(STCOV);
[d,sortidxs] = sort(diag(d));
v = v(:,sortidxs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Another stab at trying to find something 
%  interesting in the "non-classical" surround:
%  Basically, I'm going to find the classical
%  RF, pool across those pixels that are not
%  in the the classical RF, take the central
%  pixel, pool this data into one big old
%  vector (per SpikeStim) and do PCA.
%  The basic question is: are there any correlations
%  between what's going on in the center and 
%  the surround when the cell fires a spike.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maxT = 20;
alpha = 0.0001;
name = '/home/horwitz/Cortex/Data/CORN/C129.2c6';
[foo, index, raster, genparams,e1h,e1v,e1t] = readdata(name,102);
paranum = get_paradigm(name,1000);
[ExpectedValues, Variances] = StimStats(genparams, paranum);
Varvect = repmat(Variances,maxT,1);
Varvect = Varvect(:)';
if (paranum == 60)
   indices_whtnsSpatial2;
else
   indices_whtnsSpatial;
end

RFmod('init',{ExpectedValues,Variances,alpha});
[out] = getWhtnsStats(foo,index,raster,genparams,maxT,0,'spike',paranum,'RFmod');

n = out{1};
framestart = out{2};
framestop = out{3};
RFmap = out{4};
STA = out{5};
nspikes = sum(sum(~isnan(raster)));

[whichpix,whichgun,whichtime] = pickapix(STA,[],n,1);
AllStimmod('init',{whichpix,[nspikes, 1 3 maxT]});
[out] = getWhtnsStats(foo,index,raster,genparams,maxT,0,'spike',...
                    paranum,'AllStimmod');
SpikeStim = out{1};
SpikeStim = vectorizestim(permute(SpikeStim,[1 2 4 3]));
CenterSpikeStim = SpikeStim./repmat(Varvect,size(SpikeStim,1),1);

SurAvg = [];
for i = find(~RFmap)'
   i
   SpikeStimDim = [nspikes 1 3 maxT];
   AllStimmod('init',{i,SpikeStimDim});
   [out] = getWhtnsStats(foo,index,raster,genparams,maxT,0,'spike',...
                    paranum,'AllStimmod');
   SpikeStim = out{1};
   SpikeStim = vectorizestim(permute(SpikeStim,[1 2 4 3]));
   SpikeStim = SpikeStim./repmat(Varvect,size(SpikeStim,1),1);
   if (isempty(SurAvg))
      SurAvg = SpikeStim;
   else
      SurAvg = SurAvg + SpikeStim;
   end
end
SurAvg = SurAvg./sqrt(sum(sum(~RFmap)));
S = cov([projortho(CenterSpikeStim,mean(CenterSpikeStim)), SurAvg]);
[v,d] = eig(S);
[d,sortidxs] = sort(diag(d));
v = v(:,sortidxs);
v1 = v(:,end);
figure
subplot(3,1,1); hold on;
plot(v1([1:20]),'b-'); plot(v1([61:80]),'m-')
subplot(3,1,2); hold on;
plot(v1([21:40]),'b-'); plot(v1([81:100]),'m-')
subplot(3,1,3); hold on;
plot(v1([41:60]),'b-'); plot(v1([101:120]),'m-')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Looking at the effect of eye movements on
%  the Nykamp QIN index.
%
%  First, we're just going to take a cell and
%  randomly shift the spike-triggered stimuli
%  by different amounts.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

maxT = 20;
niter = 10;
alpha = 0.000001/6;
sigmas = [0:0.02:0.2];    % standard deviation of eye position in degrees
name = '/home/horwitz/Cortex/Data/KIMBALL/K092.3'; spikecd = 102;
name = '/home/horwitz/Cortex/Data/CORN/C101.2';  spikecd = 102;
name = '/home/horwitz/Cortex/Data/CORN/C174.3';  spikecd = 102;
name = '/home/horwitz/Cortex/Data/KIMBALL/K068.6';  spikecd = 102;

paranum = get_paradigm(name,1000);
indices_whtnsSpatial2;
[foo, index, raster, genparams,e1h,e1v,e1t] = readdata (name,spikecd);
eyesamplerate = get_msperADsample(name);
[ExpectedValues, Variances] = StimStats(genparams, paranum);
data = {};
for j = 1:niter
   for i = 1:length(sigmas)
      e1hrand = e1h;  e1vrand = e1v;
      if (sigmas(i) == 0)
         e1hrand(~isnan(e1h)) = 0;
         e1vrand(~isnan(e1h)) = 0;
      else
         e1hrand(~isnan(e1h)) = normrnd(0,sigmas(i),sum(sum(~isnan(e1h))),1);
         e1vrand(~isnan(e1h)) = normrnd(0,sigmas(i),sum(sum(~isnan(e1v))),1);
      end
%      RFmod('init',{ExpectedValues,Variances,alpha});
%      [out] = getWhtnsStats2(foo,index,raster,genparams,e1hrand,e1vrand,e1t,eyesamplerate,maxT,0,'spike',paranum,'RFmod');
%      n = out{1};
%      framestart = out{2};
%      framestop = out{3};
%      RFmap = out{4};
%      STA = out{5};

      [out] = getWhtnsStats2(foo,index,raster,genparams,e1hrand,e1vrand,e1t,eyesamplerate,maxT,0,'spike',paranum,'Momentmod');

      nframes = foo(:,I_NFRAMES1)+foo(:,I_NFRAMES2).*2^8;
      nframes = max(nframes - maxT + 1, 0);
      nframes = sum(nframes(:));
      in2{1} = out{5};  % n
      in2{2} = 4;
      in2{3} = 12;
      in2{4} = ones(8,8);

      [QIN, v1, v2] = NykampAlpha2(out, in2, Variances, nframes);
      data{i}.i = sigmas(i);
      if (isfield(data{i},'QINs'))
         data{i}.QINs(length(data{i}.QINs)+1) = QIN;
      else
         data{i}.QINs = QIN;
      end
   end
end

tmp = [];
for j = 1:length(data)
   tmp = [tmp; data{j}.i std(data{j}.QINs) mean(data{j}.QINs)];
end

figure; axes; hold on;
for i = 1:size(tmp,1)
   plot(tmp(i,1),tmp(i,3),'k*')
   plot([tmp(i,1) tmp(i,1)],[tmp(i,3)+tmp(i,2) tmp(i,3)-tmp(i,2)],'k-')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Trying to create a confidence ellipsoid
%  around various RGB triplets in the STA
%  Then projecting these ellipsoids onto the
%  unit sphere to determine whether any of them
%  are oppositely tuned to the center.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lin = linspace(0,2*pi,1000);
sine = [cos(lin);sin(lin)];

name = '/home/horwitz/Cortex/Data/KIMBALL/K115.4'
spikecd = 102;
[foo, index, raster, genparams,e1h,e1v,e1t] = readdata (name,spikecd);
paranum = get_paradigm(name,1000);
[ExpectedValues, Variances] = StimStats(genparams, paranum);
if (paranum == 31)
   indices_whtnsSpatial;
else
   indices_whtnsSpatial2;
end

RFmod('init',{ExpectedValues,Variances,0.0000001});
[out] = getWhtnsStats(foo,index,raster,genparams,20,0,'spike',paranum,'RFmod');
RFmap = out{4};
[whichpix, whichgun, whichframe] = pickapix(out{5},[],[],1);

feval('STCOVmod4','init',[genparams(GEN_STIMNELEMX)*genparams(GEN_STIMNELEMY) 3 whichframe whichframe]);
[out] = getWhtnsStats3(foo,index,raster,genparams,whichframe,0,'spike',paranum,'STCOVmod4');

STA = out{1};
STCOV = out{2};
n = out{3};
crit = finv(1-0.01,3,n-3);
crit = crit*((n-1)*3)/(n-3);

figure; axes; hold on;
set(gca,'XLim',[-2 2],'YLim',[-2 2],'Zlim',[-2 2]);
axis square;
axis vis3d;
view(45,45);
for whichpix = find(RFmap)'
   rgbindex = [1:3]+(whichpix-1)*3;
   mu = STA(rgbindex);
   S = STCOV(rgbindex,rgbindex)./n;
   Sproj = projortho(S,mu,'cov');
   [v,d] = eig(Sproj);
   [sorted,order] = sort(diag(d));
   v = v(:,order);
   newd = diag(sorted); newv = v;
   newd(1,:) = []; newd(:,1) = []; 
   newv(:,1) = [];

   greg = newv*sqrt(crit*newd)*sine+repmat(mu,1,size(sine,2));

   [theta,phi,rho] = cart2sph(greg(1,:),greg(2,:),greg(3,:));
%   plot3(greg(1,:),greg(2,:),greg(3,:),'k-');
   g = mkbasis(greg);
   plot3(g(1,:),g(2,:),g(3,:),'k-','LineWidth',ceil(norm(mu*100)));
end

nframes = 35;
M = moviein(nframes);
angles = linspace(0,360,nframes+1);
angles(end) = [];
for i = 1:length(angles)
   view(angles(i),30);
   M(:,i) = getframe;
end
movie(M);

%%% Are there any pixels that are significant and could have
%%% a color tuning that is opposite of the center?
centerpix = ceil(find(STA == max(STA))/3)
rgbindex = [1:3]+(centerpix-1)*3;
mu = mkbasis(STA(rgbindex));

figure; axes; hold on;
set(gca,'XLim',[-2 2],'YLim',[-2 2],'Zlim',[-2 2]);
plot3(mu(1),mu(2),mu(3),'r*');
plot3(-mu(1),-mu(2),-mu(3),'k*');
plot3(0,0,0,'b*');
for whichpix = 1:64
   rgbindex = [1:3]+(whichpix-1)*3;
   mu = STA(rgbindex);
   S = STCOV(rgbindex,rgbindex)./n;
   if (mu'*inv(S)*mu > crit)
      Sproj = projortho(S,mu,'cov');
      [v,d] = eig(Sproj);
      [sorted,order] = sort(diag(d));
      v = v(:,order);
      newd = diag(sorted); newv = v;
      newd(1,:) = []; newd(:,1) = []; 
      newv(:,1) = [];
      greg = newv*sqrt(crit*newd)*sine+repmat(mu,1,size(sine,2));   
      g = mkbasis(greg);
      plot3(g(1,:),g(2,:),g(3,:),'k-','LineWidth',ceil(norm(mu*100)));
   end
end

%%%%%%%%%
%%
%% Trying a few methods to see if the confidence ellipsoid
%% in cone space intersects one of the cardinal planes.
%% 
%%%%%%%%%
alpha = 0.01;
fundamentals = getGlobalParams('fundamentals');
name = '/home/horwitz/Cortex/Data/CORN/C111.4c5'
spikecd = 102;
[foo, index, raster, genparams,e1h,e1v,e1t] = readdata (name,spikecd);
paranum = get_paradigm(name,1000);
if (paranum == 31)
   indices_whtnsSpatial;
else
   indices_whtnsSpatial2;
end

monserial = genparams(GEN_MONSERIAL1)+...
            genparams(GEN_MONSERIAL2)*2^8+...
            genparams(GEN_MONSERIAL3)*2^16;
cal = getGlobalParams('cal',monserial);
M = fundamentals*cal.P_device;
[out] = getWhtnsStats(foo,index,raster,genparams,20,0,'spike',paranum,'STAmod');
STA = out{1};
[centerpix, whichgun, whichframe] = pickapix(STA,[],[],1);

feval('STCOVmod4','init',[genparams(GEN_STIMNELEMX)*genparams(GEN_STIMNELEMY) 3 whichframe whichframe]);
[out] = getWhtnsStats3(foo,index,raster,genparams,whichframe,0,'spike',paranum,'STCOVmod4');

STCOV = out{2};
n = out{3};

data = [];
for whichpix = 1:size(STA,1)
   STArgb = squeeze(STA(whichpix,:,whichframe));
   STAlms = inv(M')*STArgb';
   STCrgb = STCOV(3*(whichpix-1)+[1 2 3],3*(whichpix-1)+[1 2 3]);
   STCrgb = STCrgb./n;
   STClms = inv(M')*STCrgb*inv(M);

   crit = finv(1-alpha,3,n-3);
   crit = crit*((n-1)*3)/(n-3);
   bounds = FindEllipsoidExtrema(STClms,crit);

   statistic = STAlms'*inv(STClms)*STAlms;
   p = 1-fcdf(statistic * (n-3)/((n-1)*3),3, n-3);
   data = [data; p (abs(STAlms') > bounds) .* sign(STAlms')];
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Testing
%%%
%%% Note: this eigenvector decomposition
%%% thing only works for well-conditioned
%%% covariance matrices.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
crit = finv(1-alpha,3,n-3);
crit = crit*((n-1)*3)/(n-3);

% [v,d] = eig(STClms);
% scalefactors = sqrt(diag(d).*crit);
% newd = diag(scalefactors);
[u,s,v] = svd(STClms);
news = diag(sqrt(diag(s)*crit))
points = [];
for i = 1:2000
    vect = mkbasis(normrnd(0,1,3,1));
%    vect = normrnd(0,1,3,1);
%    points = [points; (v*newd*vect)'];
    points = [points; ((u*news*v)*vect)'];
end
bounds = FindEllipsoidExtrema(STClms,crit);
figure;
plot3(points(:,1),points(:,2),points(:,3),'k.');
hold on;
plot3(bounds(1),0,0,'m*'); plot3(-bounds(1),0,0,'m*');
plot3(0,bounds(2),0,'g*'); plot3(0,-bounds(2),0,'g*');
plot3(0,0,bounds(3),'b*'); plot3(0,0,-bounds(3),'b*');
xlabel('x'); ylabel('y'); zlabel('z');
max(points)
bounds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L = data(:,1) < alpha;
data(L,:)

RF = zeros(8,8)
RF(find(L)) = 1;
figure;
imagesc(RF);

%% Definition of single-opponency
%% There has to be a sign flip in at least one cone type
%% at one pixel.

centersigns = data(centerpix,[2:4]);

testmat = repmat(centersigns,sum(L),1).*data(L,[2:4])
noppositepix = sum(any(testmat == -1,2))
npix = sum(L)-1   %% '-1' because we don't want to count the center pixel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Simulating a grating experiment using the 
%  "ratio of two very high dimensional Gaussians"
%  on a single cell.  The question is: can we use
%  something like this to sort simple and complex
%  cells.
%
%  The problem with this approach is inverting the 
%  spike-triggered covariance matrix (which you 
%  need to do to solve for Pr(stim|spike), assuming
%  a Gaussian density) doesn't seem to work!
%  Numerical error?
%
%  Another approach -- chucking elements that are
%  "not in the RF" and combining guns.  Roundoff 
%  errors are still killing this.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha = 0.00001;
maxT = 20;
path = '/home/horwitz/Cortex/Data';
filename = 'CORN/C191.1';
spikecd = 1;
disp([filename,': ',num2str(spikecd)]);
paranum = get_paradigm([path,'/',filename],1000);
[foo, index, raster, genparams] = readdata([path,'/',filename],spikecd);
[ExpectedValues, Variances] = StimStats(genparams, paranum);
if (paranum == 60)
   indices_whtnsSpatial2;
else
   indices_whtnsSpatial;
end

RFmod('init',{ExpectedValues,Variances,alpha});
[out] = getWhtnsStats(foo,index,raster,genparams,maxT,0,'spike',paranum,'RFmod');
nspikes = out{1};
startframe = out{2};
endframe = out{3};
nframes = endframe-startframe+1;
RF = out{4};

% Getting rid of stuff that's "not in RF"

mask = logical([repmat(zeros(64,1), maxT-endframe,1);...
                repmat(RF(:), endframe-startframe+1,1);...
                repmat(zeros(64,1), startframe-1,1)]);

gunmat = toeplitz([1; zeros(64*3*maxT-1,1)],[ones(1,3) zeros(1,64*3*maxT-3)]);
gunmat = gunmat([1:3:end],:);

feval('STCOVmod3','init',[genparams(GEN_STIMNELEMX)*genparams(GEN_STIMNELEMY) 3 maxT]);
[out] = getWhtnsStats3(foo,index,raster,genparams,endframe,0,'spike',paranum,'STCOVmod3');
STA_ns = gunmat*out{1};
STA_ns = STA(mask);
STCOV_ns = gunmat*out{2}*gunmat';
STCOV_ns = STCOV(mask,mask);

feval('STCOVmod3','init',[genparams(GEN_STIMNELEMX)*genparams(GEN_STIMNELEMY) 3 maxT]);
[out] = getWhtnsStats3(foo,index,raster,genparams,endframe,0,'all',paranum,'STCOVmod3');
STA_na = gunmat*out{1};
STA_na = STA(mask);
STCOV_na = gunmat*out{2}*gunmat';
STCOV_na = STCOV(mask,mask);


grating = mkgrating(8,8,pi/4,3,pi/2);
grating = repmat(grating(:),maxT,1);

%% Scaling all the guns to have variance 1
 
varmat = repmat(sqrt(sum(Variances)),maxT*64,1)*repmat(sqrt(sum(Variances)),maxT*64,1)';
STCOV = STCOV./varmat(mask,mask);
STA = STA./sqrt(varmat(mask,1));
invSTC = inv(STCOV);

a = (exp(-.5*(grating(mask)-STA)'*invSTC*(grating(mask)-STA)))/sqrt(det(STCOV))
b = exp(-.5*(grating(mask))'*eye(size(STCOV))*(grating(mask)))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Experimenting with non-causal estimation of 
%  1st and 2nd order kernels.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maxT = 20;
path = '/home/horwitz/Cortex/Data';
filename = 'CORN/C101.2';
spikecd = 102;
paranum = get_paradigm([path,'/',filename],1000);
[foo, index, raster, genparams] = readdata([path,'/',filename],spikecd);
[ExpectedValues, Variances] = StimStats(genparams, paranum);
if (paranum == 60)
   indices_whtnsSpatial2;
else
   indices_whtnsSpatial;
end

[nsp, dt] = countspikes(raster, index, SIN_STIMON, 0, SIN_STIMOFF,0);
fr = sum(nsp)./sum(dt)*10;
feval('STCOVmod7','init',[genparams(GEN_STIMNELEMX)*genparams(GEN_STIMNELEMY) 3 4 12 fr]);
[out] = getWhtnsStats3(foo,index,raster,genparams,0,0,'spike',paranum,'STCOVmod7');
v1 = out{1}*out{1}';
v1 = v1(:);
v2 = out{2}(:);

feval('STCOVmod7','init',[genparams(GEN_STIMNELEMX)*genparams(GEN_STIMNELEMY) 3 -12 -6 fr]);
[out] = getWhtnsStats3(foo,index,raster,genparams,0,0,'spike',paranum,'STCOVmod7');
noise = out{2}(:);

scalefactor = v1\v2;
k = length(v1);
chi2resid = sum(((v2-scalefactor*v1)/std(noise)).^2);
Zresid = (chi2resid - k)/sqrt(2*k);
chi2sig = sum((v2/std(noise)).^2);
Zsig = (chi2sig - k)/sqrt(2*k);
Rln  = (Zsig-Zresid)/Zsig;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The direction that is orthogonal to the "luminance"
% and "L-M" directions of complex cells.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = [0.0518 0.1183 0.0277; 0.0184 0.1215 0.0401; 0.0018 0.0089 0.1451];
% In RGB space 
v1 = [-0.3865  -0.3865 0.8374]';
v2 = [0.2476  0.8666 0.4333]';
v3 = cross(v1,v2);
v3lms = M*v3;
v3lms = v3lms./norm(v3lms)

% Doing it another way
v11 = inv(M')*v1;
v22 = inv(M')*v2;
v33 = cross(v11,v22);
v33 = v33./norm(v33)

%%%%%%%%%%%%%%%%%%%%%%%%
%
% Testing the new getWhtnsStim
%
%%%%%%%%%%%%%%%%%%%%%%%%
datapath = '/home/horwitz/Cortex/Data'
filename = 'DONDERS/D055.1';
spikecd = 1;
[foo, index, raster, genparams] = readdata([datapath,'/',filename],spikecd);
paranum = get_paradigm([datapath,'/',filename],1000);

maxT = 20;
[SpikeStim,TrialVect,TimeVect,Lspike] = getWhtnsStim(foo,index,raster,genparams,maxT,0,'spike',paranum,[6 4]);

maxT = 20;
[SpikeStim,TrialVect,TimeVect,Lspike] = getWhtnsStim(foo,index,raster,genparams,maxT,0,'spike',paranum,[]);

SpikeStim = vectorizestim(permute(SpikeStim,[1 2 4 3]));
% For getting SpikeStim with repeats
for k = 2:max(Lspike)
   SpikeStim = [SpikeStim; SpikeStim(Lspike >= k,:)];
end
imboost(reshape(mean(SpikeStim),[8 8 3]))
S = cov(SpikeStim);
figure;
imagesc(S);

feval('STCOVmod4','init',[genparams(GEN_STIMNELEMX)*genparams(GEN_STIMNELEMY) 3 whichframe whichframe]);
[out] = getWhtnsStats3(foo,index,raster,genparams,whichframe,0,'spike',paranum,'STCOVmod4');
STCOV = out{2};
figure;
imagesc(STCOV);

figure;
imagesc(S-STCOV);

figure; axes; hold on;
plot(mean(SpikeStim),'k-')
tmp = permute(reshape(out{1},3,8,8),[2 3 1]);
plot(tmp(:),'mo')

[v1,d1,flag] = eigs(S,1,'LM',OPTIONS);
[v2,d1,flag] = eigs(STCOV,1,'LM',OPTIONS);
figure; axes; hold on;
plot(v1,'k-')
tmp = permute(reshape(v2,3,8,8),[2 3 1]);
plot(tmp(:),'mo')

imboost(reshape(v1,[8 8 3]))

%%%%%%%%%%%%%%%%%%%%%%%%
%
% Testing out the paradigm on Jude's rig
%
%%%%%%%%%%%%%%%%%%%%%%%%

name = '/home/horwitz/Cortex/Data/MAX/M004.1'
paranum = get_paradigm(name,1000);
[foo, index, raster, genparams,e1h,e1v,e1t] = readdata (name,3);
%% First just looking at the spike rate
indices_whtnsSpatial2;

[nsp1, dt1] = countspikes(raster, index, SIN_FIX, 0, SIN_STIMON,0);
[nsp2, dt2] = countspikes(raster, index, SIN_STIMON, 0, SIN_STIMOFF,0);
fr1 = nsp1./dt1*1000;
fr2 = nsp2./dt2*1000;
figure; axes; hold on;
plot(fr1, fr2,'k.');
plot([0 30],[0 30],'k-');
axis equal
xlabel('Baseline firing rate (sp/sec)');
ylabel('Driven firing rate (sp/sec)');

%framerate = 14.31;  % Xin's monitor
framerate = 8.327;  % Jude's monitor
stimdurations = index(:,SIN_STIMOFF)-index(:,SIN_STIMON);
nframes = foo(:,I_NFRAMES1)+foo(:,I_NFRAMES2).*2^8;
expectednframes = stimdurations/framerate;
residnframes = nframes-expectednframes;
figure
plot(nframes,residnframes,'k.')

[out] = getWhtnsStats(foo,index,raster,genparams,20,0,'spike',...
                      paranum,'STAmod');

STA = out{1};
STV = out{2};
figure;
subplot(2,1,1);
plot(STA(:));
subplot(2,1,2);
plot(STV(:));

%% Looking at frame skips
[n, nm, nl, EventTimes, names, EventCodes] = nex_marker2(name, 'Strobed');
L = EventCodes == 119;
starts = find(diff(L) == 1);
ends = find(diff(L) == -1);
data = [];
figure; axes; hold on;
for i = 1:length(starts)
   [EventCodes(starts(i)+1:ends(i))];
   delta_t = diff(EventTimes(starts(i)+1:ends(i)));
   if (~isempty(delta_t))
      plot(delta_t*1000);
      data = [data; delta_t'];
   end
end

sum(abs(data-median(data)) < .001)
sum(abs(data-median(data)) > .001)

1/median(data)  % Frequency in Hz

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Running an L-N-P model
%  that tests the hypothesis that
%  color/time separability (in the absence
%  of color/space separability) leads to
%  color time inseparability when tested
%  with a spatially-uniform stimulus.
%
%  Can't use real data because the naive separable
%  model fit has a very noisy STA when averaged 
%  across pixels.
%
%  This is kind of trivial, but at least it makes
%  the point.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CONSTANTS
maxT = 20;
centerRGB = [.2 .2 -2];
surroundRGB = [.2 .3 -.1];
centerRGB = centerRGB./norm(centerRGB);
surroundRGB = surroundRGB./norm(surroundRGB);
centersigma = 0.5; % totally irrelevant
surroundsigma = 1; % totally irrelevant
surroundweight = 1;
centertimeparams = [3,1.5];
surroundtimeparams = [3,3];
nblocks = 10;
nstim = 10;

center = normpdf(linspace(-3,3,8),0,centersigma);
centermat = center'*center;
centertimecourse = gampdf([0:maxT-1],centertimeparams(1),centertimeparams(2));
centervolume = reshape(centermat(:)*centertimecourse,[8,8,maxT]);

surround = normpdf(linspace(-3,3,8),0,surroundsigma);
surroundmat = surround'*surround;
surroundtimecourse = gampdf([0:maxT-1],surroundtimeparams(1),surroundtimeparams(2));
surroundvolume = reshape(surroundmat(:)*surroundtimecourse,[8,8,maxT]);

center_t = squeeze(sum(sum(centervolume)))*centerRGB;
surround_t = squeeze(sum(sum(surroundvolume)))*surroundRGB;

%% Finding the theoretical STA
STA_theory = center_t-surroundweight*surround_t;
[u,s,v] = svd(STA_theory);
figure;
subplot(2,1,1);
plot(STA_theory);
subplot(2,1,2);
plot(u(:,1)*s(1)*v(:,1)');

center_v = center_t(:);
surround_v = surround_t(:);

STS = zeros(size(center_v));
for i = 1:nblocks
   stimulus = normrnd(0,1,maxT*3,nstim);
   center_dotprod = center_v'*stimulus;
   surround_dotprod = surround_v'*stimulus;
   dotprods = center_dotprod-(surroundweight*surround_dotprod);
   dotprods(dotprods < 0) = 0;  % halfwave rectification
   STS = STS + stimulus*dotprods';
end

figure;
plot(reshape(STS,[maxT, 3]));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Doing a L-N-P model fit to determine
%  whether Johnson et al. would have classified
%  a given cell as "luminance", "color-luminance",
%  or "color".
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

maxT = 20;
nbins = 6;
niter = 20;
CVprop = .5;
framedur = getGlobalParams('framedur');
datapath = '/home/horwitz/Cortex/Data';
fundamentals = getGlobalParams('fundamentals');

filename = 'CORN/C101.2'; spikecd = 102;
filename = 'KIMBALL/K092.3'; spikecd = 102;

[foo, index, raster, genparams] = readdata([datapath,'/',filename],spikecd);
paranum = get_paradigm([datapath,'/',filename],1000);
[ExpectedValues, Variances] = StimStats(genparams, paranum);
if (paranum == 60)
   indices_whtnsSpatial2;
else
   indices_whtnsSpatial;
end
[out] = getWhtnsStats(foo,index,raster,genparams,maxT,0,'spike',paranum,'STAmod');
STA = out{1};
STV = out{2};
n = out{3};

[whichpix,whichgun,whichframe] = pickapix(STA,STV,n,6,Variances);

bins = linspace(-.02,.02,20);
Ns = []; As = [];
for iter = 1:niter
   iter
   shuffledtrials = randperm(size(foo,1));
   basistrials = shuffledtrials(1:floor(size(foo,1)*CVprop));
   projtrials = shuffledtrials([ceil(size(foo,1)*CVprop):end]);
   [out] = getWhtnsStats(foo(basistrials,:),index(basistrials,:),raster(basistrials,:),genparams,maxT,0,'spike',paranum,'STAmod');
   STA = squeeze(out{1}(:,:,whichframe));
   STA = STA';

   stimdurations = index(projtrials,SIN_STIMOFF) - index(projtrials,SIN_STIMON);
   STPROJmod('init',STA(:),whichframe-1,sum(stimdurations)/10,[64 3 1])
   [out] = getWhtnsStats3(foo(projtrials,:),index(projtrials,:),raster(projtrials,:),genparams,whichframe,0,'spike',paranum,'STPROJmod');
   proj = out{1}; 
   Lspike = out{2};
   [A,X] = hist(proj,bins);
   [N,X] = hist(proj(Lspike==1,:),bins);
   for j = 2:max(Lspike)
      [tmp1,X] = hist(proj(Lspike==j,:),bins);
      N = N+j*tmp1;
   end
   Ns(:,iter) = N;
   As(:,iter) = A;
end

meanNs = sum(Ns,2);
meanNa = sum(As,2);
fr = (meanNs./meanNa)*1000/framedur;  % Spikes per second

%%% Generating luminance and red/green isoluminant stimuli
vlambda = getGlobalParams('vlambda');
monserial = genparams(GEN_MONSERIAL1)+...
            genparams(GEN_MONSERIAL2)*2^8+...
            genparams(GEN_MONSERIAL3)*2^16;
cal = getGlobalParams('cal',monserial);
gunlum = vlambda*cal.P_device;
gamma = cal.gammaTable;
bkgndRGB = genparams([GEN_BKGNDR; GEN_BKGNDG; GEN_BKGNDB]);
bkgndrgb = [gamma(bkgndRGB(1)-1,1), gamma(bkgndRGB(2)-1,2), gamma(bkgndRGB(3)-1,3)];
M = fundamentals*cal.P_device;
bkgndlms = M*bkgndrgb';
bkgndlum = bkgndlms(1)+bkgndlms(2);  % Assuming luminance is given by L+M
lumrgb = [bkgndrgb.*.1+bkgndrgb; bkgndrgb.*-.1+bkgndrgb];
% Two rows indicate peak and trough of the grating, respectively

%%% Testing
%tmp = gunlum*lumrgb';
%(tmp(1)-tmp(2))/mean([tmp(1) tmp(2)])

% Maxing the red gun out and adjusting green to make an nice-oluminant stimulus
lumratio = gunlum(1)/gunlum(2);
newr = [.18; -.18]+bkgndrgb(1);
newg = bkgndrgb(2)-((newr-bkgndrgb(1)).*lumratio);
equirgb = [newr newg [bkgndrgb(3); bkgndrgb(3)]];

%%% Testing 
%gunlum*equirgb'
%equilms = M*equirgb';
%(equilms(:,1)-equilms(:,2))./bkgndlms

%%%% Gotta get the orientation and spatial frequency
[u,s,v] = svd(STA);
RF = reshape(v(:,1),[8 8]);
fftim = fft2(RF)
tmp = find(abs(fftim) == max(abs(fftim(:))))
[ysf,xsf] = ind2sub([8 8],min(tmp));
phase = angle(fftim(min(tmp)));

xinc = (2*pi/8)*(xsf-1);
yinc = (2*pi/8)*(ysf-1);
[xramp, yramp] = meshgrid(xinc*([0:7]), yinc*([0:7]));
grating = cos(xramp+yramp+phase);

%%%% Now doing the actual simulation
lumprojs = lumrgb-repmat(bkgndrgb,2,1);
equiprojs = equirgb-repmat(bkgndrgb,2,1);

lumgrating = grating(:)*lumprojs(1,:);
equigrating = grating(:)*equiprojs(1,:);

lumgensig = sum(sum(lumgrating.*STA'))
equigensig = sum(sum(equigrating.*STA'))

nsubbins = 100;
nbins = length(bins);
zi = interp1(bins,fr,linspace(min(bins),max(bins),nsubbins));
subbins = linspace(min(bins),max(bins),nsubbins)

lumresp1 = zi(min(find(subbins > lumgensig)));
lumresp2 = zi(max(find(subbins < -lumgensig)));

equiresp1 = zi(min(find(subbins > equigensig)));
equiresp2 = zi(max(find(subbins < -equigensig)));

abs(equiresp1-equiresp2)/abs(lumresp1-lumresp2)

%%% Johnson's cell sorting criteria
%%% < 0.5 is luminance
%%% > 0.5  and < 2 is color-luminance
%%% > 2 is color

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Looking at conditional STAs for a blue-yellow
%  cell (positive versus negative projections 
%  onto PC1).  Testing the hypothesis that 
%  *really* what drives the cell is a contrast
%  edge with blue consistently on one side or the other.
%  This hypothesis seems to be wrong: the conditional
%  STAs look a lot like the regular STAs.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

name = '/home/horwitz/Cortex/Data/DONDERS/D034.2c6'
spikecd = 1;
paranum = get_paradigm(name,1000);
[foo, index, raster, genparams,e1h,e1v,e1t] = readdata (name,spikecd);
if (paranum == 60)
   indices_whtnsSpatial2;
else
   indices_whtnsSpatial;
end

whichframe = 7;
nframes = ceil(sum(index(:,SIN_STIMOFF)-index(:,SIN_STIMON))/10);
AllStimDim =  [nframes 64 3 1];
AllStimmod2('init',{whichframe,AllStimDim});

[out] = getWhtnsStats(foo,index,raster,genparams,20,0,'spike',...
                      paranum,'AllStimmod2');

SpikeStim = out{1};
Lspike = out{2};
SpikeStim = reshape(SpikeStim,size(SpikeStim,1), 64*3);
STA = mean(SpikeStim);
STC = cov(SpikeStim);
[pc1,d] = eigs(projortho(STC,STA,'cov'),1);

% Debugging
imboost(reshape(STA,8,8,3))
imboost(reshape(pc1,8,8,3))
projs = SpikeStim*pc1;
L = logical(projs > 0);
STA1 = mean(SpikeStim(L,:));
STA2 = mean(SpikeStim(~L,:));

figure;
subplot(2,1,1);
imboost(reshape(projortho(STA1,pc1),8,8,3));
subplot(2,1,2);
imboost(reshape(projortho(STA2,pc1),8,8,3));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculating the distribution of the 
% inseparability index under Ho.
% Take STA, fit separable model,
% calculate distance.  Then, use
% the separable model fit as the new 
% STA, translate the error ellipsoid
% over to it, and pick new "separable
% fits" and calculate their distances.
%
% To do this "for real" would require 
% using the correct "whichframe".
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

name = '/home/horwitz/Cortex/Data/DONDERS/D034.2c6'
spikecd = 1;
maxT = 20;
[foo, index, raster, genparams,e1h,e1v,e1t] =...
             readdata(name, spikecd);
paranum = get_paradigm(name,1000);
if (paranum == 60)
   indices_whtnsSpatial2;
elseif (paranum == 31)
   indices_whtnsSpatial;
else %  (paranum == 61)   % S part of U and S white noise
   indices_UandSwhtns;
   L = logical(foo(:,I_STIMTYPE) == 0);
   foo = foo(L,:); index = index(L,:); raster = raster(L,:);
   paranum = 62;
end
         
[out] = getWhtnsStats(foo,index,raster,genparams,maxT,0,'spikes',paranum,'STAmod');
STA = out{1};
STV = out{2};
nspikes = out{3};
[whichpix,whichgun,whichframe] = pickapix(STA,STV,nspikes,1);
STA_st = squeeze(STA(whichpix,:,:));
STA_ss = squeeze(STA(:,:,whichframe));

%%% Getting the Mahalanobis distance for the
%%% checkerboard stimulus in the spatial domain
ndims = 3*genparams(GEN_STIMNELEMX)*genparams(GEN_STIMNELEMY);
feval('STCOVmod4','init',[genparams(GEN_STIMNELEMX)*genparams(GEN_STIMNELEMY) 3 whichframe whichframe]);
[out] = getWhtnsStats3(foo,index,raster,genparams,whichframe,0,'spike',paranum,'STCOVmod4');
meanvect = out{1}';
nstim = out{3};
STCOV = out{2}/nstim;

meanmat = reshape(meanvect,3,genparams(GEN_STIMNELEMX)*genparams(GEN_STIMNELEMY));
[u,s,v] = svd(meanmat);
predmat = u(:,1)*s(1,1)*v(:,1)';
predvect = reshape(predmat,1,ndims);
dist_ss = sqrt((meanvect-predvect)*inv(STCOV)*(meanvect-predvect)')

niter = 2000;
%[v,d] = eig(STCOV);
%transform = v*sqrt(d);
transform = sqrtm(STCOV);
data = [];
STAs = [];
for i = 1:niter
   tmp = normrnd(0,1,ndims,1);
   noise = transform*tmp;
   newSTA = noise'+predvect;
   meanmat = reshape(newSTA,3,genparams(GEN_STIMNELEMX)*genparams(GEN_STIMNELEMY));
   [u,s,v] = svd(meanmat);
   simpredmat = u(:,1)*s(1,1)*v(:,1)';
   simpredvect = reshape(predmat,1,ndims);
   data = [data; sqrt((newSTA-simpredvect)*inv(STCOV)*(newSTA-simpredvect)')];
   STAs = [STAs; newSTA];
end

hist(data./dist_ss)

TestCOV = cov(STAs);
imagesc(STCOV-TestCOV)
plot(STCOV(:)-TestCOV(:))
plot(diag(STCOV)); hold on; plot(diag(TestCOV),'m-');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculating number of cones per pixel
% in the stimulus.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
conespermm2 = 10000;  % at ~5-10 deg ecc.
degperpix = .22;
degpermm = 4;
mmperdeg = 1/degpermm;
mmperpix = mmperdeg*degperpix;
mm2perpix = mmperpix^2;
ncones = conespermm2*mm2perpix
propS = .07;
propL = .6*(1-propS);
propM = .4*(1-propS);

ExpectedL = propL*ncones
ExpectedM = propM*ncones
ExpectedS = propS*ncones

% Assuming a fixed number of S cones and total cones

LandMcones = round(ncones)-round(ExpectedS)

% Probability as a function of L cone proportion
p = [];
x =  [0:LandMcones]
for i = x
   p = [p; binopdf(i,LandMcones,propL/(1-propS))];
end
plot(x/LandMcones,p,'r-');
xlabel('prop. L cones')

nL = binoinv([.025 .975],LandMcones,propL/(1-propS))
nM = LandMcones-nL
nL(1)/(nL(1)+nM(1))
nL(2)/(nL(2)+nM(2))
% The bottom line seems to be that the ratio of L to M cones
% in the pixel is not terribly variable (42% to 79% L - for 95% of pixels).
% The expected number of S cones is small ~2, so there could be an issue here.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Investigating separability
% within and between subunits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load /home/horwitz/matlab/PopAnal/Data/SeparabilityLatencies
datapath = '/home/horwitz/Cortex/Data/';

names = {'KIMBALL/K025.4c5','KIMBALL/K128.4','TIGER/T025.4'};
for j = 1:length(names)

name = char(names{j});
whichframe = 0;
for i = 1:length(data)
   if (strcmp(name, data{i}.filename))
      whichframe = data{i}.whichframe;
      spikecd = data{i}.spikecode;
   end
end

maxT = 20;
[foo, index, raster, genparams,e1h,e1v,e1t] =...
             readdata([datapath,name], spikecd);
paranum = get_paradigm([datapath,name],1000);
if (paranum == 60)
   indices_whtnsSpatial2;
elseif (paranum == 31)
   indices_whtnsSpatial;
else %  (paranum == 61)   % S part of U and S white noise
   indices_UandSwhtns;
   L = logical(foo(:,I_STIMTYPE) == 0);
   foo = foo(L,:); index = index(L,:); raster = raster(L,:);
   paranum = 62;
end
         
[out] = getWhtnsStats(foo,index,raster,genparams,maxT,0,'spikes',paranum,'STAmod');
STA = out{1};
STV = out{2};
nspikes = out{3};
STA_ss = squeeze(STA(:,:,whichframe));
figure
plot3(STA_ss(:,1),STA_ss(:,2),STA_ss(:,3),'k.')

labels = {'R weight','G weight','B weight'};
order = [1 2; 2 3; 1 3];
figure;
for i = 1:size(order,1);
   subplot(2,2,i); hold on;
   tmp = STA_ss(:,order(i,:));
   plot(tmp(:,1),tmp(:,2),'k.');
   lims = [min(min(tmp(:))) max(max(tmp(:)))];
   axis square; set(gca,'XLim',lims','YLim',lims);
   plot([0 0],[lims(1) lims(2)],'k:');
   plot([lims(1) lims(2)],[0 0],'k:');
   xlabel(labels(order(i,1)));
   ylabel(labels(order(i,2)));
   if (i == 1)
      title(name)
   end
end
subplot(2,2,4)
imboost(reshape(STA_ss,8,8,3))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Doing an analysis of *unsigned* color tuning for each pixel.
%%%% Projecting the RGB of all the pix onto the RGB for the center
%%%% pixel, removing the signs of the projections, and seeing 
%%%% what's where in the spatial map.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


load /home/horwitz/matlab/PopAnal/Data/CenSurr.mat
datapath = '/home/horwitz/Cortex/Data/';

name = 'KIMBALL/K025.4c5'
spikecd = 102;

whichframe = 0;
for i = 1:length(data)
   if (strcmp(name, data{i}.filename))
      whichframe = data{i}.whichframe;
      celldata = data{i};
   end
end

maxT = 20;
[foo, index, raster, genparams,e1h,e1v,e1t] =...
             readdata([datapath,name], spikecd);
paranum = get_paradigm([datapath,name],1000);
if (paranum == 60)
   indices_whtnsSpatial2;
elseif (paranum == 31)
   indices_whtnsSpatial;
else %  (paranum == 61)   % S part of U and S white noise
   indices_UandSwhtns;
   L = logical(foo(:,I_STIMTYPE) == 0);
   foo = foo(L,:); index = index(L,:); raster = raster(L,:);
   paranum = 62;
end
         
[out] = getWhtnsStats(foo,index,raster,genparams,maxT,0,'spikes',paranum,'STAmod');
STA = out{1};
STV = out{2};
nspikes = out{3};
STA_ss = squeeze(STA(:,:,whichframe));

centerpixRGB = STA_ss(find(sum(STA_ss.^2,2) == max(sum(STA_ss.^2,2))),:)
centerpixRGB = centerpixRGB./norm(centerpixRGB);
projs = STA_ss*centerpixRGB';
correlations = projs./sqrt(sum(STA_ss.^2,2))
centercors = correlations(find(celldata.Center),:);
surroundcors = correlations(find(celldata.Surround),:);

bins = linspace(.8,1,20);
[n1,x] = hist(centercors(centercors<1),bins);
[n2,x] = hist(-surroundcors,bins);

figure;
subplot(2,1,1);
bar(x,[n1+n2],'r');   % Red is for center
hold on;
bar(x,n2,'k');        % Black is for surround
axis square;
set(gca,'XLim',[min(bins) 1])
xlabel('Correlation with center pixel');
ylabel('Count');

subplot(2,1,2);
imboost(reshape(STA_ss,8,8,3));
hold on;
[i,j] = ind2sub([8 8],find(celldata.Center))
plot(j,i,'rx')
[i,j] = ind2sub([8 8],find(celldata.Surround))
plot(j,i,'bo')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Comparing unsigned cone weights between center and surround
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load /home/horwitz/matlab/PopAnal/Data/CenSurr.mat
datapath = '/home/horwitz/Cortex/Data/';

name = 'KIMBALL/K104.5'

for i = 1:length(data)
   if (strcmp(name, data{i}.filename))
      whichframe = data{i}.whichframe;
      celldata = data{i};
   end
end

cenweights = celldata.coneweights(find(celldata.Center),:)
surweights = celldata.coneweights(find(celldata.Surround),:)
cenweights = cenweights(:,10)./sum(abs(cenweights(:,[8 9 10])),2);
surweights = surweights(:,10)./sum(abs(surweights(:,[8 9 10])),2);
%cenweights = cenweights(:,8)./sum(abs(cenweights(:,9)),2);
%surweights = surweights(:,8)./sum(abs(surweights(:,9)),2);

bins = linspace(min([cenweights; -surweights]),...
                max([cenweights; -surweights]),20);
[n1,x] = hist(cenweights,bins)
[n2,x] = hist(-surweights,bins)
figure; axes; hold on;
bar(x,n1+n2,'r')
bar(x,n1,'k')
[h,p] = ttest2(cenweights,surweights);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Looking at the distribution of directions of 
% residuals from the separable model fit.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load /home/horwitz/matlab/PopAnal/Data/CenSurr.mat
datapath = '/home/horwitz/Cortex/Data/';

name = 'KIMBALL/K025.4c5'
spikecd = 102;

whichframe = 0;
for i = 1:length(data)
   if (strcmp(name, data{i}.filename))
      whichframe = data{i}.whichframe;
      celldata = data{i};
   end
end

maxT = 20;
[foo, index, raster, genparams,e1h,e1v,e1t] =...
             readdata([datapath,name], spikecd);
paranum = get_paradigm([datapath,name],1000);
if (paranum == 60)
   indices_whtnsSpatial2;
elseif (paranum == 31)
   indices_whtnsSpatial;
else %  (paranum == 61)   % S part of U and S white noise
   indices_UandSwhtns;
   L = logical(foo(:,I_STIMTYPE) == 0);
   foo = foo(L,:); index = index(L,:); raster = raster(L,:);
   paranum = 62;
end

ndims = 3*genparams(GEN_STIMNELEMX)*genparams(GEN_STIMNELEMY);
feval('STCOVmod4','init',[genparams(GEN_STIMNELEMX)*genparams(GEN_STIMNELEMY) 3 whichframe whichframe]);
[out] = getWhtnsStats3(foo,index,raster,genparams,whichframe,0,'spike',paranum,'STCOVmod4');
meanvect = out{1}';
nstim = out{3};
STCOV = out{2}/nstim;

meanmat = reshape(meanvect,3,genparams(GEN_STIMNELEMX)*genparams(GEN_STIMNELEMY));
[u,s,v] = svd(meanmat);
predmat = u(:,1)*s(1,1)*v(:,1)';
predvect = reshape(predmat,1,ndims);
dist_ss = sqrt((meanvect-predvect)*(meanvect-predvect)')

niter = 10000;
simdata = zeros(niter, 1);
resids = zeros(niter, 192);
STCOVsqrt = sqrtm(STCOV);
for i = 1:niter
   residvect = (STCOVsqrt*normrnd(0,1,length(predvect),1))';
   meanvect = predvect+residvect;  % Predvect comes from the real data
   meanmat = reshape(meanvect,3,genparams(GEN_STIMNELEMX)*genparams(GEN_STIMNELEMY));
   [u,s,v] = svd(meanmat);
   simpredmat = u(:,1)*s(1,1)*v(:,1)';
   simpredvect = reshape(simpredmat,1,ndims);
   residvect = meanvect-simpredvect;
   dist = sqrt((residvect)*(residvect)');
   simdata(i) = dist;
   resids(i,:) = residvect;
end

hist(simdata)
S = cov(resids);
imagesc(S)
plot(diag(S))
[v,d] = eig(S);
sum(diag(d) > mean(diag(d)))
% Interesting.  All the residuals are confined to a 126-D subspace.
% presumably because the separable model fit has 63*2+1 = 127 df.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Trying linear transformations of the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AllStimmod('init',{celldata.centerpix,[sum(sum(~isnan(raster))) 1 3 20]});
[out] = getWhtnsStats(foo,index,raster,genparams,maxT,0,'spike',...
              paranum,'AllStimmod');
SpikeStim = out{1};
SpikeStim = squeeze(vectorizestim(permute(SpikeStim,[1 2 4 3])));
nstim = size(SpikeStim,1);
S = cov(SpikeStim)./nstim;   % Covariance of the mean
meanvect = mean(SpikeStim);
meanmat = reshape(meanvect,maxT,3);
[u,s,v] = svd(meanmat);
predmat = u(:,1)*s(1,1)*v(:,1)';
predvect = reshape(predmat,1,3*maxT);

out = [];
for i = 1:100
   i
   if (i < 50)
      Afact = nan;
      A = eye(60);
   else
      Afact = unifrnd(0,2);
      A = normrnd(0,Afact,length(meanvect),length(meanvect));
   end
   Amean = (A*meanvect')';
   Apred = (A*predvect')';
   realdist = sqrt((Amean-Apred)*inv(A*S*A')*(Amean-Apred)');
%   realdist = sqrt((Amean-Apred)*(Amean-Apred)');

   sS = sqrtm(S);
   invSA = inv(A*S*A');
   tmpdist = [];
   for i = 1:niter
      noise = normrnd(0,1,length(Apred),1);
      sta_mat = reshape(predvect'+sS*noise,length(predvect)/3,3);
      [u,s,v] = svd(sta_mat);
      pred_mat = u(:,1)*s(1,1)*v(:,1)';
      %tmpdist(i) = sqrt((sta_mat(:)-pred_mat(:))'*(sta_mat(:)-pred_mat(:)));
      %tmpdist(i) = sqrt((sta_mat(:)-pred_mat(:))'*invSA*(sta_mat(:)-pred_mat(:)));

      % This does not work (below)
      % Amean = (A*sta_mat(:))';
      % Apred = (A*pred_mat(:))';
      % tmpdist(i) = sqrt((Amean-Apred)*(Amean-Apred)');

      % This works (below)
      Amean = (A*sta_mat(:))';
      Apred = (A*pred_mat(:))';
      tmpdist(i) = sqrt((Amean-Apred)*invSA*(Amean-Apred)');
   end

   stat = (realdist-mean(tmpdist))/std(tmpdist);
   out = [out; Afact realdist mean(tmpdist) std(tmpdist) stat];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Comparing STAs recorded under multiple conditions.
%  Binconen.stt.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filename = '/home/horwitz/Cortex/Data/ROLLAND/RC045.1';
spikecd = 1;
[foo, index, raster, genparams,e1h,e1v,e1t] =...
             readdata(filename, spikecd);
indices_bincone;
monserial = genparams(GEN_MONSERIAL1)+...
            genparams(GEN_MONSERIAL2)*2^8+...
            genparams(GEN_MONSERIAL3)*2^16;
cal = getGlobalParams('cal',monserial);
fundamentals = getGlobalParams('fundamentals');
M = fundamentals*cal.P_device;
conestats = [];
for i = 2:3
   L = foo(:,I_COND) == i;

   tmp = [mean(foo(L,I_MEANL)); mean(foo(L,I_SDL));...
          mean(foo(L,I_MEANM)); mean(foo(L,I_SDM));...
          mean(foo(L,I_MEANS)); mean(foo(L,I_SDS))];

   bkgndrgb = [cal.gammaTable(mean(foo(L,I_BKGNDR))+1,1);
               cal.gammaTable(mean(foo(L,I_BKGNDG))+1,2);
               cal.gammaTable(mean(foo(L,I_BKGNDB))+1,3)];
   bkgndlms = M*bkgndrgb;
   tmp = [tmp; bkgndlms*1000];
   conestats = [conestats, tmp];
end
conestats = conestats./1000;
conestats(:,3) = conestats(:,1)./conestats(:,2);

L = foo(:,I_COND) == 2;
[nsp1, dt1] = countspikes(raster(L,:), index(L,:), SIN_STIMON, 0, SIN_STIMOFF, 0);
[nsp2, dt2] = countspikes(raster(~L,:), index(~L,:), SIN_STIMON, 0, SIN_STIMOFF, 0);
[h,p] = ttest2(nsp1./dt1, nsp2./dt2);
[p mean(nsp1./dt1)*1000 mean(nsp2./dt2)*1000]

maxT = 10;
STAs = [];
ns = [];
for i = 2:3
   L = logical(foo(:,I_COND) == i);
   [out] = getWhtnsStats5(foo(L,:),index(L,:),raster(L,:),genparams,maxT,0,'spikes',90,'STAmod');
   STA = out{1};
   STAs = [STAs, STA(:)];
   ns = [ns; out{3}];
end

STA1 = reshape(STAs(:,1),64,3,maxT);
STA2 = reshape(STAs(:,2),64,3,maxT);
diffSTA = STA1-STA2;

se = sqrt((1/ns(1)+(1/ns(2))));
zthresh = 4;
STAscalefactor = max(max(abs([STA1(:) STA2(:) diffSTA(:)])));
figure;
colormap(gray(255));
for i = 1:10
   tmp = reshape(STA1(:,:,i),8,8,3);
   z = tmp.*sqrt(ns(1));
   tmp = tmp./(2*STAscalefactor)+0.5;
   for j = 1:3
      subplot(9,maxT,maxT*(j-1)+i); hold on;
      image(squeeze(tmp(:,:,j))*255);
      axis image;
      set(gca,'visible','off')
      idxs = find(abs(z(:,:,j))>zthresh);
      if (any(idxs))
         [y,x] = ind2sub([8 8], idxs);
         plot(x,y,'m.');
      end
   end
   tmp = reshape(STA2(:,:,i),8,8,3);
   z = tmp.*sqrt(ns(2));
   tmp = tmp./(2*STAscalefactor)+0.5;
   for j = 4:6
      subplot(9,maxT,maxT*(j-1)+i); hold on;
      image(squeeze(tmp(:,:,j-3))*255);
      axis image;
      set(gca,'visible','off')
      idxs = find(abs(z(:,:,j-3))>zthresh);
      if (any(idxs))
         [y,x] = ind2sub([8 8], idxs);
         plot(x,y,'m.');
      end
   end
   tmp = reshape(diffSTA(:,:,i),8,8,3);
   z = tmp./se;
   tmp = tmp./(2*STAscalefactor)+0.5;
   for j = 7:9
      subplot(9,maxT,maxT*(j-1)+i); hold on;
      image(squeeze(tmp(:,:,j-6))*255);
      axis image;
      set(gca,'visible','off');
      idxs = find(abs(z(:,:,j-6))>zthresh);
      if (any(idxs))
         [y,x] = ind2sub([8 8], idxs);
         plot(x,y,'m.');
      end
   end
end

[whichpix, whichcone, whichframe] = pickapix(STA1+STA2,[], 1, [1 1 1]);
STA1_t = squeeze(STA1(whichpix,:,:))';
STA2_t = squeeze(STA2(whichpix,:,:))';

STA1_t_SDnorm = STA1_t./repmat(conestats([2 4 6],1)',size(STA1,3),1);
STA2_t_SDnorm = STA2_t./repmat(conestats([2 4 6],2)',size(STA2,3),1);
STA1_t_contrastnorm = STA1_t./(repmat(conestats([2 4 6],1)',size(STA1,3),1)./...
                        repmat(conestats([1 3 5],1)',size(STA1,3),1));
STA2_t_contrastnorm = STA2_t./(repmat(conestats([2 4 6],2)',size(STA2,3),1)./...
                        repmat(conestats([1 3 5],2)',size(STA1,3),1));

figure;
set(gcf,'DefaultAxesColorOrder',[1 0 0; 0 1 0; 0 0 1]);
subplot(3,3,1); 
plot(STA1_t);
title('Raw STA');
subplot(3,3,2);
plot(STA1_t_SDnorm);
title('SD normalized');
subplot(3,3,3); 
plot(STA1_t_contrastnorm);
title('Contrast normalized');
subplot(3,3,4);
plot(STA2_t);
subplot(3,3,5);
plot(STA2_t_SDnorm);
subplot(3,3,6);
plot(STA2_t_contrastnorm);
subplot(3,3,7);
plot(STA1_t-STA2_t);
subplot(3,3,8);
plot(STA1_t_SDnorm-STA2_t_SDnorm);
subplot(3,3,9);
plot(STA1_t_contrastnorm-STA2_t_contrastnorm);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Power analysis on STA differences
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
niter = 2000;
x1 = squeeze(STA1(whichpix,:,whichframe)).*ns(1);
x2 = squeeze(STA2(whichpix,:,whichframe)).*ns(2);
p1 = (x1./ns(1))/2+.5;
p2 = (x2./ns(2))/2+.5;

CIs = [];
pointest = STA1(whichpix,:,whichframe)./STA2(whichpix,:,whichframe);
for cone = 1:3
   stim1 = binornd(ns(1), p1(cone), 1, niter);
   stim2 = binornd(ns(2), p2(cone), 1, niter);
   xx1 = (stim1-(ns(1)-stim1))/ns(1);
   xx2 = (stim2-(ns(2)-stim2))/ns(2);
   ratios = xx1./xx2;
   CIs = [CIs; pointest(cone) prctile(ratios, [2.5 97.5])]
end

figure; axes; hold on; 
colors = [1 0 0; 0 1 0; 0 0 1];
for i = 1:3
   h = bar(i,CIs(i,1));
   set(h,'FaceColor',colors(i,:));
   plot([i i],[CIs(i,2),CIs(i,3)],'k-','LineWidth',2);
end
for i = 1:3
   h = bar(i+3,conestats(2*i-1,3))
   set(h,'FaceColor',colors(i,:));
end

set(gca,'XTick',[1 2 3 4 5 6],'XTickLabel',{'L wt ratio','M wt ratio','S wt ratio','L','M','S'})
title(filename)

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Trying to move to RGB
%%%%%%%%%%%%%%%%%%%%%%%%%%

monserial = genparams(GEN_MONSERIAL1)+...
            genparams(GEN_MONSERIAL2)*2^8+...
            genparams(GEN_MONSERIAL3)*2^16;
cal = getGlobalParams('cal',monserial);
fundamentals = getGlobalParams('fundamentals');
M = fundamentals*cal.P_device;
invM = inv(M);
newSTA1 = permute(STA1,[2 1 3]);
newSTA1 = reshape(newSTA1,3, size(newSTA1,2)*size(newSTA1,3));
means = mean([foo(L,I_MEANL) foo(L,I_MEANM) foo(L,I_MEANS)])/10000;
stds = mean([foo(L,I_SDL) foo(L,I_SDM) foo(L,I_SDS)])/10000;
newSTA1 = newSTA1./repmat(stds',1,size(newSTA1,2));
% Normally to go from cones to guns we go through invM
% so now we go through M'
rgbSTA1 = M'*newSTA1;
meanrgb = invM*means';
rgbSTA1 = rgbSTA1+repmat(meanrgb, 1, size(rgbSTA1,2));
rgbSTA1 = reshape(rgbSTA1,3,64,10);
rgbSTA1 = permute(rgbSTA1,[2 1 3]);

figure;
scalefactor = 1/(2*range(rgbSTA1(:)));
set(gcf,'DefaultAxesUnits','inches')
for i = 1:maxT
   tmp = squeeze(rgbSTA1(:,:,i));
   tmp = tmp-repmat(meanrgb', size(tmp,1), 1);
   tmp = tmp*scalefactor;
   tmp = tmp+repmat(meanrgb', size(tmp,1), 1);
   tmp = reshape(tmp,8,8,3);

   axes('position',[.6*(i-1) .5 .5 .5]);
   image(tmp)
   set(gca,'visible','off')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Looking at PCs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OPTIONS.issym = 1;
OPTIONS.disp = 1;

pcs = [];
STAs = [];
for whichframe = 1:10
   feval('STCOVmod4','init',[genparams(GEN_STIMNELEMX)*genparams(GEN_STIMNELEMY) 3 whichframe whichframe]);
   [out] = getWhtnsStats6(foo,index,raster,genparams,whichframe,0,'spike',90,'STCOVmod4');
   meanvect = out{1}';
   nstim = out{3};
   STCOV = out{2}/nstim;
   [v1,d1,flag] = eigs(STCOV,3,'LM',OPTIONS);
   pcs = [pcs, v1(:,1)];
   STAs = [STAs, out{1}];
end

STAscalefactor = max(abs(STAs(:)));
figure;
colormap(gray(255));
for i = 1:10
   tmp = reshape(STAs(:,i),3,8,8);
   tmp = tmp./(2*STAscalefactor)+0.5;
   for j = 1:3
      subplot(6,10,10*(j-1)+i)
      image(squeeze(tmp(j,:,:))*255);
      axis square;
      set(gca,'visible','off')
   end

   tmp = reshape(pcs(:,i),3,8,8);
   tmp = tmp./(2*max(abs(tmp(:))))+0.5;
   for j = 4:6
      subplot(6,10,10*(j-1)+i)
      image(squeeze(tmp(j-3,:,:))*255);
      axis square;
      set(gca,'visible','off')
   end
end

[whichpix, whichcone, whichframe] = pickapix(STA1+STA2,[],8,[1 1 1])
L = logical(foo(:,I_COND) == 2);
[j,i]= ind2sub([8 8],whichpix);

[SpikeStim1,TrialVect,TimeVect,Lspike] = getWhtnsStim1(foo(L,:),index(L,:),raster(L,:),genparams,maxT,0,'spike',90,[i j]);
for k = 2:max(Lspike)
   SpikeStim1 = cat(1,SpikeStim1, SpikeStim1(Lspike >= k,:,:,:));
end
[SpikeStim2,TrialVect,TimeVect,Lspike] = getWhtnsStim1(foo(~L,:),index(~L,:),raster(~L,:),genparams,maxT,0,'spike',90,[i j]);
for k = 2:max(Lspike)
   SpikeStim2 = cat(1,SpikeStim2, SpikeStim2(Lspike >= k,:,:,:));
end
[pcs1,score,d1] = princomp(vectorizestim(squeeze(SpikeStim1)));
[pcs2,score,d2] = princomp(vectorizestim(squeeze(SpikeStim2)));

figure; set(gcf,'DefaultAxesColorOrder',[1 0 0; 0 1 0; 0 0 1]);
subplot(2,3,1);
plot(squeeze(mean(SpikeStim1))');
subplot(2,3,2);
plot(reshape(pcs1(:,1),3,maxT)');
subplot(2,3,3);
plot(reshape(pcs1(:,2),3,maxT)');
subplot(2,3,4);
plot(squeeze(mean(SpikeStim2))');
subplot(2,3,5);
plot(reshape(pcs2(:,1),3,maxT)');
subplot(2,3,6);
plot(reshape(pcs2(:,2),3,maxT)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Trying to get a space/time/color
%%% RF for a complex cell by 1) calculating
%%% PC1s independently on a frame by frame
%%% basis and 2) flipping the signs 
%%% such the dot product between consecutive
%%% frames is positive. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filename = '/home/horwitz/Cortex/Data/KIMBALL/K069.3'
spikecd = 102;
[foo, index, raster, genparams] = readcort2(filename,spikecd);
paranum = get_paradigm(filename,1000);
if (paranum == 60)
   indices_whtnsSpatial2;
else
   indices_whtnsSpatial;
end
STAs = [];
PC1s = [];
for whichframe = 1:20
   feval('STCOVmod4','init',[genparams(GEN_STIMNELEMX)*genparams(GEN_STIMNELEMY) 3 whichframe whichframe]);
   [out] = getWhtnsStats3(foo,index,raster,genparams,whichframe,0,'spike',paranum,'STCOVmod4');
   STAs(whichframe,:) = out{1};
   [v,d] = eig(out{2});
   PC = v(:,end)*sqrt(d(end));
   if (whichframe == 1)
      PC1s(1,:) = PC';
   elseif(PC1s(whichframe-1,:)*PC < 0)
      PC1s(whichframe,:) = -PC';
   else
      PC1s(whichframe,:) = PC';
   end
end

figure;
for whichframe = 1:size(STAs,1)
   subplot(2,20,whichframe);
   im = reshape(STAs(whichframe,:),[3 8 8]);
   imboost(permute(im,[2 3 1]));
   subplot(2,20,whichframe+20);
   im = reshape(PC1s(whichframe,:),[3 8 8]);
   imboost(permute(im,[2 3 1]));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% Trying to calculate Pillowstat in the
%%% high dimensional space.
%%%
%%% This didn't work - even the most complex
%%% cells come out looking simple.  I'm
%%% not sure why.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

minT = 4;
maxT = 12;
datapath = '/home/horwitz/Cortex/Data';
monks = {'C','CORN';'K','KIMBALL';'T','TIGER';'D','DONDERS'};
[filelist, spikecodes] = readfilelist('/home/horwitz/Cortex/lists/whtnsSpatiallist4',monks);
path = '/home/horwitz/Cortex/Data';
data = {};
for whichcell = 1:length(filelist)
   % First, loading files
   filename = char(filelist{whichcell});
   spikecd = str2num(spikecodes{whichcell});
   disp([filename,': ',num2str(spikecd)]);
   paranum = get_paradigm([datapath,'/',filename],1000);
   [foo, index, raster, genparams] = readdata([path,'/',filename],spikecd);

   if (paranum == 60)
      indices_whtnsSpatial2;
   elseif (paranum == 31)
      indices_whtnsSpatial;
   else %  (paranum == 61)   % S part of U and S white noise
      indices_UandSwhtns;
      L = logical(foo(:,I_STIMTYPE) == 0);
      foo = foo(L,:); index = index(L,:); raster = raster(L,:);
      paranum = 62;
   end

   trials1 = [1:floor(size(foo,1)*.5)];
   trials2 = [trials1(end)+1:size(foo,1)];  % For cross-validation
   feval('STCOVmod4','init',[genparams(GEN_STIMNELEMX)*genparams(GEN_STIMNELEMY) 3 minT maxT]);
   [out] = getWhtnsStats3(foo(trials1,:),index(trials1,:),raster(trials1,:),genparams,maxT,0,'spike',paranum,'STCOVmod4');
   STA_spike = out{1};
   STCOV_spike = out{2};

   feval('STCOVmod4','init',[genparams(GEN_STIMNELEMX)*genparams(GEN_STIMNELEMY) 3 minT maxT]);
   [out] = getWhtnsStats3(foo(trials1,:),index(trials1,:),raster(trials1,:),genparams,maxT,0,'all',paranum,'STCOVmod4');
   STA_all = out{1};
   STCOV_all = out{2};

   [b, vals, GaussParams] = compiSTAC(STA_spike, STCOV_spike, STA_all, STCOV_all, 1);

   stimdurations = index(trials1,SIN_STIMOFF) - index(trials1,SIN_STIMON);
   STPROJmod('init',STA_spike,minT-1, sum(stimdurations)/10,[64 3 9]);
   out = getWhtnsStats3(foo(trials1,:),index(trials1,:),raster(trials1,:),genparams,maxT,0,'spike',paranum,'STPROJmod');
   proj = out{1};
   Lspike = out{2};

   %%% Getting rid of outliers
   lowerbound = prctile(proj,1);
   upperbound = prctile(proj,99);
   L = logical(proj < lowerbound | proj > upperbound);
   proj(L) = [];
   Lspike(L) = [];

   bins = linspace(min(proj), max(proj),15)';
   [Na,Xa] = hist(proj,bins);
   [Ns,Xs] = hist(proj(Lspike > 0),bins);
   fr = Ns./Na;

   [B_lin, BINT_lin, R, RINT, STATS_lin] = regress(fr', [ones(length(bins),1), bins]);
   [B_quad, BINT_quad, R, RINT, STATS_quad] = regress(fr', [ones(length(bins),1), bins.^2]);
   [B_both, BINT_both, R, RINT, STATS_both] = regress(fr', [ones(length(bins),1), bins, bins.^2]);

   Pillowstat = (STATS_quad(1)-STATS_lin(1))/STATS_both(1);

   data{whichcell}.MaxInfoVect = b;
   data{whichcell}.STA_spike = STA_spike;
   data{whichcell}.fr = fr;
   if (Pillowstat > 0 & B_quad(2) < 0)
      data{whichcell}.Pillowstat = nan;
   else
      data{whichcell}.Pillowstat = Pillowstat;
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Experimenting with using SVD
% to look at separability in color/
% space/time all together.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

name = '/home/horwitz/Cortex/Data/DONDERS/D031.1c5'
spikecd = 1;
maxT = 20;
[foo, index, raster, genparams,e1h,e1v,e1t] =...
             readdata(name, spikecd);
paranum = get_paradigm(name,1000);
if (paranum == 60)
   indices_whtnsSpatial2;
elseif (paranum == 31)
   indices_whtnsSpatial;
else %  (paranum == 61)   % S part of U and S white noise
   indices_UandSwhtns;
   L = logical(foo(:,I_STIMTYPE) == 0);
   foo = foo(L,:); index = index(L,:); raster = raster(L,:);
   paranum = 62;
end
         
[out] = getWhtnsStats(foo,index,raster,genparams,maxT,0,'spikes',paranum,'STAmod');
STA = out{1};
STV = out{2};
nspikes = out{3};

STAmat = reshape(permute(STA,[1 3 2]),64*maxT,3);
[u,s1,v] = svd(STAmat);
rgb_v = v(:,1);
predmat = v(:,1)*s(1,1)*u(:,1)';

spacetimemat = reshape(u(:,1),[64 maxT]);
[u,s2,v] = svd(spacetimemat);
predmat = v(:,1)*s(1,1)*u(:,1)';
time_v = v(:,1);
space_v = u(:,1);

spacergbmat = space_v*rgb_v';
timeweighting = permute(repmat(time_v,[1 64 3]),[2 3 1]);
fullpredvol1 = repmat(spacergbmat,[1 1 maxT]).*timeweighting.*s1(1,1).*s2(1,1);

%%% This way seems to have worked well.  Now let's see what happens
%%% when we do the separating in a different order (peeling off time first)

STAmat = reshape(STA,64*3,maxT);
[u,s1,v] = svd(STAmat);
time_v = v(:,1);
spacecolormat = reshape(u(:,1),[64 3]);

[u,s2,v] = svd(spacecolormat);
space_v = u(:,1);
rgb_v = v(:,1);

spacergbmat = space_v*rgb_v';
timeweighting = permute(repmat(time_v,[1 64 3]),[2 3 1]);
fullpredvol2 = repmat(spacergbmat,[1 1 maxT]).*timeweighting.*s1(1,1).*s2(1,1);

corrcoef([fullpredvol1(:) fullpredvol2(:)])

figure; axes; hold on;
plot(STA(:),'b-');
plot(fullpredvol1(:),'g-');
plot(fullpredvol2(:),'m-');

%%% These two ways of separating the STA into space, time, and color 
%%% (either peeling off color first, or time first) resulted in slightly
%%% different results, but they're pretty close.  Peeling off time first
%%% results in slightly cleaner estimates for some cells but not others.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Seeing how cone weight estimates (or RGB weight estimates) change
%%% as we incorporate more pixels into the the mix.  This might provide
%%% a more intuitive way of describing the inseparability.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

monserial = genparams(GEN_MONSERIAL1)+...
            genparams(GEN_MONSERIAL2)*2^8+...
            genparams(GEN_MONSERIAL3)*2^16;
cal = getGlobalParams('cal',monserial);
fundamentals = getGlobalParams('fundamentals');
M = fundamentals*cal.P_device;

[whichpix,whichgun,whichframe] = pickapix(STA,STV,nspikes,1);
templates = logical([]);
tmp = zeros(8);
[i,j] = ind2sub([8 8],whichpix);
tmp(i,j) = 1;
templates(:,1) = tmp(:);

tmp = zeros(8);
idxs = fullfact([3 3])-2+repmat([i j],[9 1]);
for k = 1:size(idxs,1)
   tmp(idxs(k,1),idxs(k,2)) = 1;
end
templates = [templates, tmp(:)];

tmp = zeros(8);
idxs = fullfact([5 5])-3+repmat([i j],[25 1]);
for k = 1:size(idxs,1)
   tmp(idxs(k,1),idxs(k,2)) = 1;
end
templates = [templates, tmp(:)];

%%% Getting RGB values the "mean" way
rgbs = [];
for i = 1:size(templates,2)
   tmpSTA = squeeze(STA(find(templates(:,i)),:,:));
   if (ndims(tmpSTA) > 2)
      tmpSTA = squeeze(mean(tmpSTA));
   end
   [u,s,v] = svd(tmpSTA);
   if (max(abs(tmpSTA'*u(:,1))) == max(tmpSTA'*u(:,1)))
      rgbs(:,i) = u(:,1);
   else
      rgbs(:,i) = -u(:,1);
   end
end
plot(rgbs','-')

%%% Getting RGB values the SVD way
rgbs = [];
for i = 1:size(templates,2)
   tmpSTA = squeeze(STA(find(templates(:,i)),:,:));
   if (ndims(tmpSTA) > 2)
      tmpSTAmat = reshape(permute(tmpSTA,[2 1 3]),size(tmpSTA,2),size(tmpSTA,1)*size(tmpSTA,3));
   else
      tmpSTAmat = tmpSTA;
   end
   [u,s,v] = svd(tmpSTAmat);
   if (max(abs(tmpSTAmat'*u(:,1))) == max(tmpSTAmat'*u(:,1)))
      rgbs(:,i) = u(:,1);
   else
      rgbs(:,i) = -u(:,1);
   end
end
plot(rgbs',':')

coneweights = inv(M')*rgbs;
normconeweights = coneweights./repmat(sum(abs(coneweights)),3,1)
plot(normconeweights(1,:),normconeweights(2,:),'k.-')
set(gca,'Xlim',[-1 1],'Ylim',[-1 1]);
axis square;

% Fractional change in normalized cone weights
min([normconeweights(:,1)'./normconeweights(:,3)'; normconeweights(:,3)'./normconeweights(:,1)'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Asking the question, does a combination of measurement error 
%%% and the statistical properties of the cone mosaic predict
%%% the lack of color/space separability?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% From PopAnal2, section 4.3
load /home/horwitz/matlab/PopAnal/Data/CenSurr.mat

% From Simulations.m, section 24
conemosaic_mu = [40.8357 28.4899 6.5187]; % In number of cones
conemosaic_sigma = [ 14.4962 -14.1756 0.7562;...
                    -14.1756 23.8179 -1.2000;...
                      0.7562 -1.2000 1.3488];


whichcell = 13;
datapath = '/home/horwitz/Cortex/Data/';
filename = data{whichcell}.filename;
spikecd = data{whichcell}.spikecd;
disp([filename,': ',num2str(spikecd)]);
[foo, index, raster, genparams] = readdata([datapath,'/',filename],spikecd);
paranum = get_paradigm([datapath,'/',filename],1000);
if (paranum == 60)
   indices_whtnsSpatial2;
elseif (paranum == 31)
   indices_whtnsSpatial;
else  % paranum == 61
   indices_UandSwhtns;
   L = logical(foo(:,I_STIMTYPE) == 0);
   foo = foo(L,:);  index = index(L,:);  raster = raster(L,:);
   paranum = 62;
end

monserial = genparams(GEN_MONSERIAL1)+...
            genparams(GEN_MONSERIAL2)*2^8+...
            genparams(GEN_MONSERIAL3)*2^16;
cal = getGlobalParams('cal',monserial);
fundamentals = getGlobalParams('fundamentals');
M = fundamentals*cal.P_device;

feval('STCOVmod4','init',[genparams(GEN_STIMNELEMX)*genparams(GEN_STIMNELEMY) 3 data{whichcell}.whichframe data{whichcell}.whichframe]);
[out] = getWhtnsStats3(foo,index,raster,genparams,data{whichcell}.whichframe,0,'spike',paranum,'STCOVmod4');

STA = out{1};
STC = out{2}/out{3};  % Already divided by 'n'
n = out{3};

tmpSTA = [];
tmpSTC = [];
for pixidx = ([1:64]-1)*3+1
   tmpSTA = cat(2,tmpSTA,inv(M')*STA(pixidx:pixidx+2));
   tmpSTC = cat(3,tmpSTC,inv(M')*STC([pixidx:pixidx+2],[pixidx:pixidx+2])*inv(M));
end

% Looking at the center pixel
lms_mu = tmpSTA(:,data{whichcell}.centerpix);
lms_S = squeeze(tmpSTC(:,:,data{whichcell}.centerpix));

% Assuming the model coneweight = conecount*unitconeweight
weightpercone_mu = lms_mu./conemosaic_mu';
weightpercone_S = weightpercone_mu*weightpercone_mu'.*conemosaic_sigma;
% The above step is not eniterly accurate.  To do this right I 
% would have to figure out what the additional noise introduced
% by lms_mu is (since it's random).  But to do so would be very complicated
% since weightpercone_S would have to vary continuously throughout
% the lms_S-defined confidence ellipse.  Yuck.

% Interesting - it looks like the the measurement error is on
% the order of the noise in the cone mosaic (for cell #1 anyway).
grand_S = weightpercone_S+lms_S;

%%% Plotting axes
figure; axes; hold on;
crit = finv(1-0.01,3,n-3);
crit = crit*((n-1)*3)/(n-3);
x = [0:.01:2*pi+.01];
xx = [cos(x); sin(x); zeros(1,length(x))];
orders = [1 2 3;2 3 1;3 1 2];
tmpS = cat(3,lms_S, weightpercone_S, grand_S);
colors = {'blue','magenta','black'};
for j = 1:size(tmpS,3)
   [v,d] = eig(squeeze(tmpS(:,:,j)));
   newd = diag(sqrt(diag(d)*crit));
   for i = 1:3
      tmp = (v*newd)*xx(orders(i,:),:);
      h1 = plot3(tmp(1,:)+lms_mu(1),tmp(2,:)+lms_mu(2),tmp(3,:)+lms_mu(3),'m-');
      h2 = plot3(tmp(1,:)-lms_mu(1),tmp(2,:)-lms_mu(2),tmp(3,:)-lms_mu(3),'m-');
      set(h1,'LineWidth',2);
      set(h1,'Color',colors{j});
      set(h2,'LineWidth',2);
      set(h2,'Color',colors{j});
   end
end

Lcen = find(data{whichcell}.Center);
Lsurr = find(data{whichcell}.Surround);
plot3(tmpSTA(1,Lcen),tmpSTA(2,Lcen),tmpSTA(3,Lcen),'r.');
plot3(tmpSTA(1,Lsurr),tmpSTA(2,Lsurr),tmpSTA(3,Lsurr),'k.');

% Plotting weghts (all the same norm as the center)
%normweights = tmpSTA./repmat(sqrt(sum(tmpSTA.^2)),3,1);
%normweights = normweights.*repmat(norm(tmpSTA(:,data{whichcell}.centerpix)),3,size(tmpSTA,2));
%for i = Lcen'
%   plot3([0 normweights(1,i)],[0 normweights(2,i)],[0 normweights(3,i)],'r-');
%end
%for i = Lsurr'
%   plot3([0 normweights(1,i)],[0 normweights(2,i)],[0 normweights(3,i)],'k-');
%end

set(gca,'Xlim',[-1 1],'Ylim',[-1 1],'Zlim',[-1 1]);
plot3([-1 1],[0 0],[0 0],'k-');
plot3([0 0],[-1 1],[0 0],'k-');
plot3([0 0],[0 0],[-1 1],'k-');

STArgbmat = reshape(STA,3,64); % Going back to the original RGB representation

% Using the center pixel instead of SVD
[u,s,v] = svd(STArgbmat);
u(:,1) = STArgbmat(:,data{whichcell}.centerpix);
u_lms = gsorth(inv(M')*u);

%%% Question: Does there exist a line that passes through all of these ellipsoids?  
%%% Not sure how to address this.  Here's something easier...

%%% Rotating the space so that u_lms(:,1) in the y-axis.
%%% I think I can just use u_lms' as the rotation matrix
%%% Then I'm reduced to asking "do all the ellipses in this plane
%%% contain the origin?"

%%% Plotting axes
figure; 
subplot(2,2,1); hold on;
subplot(2,2,2); hold on;
% Critical value from Zar 
crit = finv(1-0.01,2,2*n-2-1);
crit = crit*((2*n-2)*2)/(2*n-2-1);
x = [0:.01:2*pi+.01];
xx = [cos(x); sin(x)];
L = nan*ones(2,size(tmpSTA,2));
for i = 1:2  % Which error to use (measurement or measurement + cones)
   subplot(2,2,i);
   for j = 1:size(tmpSTA,2)
      projection = abs(tmpSTA(:,j)'*tmpSTA(:,data{whichcell}.centerpix))./norm(tmpSTA(:,data{whichcell}.centerpix));
      normfactor = norm(tmpSTA(:,data{whichcell}.centerpix))/projection;
      if (i == 1)
         S = squeeze(tmpSTC(:,:,j)+tmpSTC(:,:,data{whichcell}.centerpix)/normfactor.^2);
         rotSTC = u_lms'*S*u_lms;
      else
         S = squeeze(tmpSTC(:,:,j)+tmpSTC(:,:,data{whichcell}.centerpix)/normfactor.^2)+weightpercone_S;
         rotSTC = u_lms'*S*u_lms;
      end
      rotSTA = u_lms'*tmpSTA(:,j);
      % This first coordinate is the position along the u_lms(:,1) axis
      % which we don't care about.  (Except for +/-).
      polarity = sign(rotSTA(1));
      rotSTA = rotSTA([2 3]);
      rotSTC = rotSTC([2 3],[2 3]);
      [v,d] = eig(rotSTC);
      newd = diag(sqrt(diag(d)*crit));
      tmp = (v*newd)*xx;
      h = plot(rotSTA(1)+tmp(1,:),rotSTA(2)+tmp(2,:),'r-','Linewidth',2);
      %h = plot(rotSTA(1),rotSTA(2),'r*');
      if (polarity == -1)  % Part of the "surround"
         set(h,'Color','black');
      end
   
      % Does the ellipse contain the origin?
      % This method is prone to numerical errors
      tmpvect = inv(v*newd)*rotSTA;
      tmpvect = tmpvect./norm(tmpvect);
      greg = (v*newd)*tmpvect;
      mahal = (rotSTA'*inv(rotSTC)*rotSTA);
      if (norm(greg) < norm(rotSTA))
         if (mahal < crit & j ~= data{whichcell}.centerpix)
            error('Disagreement between ways of calculating Mahal');
         end
         L(i,j) = 1;
         set(h,'Linewidth',3);
      else
         if (mahal > crit & j ~= data{whichcell}.centerpix)
            error('Disagreement between ways of calculating Mahal');
         end
         L(i,j) = 0;
         set(h,'Linewidth',.5,'Linestyle',':');
      end  
   end
   plot([-1 1],[0 0],'k-');
   plot([0 0],[-1 1],'k-');
   axis square;
end
for i = 1:2
   subplot(2,2,2+i); hold on;
   imboost(permute(reshape(STA,3,8,8),[2 3 1]));
   [row,col] = ind2sub([8 8],find(L(i,:)));
   plot(col,row,'k*');
   axis image; axis ij;
end