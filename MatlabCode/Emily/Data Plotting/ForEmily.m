%Long since replaced by iterateAndPlotFiles_modularPlusDB.
% Here's some code for analyzing LMTF data for Emily
% Section 1) Plotting the data from an LMTF experiment and the fit
% Section 2) Eye position analyses
% Section 3) graphical depiction of the yoked tilted rampy trough model
%%
% Section 1
% Read the text file that contains the names of all the .nex files.
flist = flatten(fnamesFromTxt2(fullfile(nexfilepath,'nexfilelists','Emily', 'ApolloTest.txt')));
startpath = '/Volumes/NO BACKUP/NexFiles/Greg/Apollo';
data = [];
% Iterate over the list of .nex files.
for i = 1:length(flist)
    stro = notnex2stro(findfile(flist{i},startpath)); % <-- process the information in the nex file and put it into the "stro" structure
    % Below, just figuring out what information is in what column
    Lstim_idx = strcmp(stro.sum.trialFields(1,:), 'stim_idx');
    Llcc = strcmp(stro.sum.trialFields(1,:), 'lcc');
    Lmcc = strcmp(stro.sum.trialFields(1,:), 'mcc');
    Ltf = strcmp(stro.sum.trialFields(1,:), 'tf');
    %Loog = strcmp(stro.sum.trialFields(1,:), 'oog');
    
    % Getting the threshold points
    [stim_idxs,init_stim_trial_idxs] = unique(stro.trial(:,Lstim_idx),'last');
    questmodes = stro.trial(init_stim_trial_idxs, Llcc|Lmcc);
    tfs = stro.trial(init_stim_trial_idxs,Ltf);
    
    % Out of gamut checking
    funds = reshape(stro.sum.exptParams.fundamentals,length(stro.sum.exptParams.fundamentals)/3,3);
    if (size(stro.sum.exptParams.mon_spd,1) == 303)
        spds = SplineSpd([380:4:780]',reshape(stro.sum.exptParams.mon_spd,length(stro.sum.exptParams.mon_spd)/3,3),[380:5:780]');    
    else
        spds = SplineSpd([380:2:780]',reshape(stro.sum.exptParams.mon_spd,length(stro.sum.exptParams.mon_spd)/3,3),[380:5:780]');
    end
    M = funds'*spds;
    bkgndrgb = stro.sum.exptParams.bkgndrgb;
    [in_gamut,scalar] = gamutCheck([questmodes zeros(size(questmodes,1),1)]', bkgndrgb, M, 'both');
    questmodes(~in_gamut) = questmodes(~in_gamut).*scalar(~in_gamut);
    data = [data; questmodes tfs ~in_gamut'];
end
Loog = logical(data(:,end));

% At this point "data" has n rows and 4 columns. The n rows correspond to n
% different stimulus directions (L,M,TF combinations) that were tested. The
% four columns are (1) L-cone contrast, (2) M-cone contrast, (3) Temporal
% frequency, (4) out of gamut (1 or 0)
% -----------------------------------------
% Now we're transitioning from collecting the data from the various files
% to fitting the data with the surface model.

% Stuff we're going to need later
[th,r] = cart2pol(data(:,1),data(:,2)); % Remember, the first two columns of data are L- and M-cone contrasts
% This line is just re-representing them in polar coordinates.
tf = data(:,3);
LB = [0 0 1 1 .001 .001];
UB = [100 1 20 20 .03 .03];
% Make the filter into a function with a handle
% This "@" syntax is used to define a little function without having to
% defin eit in its own .m file. It only works with short functions that can
% be written on one line.
%f1 = @(fpar,omega)fpar(1)*abs(((1i*2*pi*fpar(5).*omega+1).^-fpar(3))-fpar(2)*((1i*2*pi*fpar(6).*omega+1).^-fpar(4)));
%f2 = @(fpar,omega)fpar(1+6)*abs(((1i*2*pi*fpar(5+6).*omega+1).^-fpar(3+6))-fpar(2+6)*((1i*2*pi*fpar(6+6).*omega+1).^-fpar(4+6)));

x = data(:,1);
y = data(:,2);

initparams = [40 .1 9 3 .005 .002]; % An initial guess for the parameter values
options = optimset('MaxFunEvals',50000,'MaxIter',50000,'TolFun',10^-8);
thetas = linspace(0,pi/2,12);
toterr = [];
fpars = [];
fpar = [initparams, initparams];
% Fitting the data many times (each time rotating the data by a small
% amount)
for i = 1:length(thetas)
    rotmat = [cos(thetas(i)) -sin(thetas(i)); sin(thetas(i)) cos(thetas(i))];
    % Rotating data clockwise = rotating fit axis counterclockwise.
    xytmp = [x,y]*rotmat;
    [fpar,fv] = fmincon(@(params) tf_fiterr2(params,[xytmp(:,1) xytmp(:,2) tf],Loog),fpar,...
        [],[],[],[],[LB LB],[UB UB],[],options);
    toterr(i) = fv;
    fpars(i,:) = fpar;
end
bestrotidx = find(toterr == min(toterr));
reshape(fpars(bestrotidx,:),[6 2])  % In case the user wants to see the fitted parameters
fpar = fpars(bestrotidx,:);

% These two lines below look a lot like two lines that appear above. In
% fact I think they might be exactly the same, but "fpar" has been
% redefined so the definitions of f1 and f2 are a bit different now. In
% particular, they are fits to one of the several rotations of the data
% that were tried.
f1 = @(omega)fpar(1)*abs(((1i*2*pi*fpar(5).*omega+1).^-fpar(3))-fpar(2)*((1i*2*pi*fpar(6).*omega+1).^-fpar(4)));
f2 = @(omega)fpar(1+6)*abs(((1i*2*pi*fpar(5+6).*omega+1).^-fpar(3+6))-fpar(2+6)*((1i*2*pi*fpar(6+6).*omega+1).^-fpar(4+6)));
if (f1(1) < f2(1)) % "If f1 is luminance" then exchange f1 and f2, forcing f1 to be chromatic
    fpar= fpar([[7:12]';[1:6]']);
    f1 = @(omega)fpar(1)*abs(((1i*2*pi*fpar(5).*omega+1).^-fpar(3))-fpar(2)*((1i*2*pi*fpar(6).*omega+1).^-fpar(4)));
    f2 = @(omega)fpar(1+6)*abs(((1i*2*pi*fpar(5+6).*omega+1).^-fpar(3+6))-fpar(2+6)*((1i*2*pi*fpar(6+6).*omega+1).^-fpar(4+6)));
end

% All the plotting stuff is below
% Plotting the raw data 
figure; axes; hold on; 
plot3(x(~Loog),y(~Loog),tf(~Loog),'ko','MarkerFaceColor','black','MarkerSize',5);
plot3(-x(~Loog,1),-y(~Loog),tf(~Loog),'ko','MarkerFaceColor','black','MarkerSize',5);
plot3(x(Loog,1),y(Loog),tf(Loog),'ro','MarkerFaceColor','red','MarkerSize',5);
plot3(-x(Loog,1),-y(Loog),tf(Loog),'ro','MarkerFaceColor','red','MarkerSize',5);

% Adjusting the axes so we can see everything
lim = max(max(abs([x y])));
set(gca,'Xlim',1.1*[-lim lim]);
set(gca,'Zscale','log');
set(gca,'Ylim',1.1*[-lim lim]);
axis square
xlabel('L-cone contrast');
ylabel('M-cone contrast');
zlabel('TF (Hz)');
set(gca,'View',[135 20]);
axis vis3d

% Plotting the fit
[xx,yy,zz] = meshgrid(linspace(-max(abs(x)),max(abs(x)),20),...
    linspace(-max(abs(y)),max(abs(y)),20),...
    linspace(min(tf),max(tf),20));
a = abs(f1(zz)).^-1; % chromatic threshold
b = abs(f2(zz)).^-1; % luminance threshold
thtmp = atan2(yy,xx)+thetas(bestrotidx); % clockwise rotation from [L,M] to [a,b]
rtmp = (a.*b)./sqrt((b.*cos(thtmp)).^2+(a.*sin(thtmp)).^2); % radius of ellipse - thank you, Wikipedia

V = sqrt(xx.^2+yy.^2)-rtmp;
FV = isosurface(xx,yy,zz,V,0);
h = patch(FV);
%set(h,'FaceColor','green','EdgeColor','none');
set(h,'FaceColor','green','FaceAlpha',.5','EdgeColor','none','EdgeAlpha',0);
set(gcf,'Renderer','painters');
% Good views 
set(gca,'View',[135 12]);
set(gca,'View',[225 12]);

%%
% Section 2) some eye position analysis
stro = nex2stro(findfile('E091916001.nex'));
ntrials = length(stro.trial);
fpacq_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fpacq_t'));
fpoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fpoff_t'));
targon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'targon_t'));
saccmade_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'saccmade_t'));
act_stim_x = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'act_stim_x'));
tf = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'tf'));
stim_idx = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_idx'));
correct = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'correct'));
stim_idx = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_idx'));

plotcolors = jet(max(stim_idx));
figure;
for i = 1:ntrials
    x = stro.ras{i,strcmp(stro.sum.rasterCells,'AD11')}*4096/400;
    y = stro.ras{i,strcmp(stro.sum.rasterCells,'AD12')}*4096/400;
    eyetimes = [0:length(x)-1]./stro.sum.analog.storeRates{1} + stro.ras{i,strcmp(stro.sum.rasterCells,'anlgStartTime')};
    fixL = logical(eyetimes >= fpoff_t(i) & eyetimes < saccmade_t(i)-.1);
    meanx = mean(x(fixL));
    meany = mean(y(fixL));
    sacL = logical(eyetimes >= fpoff_t(i) & eyetimes < saccmade_t(i)+.2);
    subplot(2,2,1); hold on;
    h(1) = plot(eyetimes(sacL)-fpoff_t(i), x(sacL)-meanx,'r');
    subplot(2,2,2); hold on;
    h(2) = plot(eyetimes(sacL)-fpoff_t(i), y(sacL)-meany,'r');
    subplot(2,2,3); hold on;
    h(3) = plot(x(sacL)-meanx,y(sacL)-meany,'r');

    set(h,'Color',plotcolors(stim_idx(i),:));
end
subplot(2,2,1);
xlim = get(gca,'xlim'); ylim = get(gca,'ylim');
xlabel('Time (s) from fpoff'); ylabel('DVA'); title('Horizontal');

subplot(2,2,2);
set(gca,'Xlim',xlim,'Ylim',ylim);
xlabel('Time (s) from fpoff'); ylabel('DVA'); title('Vertical');

subplot(2,2,3);
set(gca,'Xlim',[-5 5],'Ylim',[-5 5]);
xlabel('DVA'); ylabel('DVA');

%%
% Section 3
% Surface plots for the new rampy trough model
load(fullfile(fileparts(which('IsoSampOnline.m')), 'private', 'data', 'LMTF.mat'))
params = A.legacy.mode4params; % For example
xi_lum_coeffs = params(12:15); % I'm taking these from LMTF.mat but they could be fit by regression
xi_rg_coeffs = params(16:19);

[r,phi] = meshgrid(linspace(2,15,40),linspace(-pi/2,pi/2,40));
[x,y] = pol2cart(phi,r);
X = [ones(numel(phi),1),r(:),r(:).*sin(2*phi(:)),r(:).*cos(2*phi(:))]; % <---- New tilted rampy trough
pred = X*xi_lum_coeffs;

figure; colormap(jet)
subplot(2,1,1);
surf(phi,r,reshape(pred,size(phi)),ones(size(r))); 
xlabel('\phi');
ylabel('r');
zlabel('log contrast sensitivity');
set(gca,'View',[133 15]);
subplot(2,1,2);
h = surf(x,y,reshape(pred,size(phi)),ones(size(r)));
set(h,'EdgeColor','none')
axis equal;
set(gca,'View',[0 90]);
set(h,'CDataMode','auto')
xlabel('horizontal');
ylabel('vertical');



