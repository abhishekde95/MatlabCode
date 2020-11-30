% This script is for generating figure for a Vision Research paper



%%
% Pre-Defined Variables
thetaspace = pi/4;
rhospace = .5;
gl.nRnds = 3;
gl.nPres = 5;
gl.lummaxcc = .9;
gl.colmaxcc = .09;
gl.ang = pi;
params0 = [20 .1 100 2.1 3 gl.ang];

% Construct Polar Grid
if mod(gl.nRnds,2) == 0
    thetaspace = thetaspace / 2^((gl.nRnds-2)/2);
    rhospace = rhospace / 2^(gl.nRnds/2);
elseif mod(gl.nRnds,2) == 1
    thetaspace = thetaspace / 2^((gl.nRnds-1)/2);
    rhospace = rhospace / 2^((gl.nRnds-1)/2);
end
thetas = shiftdim(0:thetaspace:2*pi-thetaspace);
rhos = shiftdim(rhospace:rhospace:1);

% Ennumerate all conditions
PolRhoIdx = fullfact([numel(thetas) numel(rhos)]);
rhothetalist = repmat([thetas(PolRhoIdx(:,1)) rhos(PolRhoIdx(:,2))],gl.nPres,1);
thetas = rhothetalist(:,1);
rhos = rhothetalist(:,2);

% Transform from polar to cartesian coordinates
[tempLcc,tempMcc] = pol2cart(thetas,rhos);

% Build radially sytmmetric space
sdLcc = sqrt(var(tempLcc));
sdMcc = sqrt(var(tempMcc));
points_radsym = [tempLcc./sdLcc tempMcc./sdMcc];

% Scale and Rot matrix (radsym space to cov space)
rotMat = [cos(-pi/4) sin(-pi/4); -sin(-pi/4) cos(-pi/4)];
scaleMat = [gl.lummaxcc*sdLcc 0; 0 gl.colmaxcc*sdMcc];
covMat = rotMat * scaleMat * inv(rotMat);
points_cov = points_radsym * covMat;

% Nonlinear scaling (radsym space to LM Space)
scale = gl.lummaxcc*gl.colmaxcc./sqrt((gl.colmaxcc.*cos(thetas-pi/4)).^2 ...
    +(gl.lummaxcc.*sin(thetas-pi/4)).^2);
Lcc = tempLcc .* scale;
Mcc = tempMcc .* scale;
points_LM = [Lcc Mcc];

% Whitening matrix (LM space to white space)
M = inv(sqrtm(inv(cov(points_LM))));
points_white = points_LM * inv(M);

% correcting true angle in radsym space
trueLM_x = cos(gl.ang);
trueLM_y = sin(gl.ang);
truewhite = [trueLM_x trueLM_y] * M';
[trueang_white,~] = cart2pol(truewhite(1),truewhite(2));
truecov_x = cos(gl.ang);
truecov_y = sin(gl.ang);
trueradsym = [truecov_x truecov_y] * covMat';
[trueang_radsym,~] = cart2pol(trueradsym(1),trueradsym(2));
deltaang = trueang_white-trueang_radsym;
tweakMat = [cos(deltaang) sin(deltaang); -sin(deltaang) cos(deltaang)];
covMat = inv(tweakMat) * covMat;

% Adjust radsym points
points_radsym = points_cov * inv(covMat);

% Responses to stimuli (must do each space seperately bc of nonlinear
% transoform between radsym and white)
nsp0 = ComputeNakaRushtonJPW(params0,points_cov,'surface7');
nsp0 = poissrnd(nsp0);
nsp1 = ComputeNakaRushtonJPW(params0,points_LM,'surface7');
nsp1 = poissrnd(nsp1);


% Spike-Triggered Average Stimulus
STA_radsym = sum(points_radsym .* repmat(nsp0,1,2),1)./sum(nsp0.*2);
STA_cov = sum(points_cov .* repmat(nsp0,1,2),1)./sum(nsp0.*2);
STA_white = sum(points_white .* repmat(nsp1,1,2),1)./sum(nsp1.*2); 
STA_LM = sum(points_LM .* repmat(nsp1,1,2),1)./sum(nsp1.*2);  


% Fitting Model
% Set up some variables
angs = linspace(-pi,pi,17);
GOF = nan(numel(angs),1);
vub = [max(nsp1)   50   50   5     max(nsp1)    pi];
vlb = [.01       .01  .01   .5      0        -pi];
options = optimset('Algorithm','interior-point','MaxFunEvals',5000,...
    'MaxIter',5000,'FinDiffType','central','Hessian','lbfgs',...
    'DiffMinChange',.05,'display','off');
Aguess = max(nsp1);
expguess = 3;
blguess = min(nsp1)+1;
params = nan(numel(angs),numel(vlb));
A = [0 1 -1 0 0 0];
b = 0; 

% Fitting Model in radially symmetric Space
for rot = 1:numel(angs)
    rotMat = [cos(angs(rot)) sin(angs(rot)); -sin(angs(rot)) cos(angs(rot))];
    projpts = points_radsym * rotMat;
    sigguess = max(projpts(:,1))/2;
    paramsGuess = [Aguess sigguess 50 expguess blguess angs(rot)];
    [f1,fval] = fmincon('FitNakaRushtonFunJPW',paramsGuess,A,b,[],[],vlb,vub,[],options,points_radsym,nsp0,'surface7','Poisson');
    params(rot,:) = f1;
    GOF(rot) = fval;
end
[~,bestIdx] = min(GOF);
params_radsym = params(bestIdx,:);
ang = params_radsym(end);
rotMat = [cos(ang) sin(ang); -sin(ang) cos(ang)];
projpts = points_radsym * rotMat;
rho = max(projpts(:,1));
[modelfit_radsym(1),modelfit_radsym(2)] = pol2cart(ang,rho);

% Fitting Model in Cov Space
for rot = 1:numel(angs)
    rotMat = [cos(angs(rot)) sin(angs(rot)); -sin(angs(rot)) cos(angs(rot))];
    projpts = points_cov * rotMat;
    sigguess = max(projpts(:,1))/2;
    paramsGuess = [Aguess sigguess 50 expguess blguess angs(rot)];
    [f1,fval] = fmincon('FitNakaRushtonFunJPW',paramsGuess,A,b,[],[],vlb,vub,[],options,points_cov,nsp0,'surface7','Poisson');
    params(rot,:) = f1;
    GOF(rot) = fval;
end
[~,bestIdx] = min(GOF);
params_cov = params(bestIdx,:);
ang = params_cov(end);
rotMat = [cos(ang) sin(ang); -sin(ang) cos(ang)];
projpts = points_cov * rotMat;
rho = max(projpts(:,1));
[modelfit_cov(1),modelfit_cov(2)] = pol2cart(ang,rho);

% Fitting model in White space
for rot = 1:numel(angs)
    rotMat = [cos(angs(rot)) sin(angs(rot)); -sin(angs(rot)) cos(angs(rot))];
    projpts = points_white * rotMat;
    sigguess = max(projpts(:,1))/2;
    paramsGuess = [Aguess sigguess 50 expguess blguess angs(rot)];
    [f1,fval] = fmincon('FitNakaRushtonFunJPW',paramsGuess,A,b,[],[],vlb,vub,[],options,points_white,nsp1,'surface7','Poisson');
    params(rot,:) = f1;
    GOF(rot) = fval;
end
[~,bestIdx] = min(GOF);
params_white = params(bestIdx,:);
ang = params_white(end);
rotMat = [cos(ang) sin(ang); -sin(ang) cos(ang)];
projpts = points_white * rotMat;
rho = max(projpts(:,1));
[modelfit_white(1),modelfit_white(2)] = pol2cart(ang,rho);

% Fitting Model in LM Space
for rot = 1:numel(angs)
    rotMat = [cos(angs(rot)) sin(angs(rot)); -sin(angs(rot)) cos(angs(rot))];
    projpts = points_LM * rotMat;
    sigguess = max(projpts(:,1))/2;
    paramsGuess = [Aguess sigguess 50 expguess blguess angs(rot)];
    [f1,fval] = fmincon('FitNakaRushtonFunJPW',paramsGuess,A,b,[],[],vlb,vub,[],options,points_LM,nsp1,'surface7','Poisson');
    params(rot,:) = f1;
    GOF(rot) = fval;
end
[~,bestIdx] = min(GOF);
params_LM = params(bestIdx,:);
ang = params_LM(end);
rotMat = [cos(ang) sin(ang); -sin(ang) cos(ang)];
projpts = points_LM * rotMat;
rho = max(projpts(:,1));
[modelfit_LM(1),modelfit_LM(2)] = pol2cart(ang,rho);


% Set Up Some Figure Stuff

% Surface for cov space
x = linspace(min(points_cov(:,1)),max(points_cov(:,1)),100);
y = linspace(min(points_cov(:,2)),max(points_cov(:,2)),100);
[X_cov,Y_cov] = meshgrid(x,y);
surf = ComputeNakaRushtonJPW(params0,[X_cov(:) Y_cov(:)],'surface7');
surf_cov = reshape(surf,size(X_cov));

% Surface for rad sym space
temp = [X_cov(:) Y_cov(:)] * inv(M);
X_radsym = reshape(temp(:,1),size(X_cov));
Y_radsym = reshape(temp(:,2),size(Y_cov));
surf_radsym = surf_cov;
% x = linspace(min(points_radsym(:,1)),max(points_radsym(:,1)),100);
% y = linspace(min(points_radsym(:,2)),max(points_radsym(:,2)),100);
% [X_radsym,Y_radsym] = meshgrid(x,y);
% surf = ComputeNakaRushtonJPW(params0,[X_radsym(:) Y_radsym(:)],'surface7');
% surf_radsym = reshape(surf,size(X_radsym));

% Surface for LM space
x = linspace(min(points_LM(:,1)),max(points_LM(:,1)),100);
y = linspace(min(points_LM(:,2)),max(points_LM(:,2)),100);
[X_LM,Y_LM] = meshgrid(x,y);
surf = ComputeNakaRushtonJPW(params0,[X_LM(:) Y_LM(:)],'surface7');
surf_LM = reshape(surf,size(X_LM));

% Surface for white space
temp = [X_LM(:) Y_LM(:)] * inv(covMat);
X_white = reshape(temp(:,1),size(X_cov));
Y_white = reshape(temp(:,2),size(Y_cov));
surf_white = surf_LM;
% x = linspace(min(points_white(:,1)),max(points_white(:,1)),100);
% y = linspace(min(points_white(:,2)),max(points_white(:,2)),100);
% [X_white,Y_white] = meshgrid(x,y);
% surf = ComputeNakaRushtonJPW(params0,[X_white(:) Y_white(:)],'surface7');
% surf_white = reshape(surf,size(X_white));


% True angles
% True angle in cov space
ang = params0(end);
rotMat = [cos(ang) sin(ang); -sin(ang) cos(ang)];
projpts = points_cov * rotMat;
rho = max(projpts(:,1));
[x,y] = pol2cart(ang,rho);
true_cov(1) = x;
true_cov(2) = y;

% True angle in radsym space
temp = true_cov * covMat';
true_radsym(1) = temp(1);
true_radsym(2) = temp(2);

% True angle in LM space
true_LM = true_cov;

% True angle in white space
temp = true_LM * M';
true_white(1) = temp(1);
true_white(2) = temp(2);

% Mechanism Transforms
mecht_radsym = STA_cov * covMat';
mecht_cov = STA_radsym * inv(covMat');
mecht_white = STA_LM * M';
mecht_LM = STA_white * inv(M');


% Adjust some values to fit space
% radsym distribution
[~,r] = cart2pol(points_radsym(:,1),points_radsym(:,2));
maxrho_radsym = max(r);
[t,~] = cart2pol(STA_radsym(1),STA_radsym(2));
[STA_radsym(1),STA_radsym(2)] = pol2cart(t,maxrho_radsym);
[t,~] = cart2pol(modelfit_radsym(1),modelfit_radsym(2));
[modelfit_radsym(1),modelfit_radsym(2)] = pol2cart(t,maxrho_radsym);
[t,~] = cart2pol(true_radsym(1),true_radsym(2));
[true_radsym(1),true_radsym(2)] = pol2cart(t,maxrho_radsym);
[t,~] = cart2pol(mecht_radsym(1),mecht_radsym(2));
[mecht_radsym(1),mecht_radsym(2)] = pol2cart(t,maxrho_radsym);

% Cov distribution
[~,r] = cart2pol(points_cov(:,1),points_cov(:,2));
maxrho_cov = max(r);
[t,~] = cart2pol(STA_cov(1),STA_cov(2));
[STA_cov(1),STA_cov(2)] = pol2cart(t,maxrho_cov);
[t,~] = cart2pol(modelfit_cov(1),modelfit_cov(2));
[modelfit_cov(1),modelfit_cov(2)] = pol2cart(t,maxrho_cov);
[t,~] = cart2pol(true_cov(1),true_cov(2));
[true_cov(1),true_cov(2)] = pol2cart(t,maxrho_cov);
[t,~] = cart2pol(mecht_cov(1),mecht_cov(2));
[mecht_cov(1),mecht_cov(2)] = pol2cart(t,maxrho_cov);

% White distributeion
[~,r] = cart2pol(points_white(:,1),points_white(:,2));
maxrho_white = max(r);
[t,~] = cart2pol(STA_white(1),STA_white(2));
[STA_white(1),STA_white(2)] = pol2cart(t,maxrho_white);
[t,~] = cart2pol(modelfit_white(1),modelfit_white(2));
[modelfit_white(1),modelfit_white(2)] = pol2cart(t,maxrho_white);
[t,~] = cart2pol(true_white(1),true_white(2));
[true_white(1),true_white(2)] = pol2cart(t,maxrho_white);
[t,~] = cart2pol(mecht_white(1),mecht_white(2));
[mecht_white(1),mecht_white(2)] = pol2cart(t,maxrho_white);

% LM distribution
[~,r] = cart2pol(points_LM(:,1),points_LM(:,2));
maxrho_LM = max(r);
[t,~] = cart2pol(STA_LM(1),STA_LM(2));
[STA_LM(1),STA_LM(2)] = pol2cart(t,maxrho_LM);
[t,~] = cart2pol(modelfit_LM(1),modelfit_LM(2));
[modelfit_LM(1),modelfit_LM(2)] = pol2cart(t,maxrho_LM);
[t,~] = cart2pol(true_LM(1),true_LM(2));
[true_LM(1),true_LM(2)] = pol2cart(t,maxrho_LM);
[t,~] = cart2pol(mecht_LM(1),mecht_LM(2));
[mecht_LM(1),mecht_LM(2)] = pol2cart(t,maxrho_LM);


% Plot bubble plots with contours
% Set up figure
figure(10); clf;
set(gcf,'Pos',[300 100 700 700]);
radsymax = axes('parent',gcf,'units','normalized','pos',[.05 .525 .425 .425]);
hold on; grid on; box on;
xlabel('Lcc'); ylabel('Mcc'); title('Raially Symmetric Space')
covax = axes('parent',gcf,'units','normalized','pos',[.525 .525 .425 .425]);
hold on; grid on; box on;
xlabel('Lcc'); ylabel('Mcc'); title('Covaried Space')
whiteax = axes('parent',gcf,'units','normalized','pos',[.05 .05 .425 .425]);
hold on; grid on; box on;
xlabel('some arbitry direction'); ylabel('some arbitrary direction'); title('Whitened Space')
LMax = axes('parent',gcf,'units','normalize','pos',[.525 .05 .425 .425]);
hold on; grid on; box on;
xlabel('Lcc'); ylabel('Mcc'); title('Cone Contrast Space'); 

% Some stuff for new figures
scalefac = 20;
[uniqueradsympts,~,bRS] = unique(points_radsym,'rows');
[uniquecovpts,~,bCov] = unique(points_cov,'rows');
[uniquewhitepts,~,bW] = unique(points_white,'rows');
[uniqueLMpts,~,bLM] = unique(points_LM,'rows');

% RadSym Dist Bubble Plots
axes(radsymax); cla; hold on; grid on; box on;
for i = 1:size(uniqueradsympts,1)
    idx = find(i == bRS);
    meannsp1 = mean(nsp0(idx));
    h = plot(uniqueradsympts(i,1),uniqueradsympts(i,2),'ko');
    set(h,'MarkerFaceColor','black','MarkerSize',(((meannsp1-min(nsp0))/max(nsp0))+.2).*scalefac,'MarkerEdgeColor','white');
end
contour(X_radsym,Y_radsym,surf_radsym);
plot([0 STA_radsym(1)],[0 STA_radsym(2)],'b'); % STA
plot([0 mecht_radsym(1)],[0 mecht_radsym(2)],'c--'); % mechanism transform
plot([0 modelfit_radsym(1)],[0 modelfit_radsym(2)],'m--'); % Model fit tuning
plot([0 true_radsym(1)],[0 true_radsym(2)],'r'); % true direction
lim = max(points_radsym(:))*1.1;
xlim([-lim lim]); ylim([-lim lim]);
axis square

% Cov Dist Bubble Plots
axes(covax); cla; hold on; grid on; box on;
for i = 1:size(uniquecovpts,1)
    idx = find(i == bCov);
    meannsp1 = mean(nsp0(idx));
    h = plot(uniquecovpts(i,1),uniquecovpts(i,2),'ko');
    set(h,'MarkerFaceColor','black','MarkerSize',(((meannsp1-min(nsp0))/max(nsp0))+.2).*scalefac,'MarkerEdgeColor','white');
end
contour(X_cov,Y_cov,surf_cov);
plot([0 STA_cov(1)],[0 STA_cov(2)],'c'); % STA
plot([0 mecht_cov(1)],[0 mecht_cov(2)],'b--'); % mech transform
plot([0 modelfit_cov(1)],[0 modelfit_cov(2)],'m--'); % Model fit tuning
plot([0 true_cov(1)],[0 true_cov(2)],'r'); % true direction
lim = max(points_cov(:))*1.1;
xlim([-lim lim]); ylim([-lim lim]);
axis square

% White Dist Bubble Plots
axes(whiteax); cla; hold on; grid on; box on;
for i = 1:size(uniquewhitepts,1)
    idx = find(i == bW);
    meannsp1 = mean(nsp1(idx));
    h = plot(uniquewhitepts(i,1),uniquewhitepts(i,2),'ko');
    set(h,'MarkerFaceColor','black','MarkerSize',(((meannsp1-min(nsp1))/max(nsp1))+.2).*scalefac,'MarkerEdgeColor','white');
end
contour(X_white,Y_white,surf_white);
plot([0 true_white(1)],[0 true_white(2)],'r'); % true direction
plot([0 STA_white(1)],[0 STA_white(2)],'b'); % STA
plot([0 mecht_white(1)],[0 mecht_white(2)],'c--'); % mech transform
plot([0 modelfit_white(1)],[0 modelfit_white(2)],'m--'); % Model fit tuning
lim = max(points_white(:))*1.1;
xlim([-lim lim]); ylim([-lim lim]);
axis square

% LM Dist Bubble Plot
axes(LMax); cla; hold on; grid on; box on;
for i = 1:size(uniqueLMpts,1)
    idx = find(i == bLM);
    meannsp1 = mean(nsp1(idx));
    h = plot(uniqueLMpts(i,1),uniqueLMpts(i,2),'ko');
    set(h,'MarkerFaceColor','black','MarkerSize',(((meannsp1-min(nsp1))/max(nsp1))+.2).*scalefac,'MarkerEdgeColor','white')
end
contour(X_LM,Y_LM,surf_LM);
plot([0 STA_LM(1)],[0 STA_LM(2)],'c'); % STA
plot([0 mecht_LM(1)],[0 mecht_LM(2)],'b--'); % mech transform
plot([0 modelfit_LM(1)],[0 modelfit_LM(2)],'m--'); % Model fit tuning
plot([0 true_LM(1)],[0 true_LM(2)],'r'); % true direction
lim = max(points_LM(:))*1.1;
xlim([-lim lim]); ylim([-lim lim]);
axis square





%% L-cone isolating cell
% Figure 2

% User defined variables
nRnds = 3;
nPres = 1;
ang = 0; % L-cone iso

% Pre-Defined Variables
thetaspace = pi/4;
rhospace = .5;

% Construct Polar Grid
if mod(nRnds,2) == 0
    thetaspace = thetaspace / 2^((nRnds-2)/2);
    rhospace = rhospace / 2^(nRnds/2);
elseif mod(nRnds,2) == 1
    thetaspace = thetaspace / 2^((nRnds-1)/2);
    rhospace = rhospace / 2^((nRnds-1)/2);
end
thetas = shiftdim(0:thetaspace:2*pi-thetaspace);
rhos = shiftdim(rhospace:rhospace:1);

% Ennumerate all conditions
PolRhoIdx = fullfact([numel(thetas) numel(rhos)]);
rhothetalist = repmat([thetas(PolRhoIdx(:,1)) rhos(PolRhoIdx(:,2))],nPres,1);
thetas = rhothetalist(:,1);
rhos = rhothetalist(:,2);

% Transform from polar to cartesian coordinates
[RSLcc,RSMcc] = pol2cart(thetas,rhos);

% Scale RS space into realistic units
newpts = [RSLcc RSMcc] * [.1 0; 0 .1];
RSLcc = newpts(:,1);
RSMcc = newpts(:,2);

% Scale Cone Contrast Units nonlinearly (NonLinearlyScaled)
%scale = lummaxcc*colmaxcc./sqrt((colmaxcc.*cos(thetas-pi/4)).^2 ...
%    +(lummaxcc.*sin(thetas-pi/4)).^2);
%NLSLcc = RSLcc .* scale;
%NLSMcc = RSMcc .* scale;

% Scale cc units linearly (LinearlyScaled)
scalemat = [4 0; 0 1] * [cos(pi/4) sin(pi/4); -sin(pi/4) cos(pi/4)];
newpts = [RSLcc RSMcc] * scalemat;
LSLcc = newpts(:,1);
LSMcc = newpts(:,2);

% RS sigma (set to center of range)
rotMat = [cos(ang) sin(ang); -sin(ang) cos(ang)];
projpts = [RSLcc RSMcc] * rotMat;
RSsig = max(projpts(:,1))/2;

% LS sigma (set to center of range)
rotMat = [cos(ang) sin(ang); -sin(ang) cos(ang)];
projpts = [LSLcc LSMcc] * rotMat;
LSsig = max(projpts(:,1))/2.5;

% Set up some variables
%RSparams = [50 RSsig 100 3 .1 ang];
LSparams = [50 LSsig 100 3 .1 ang];
RSparams = LSparams;

% Responses to stimuli
RSnsp = ComputeNakaRushtonJPW(RSparams,[RSLcc RSMcc],'surface7');
LSnsp = ComputeNakaRushtonJPW(LSparams,[LSLcc LSMcc],'surface7');
minnsp = min(RSnsp);
maxnsp = max(RSnsp);

% Whitening matrix
whitemat = (sqrtm(inv(cov([LSLcc LSMcc]))));
whitepts = [LSLcc LSMcc] * whitemat';
WLcc = whitepts(:,1);
WMcc = whitepts(:,2);

%%%%%% Plotting %%%%%
%%% RS Bubble Plot
scalefac = 10;
blval = 5;
figure(1); clf; hold on; box on;
for i = 1:size(RSLcc,1)
    h(i) = plot(RSLcc(i),RSMcc(i),'ko');
    set(h(i),'MarkerFaceColor','black','MarkerSize',((RSnsp(i)-minnsp).*scalefac./maxnsp)+blval,'MarkerEdgeColor','white')
end
axis square equal tight;

% Display RS surface and pts
x = linspace(-max(RSLcc)*2,max(RSLcc)*2,50);
y = linspace(-max(RSMcc)*2,max(RSMcc)*2,50);
[xx,yy] = meshgrid(x,y);
surface = ComputeNakaRushtonJPW(RSparams,[xx(:) yy(:)],'surface7');
surface = reshape(surface,size(xx));
contour(xx,yy,surface);

% Plot preferred direction in RS space
[pdLc,pdMc] = pol2cart(ang,max(x));
plot([0 pdLc],[0 pdMc],'r-')

% Plot Response Weighted Average 
weightresp = repmat(RSnsp,[1 2]) .* [RSLcc RSMcc];
RSrta = sum(weightresp)./numel(RSLcc);
scale = .125;
plot(RSrta(1)*scale,RSrta(2)*scale,'*r')
xlim([min(RSLcc)*1.05 max(RSLcc)*1.05])
ylim([min(RSMcc)*1.05 max(RSMcc)*1.05])
set(gca,'XTick',linspace(min(RSLcc),max(RSLcc),5))
set(gca,'YTick',linspace(min(RSMcc),max(RSMcc),5))

% Formatting
set(gcf,'PaperPositionMode','auto')
print('-depsc','Fig 2a')
disp('Fig 2a done.')


%%% Linearly Scaled Bubble Plot %%%
figure(2); clf; hold on; box on;
for i = 1:size(LSLcc,1)
    h(i) = plot(LSLcc(i),LSMcc(i),'ko');
    set(h(i),'MarkerFaceColor','black','MarkerSize',((LSnsp(i)-minnsp).*scalefac./maxnsp)+blval,'MarkerEdgeColor','white')
end
axis square equal tight;

% Display LS surface and pts
x = linspace(-max(LSLcc)*2,max(LSLcc)*2,50);
y = linspace(-max(LSMcc)*2,max(LSMcc)*2,50);
[xx,yy] = meshgrid(x,y);
surface = ComputeNakaRushtonJPW(LSparams,[xx(:) yy(:)],'surface7');
LSsurface = reshape(surface,size(xx));
contour(xx,yy,LSsurface);

% Plot preferred direction in LS Space
[pdLc,pdMc] = pol2cart(ang,max(x));
plot([0 pdLc],[0 pdMc],'r-')

% Plot Response Weighted Average 
normLSnsp = LSnsp./max(LSnsp);
weightresp = repmat(normLSnsp,[1 2]) .* [LSLcc LSMcc];
LSrta = mean(weightresp);
scale = .125;
plot(LSrta(1)*scale,LSrta(2)*scale,'*r')
xlim([min(LSLcc)*1.05 max(LSLcc)*1.05])
ylim([min(LSMcc)*1.05 max(LSMcc)*1.05])
set(gca,'XTick',linspace(min(LSLcc),max(LSLcc),7))
set(gca,'XTickLabel',rndofferr(linspace(min(LSLcc),max(LSLcc),7),2))
set(gca,'YTick',linspace(min(LSMcc),max(LSMcc),7))
set(gca,'YTickLabel',rndofferr(linspace(min(LSMcc),max(LSMcc),7),2))

% Formatting
set(gcf,'PaperPositionMode','auto')
print('-depsc','Fig 2b')
disp('Fig 2b done.')

%%% Whitened Bubble Plot %%%
figure(3); clf; hold on; box on;
for i = 1:size(WLcc,1)
    h(i) = plot(WLcc(i),WMcc(i),'ko');
    set(h(i),'MarkerFaceColor','black','MarkerSize',((LSnsp(i)-minnsp).*scalefac./maxnsp)+blval,'MarkerEdgeColor','white')
end
axis square equal tight;

% Display White contours
x = linspace(-max(LSLcc)*2,max(LSLcc)*2,50);
y = linspace(-max(LSMcc)*2,max(LSMcc)*2,50);
[xx,yy] = meshgrid(x,y);
Wpts = [xx(:) yy(:)] * whitemat';
Wxx = reshape(Wpts(:,1),size(xx));
Wyy = reshape(Wpts(:,2),size(yy));
contour(Wxx,Wyy,LSsurface);

% Plot preferred direction in White space
[pdLc,pdMc] = pol2cart(ang,max(x)*10);
whitepd = [pdLc pdMc] * inv(whitemat);
plot([0 whitepd(1)*10],[0 whitepd(2)*10],'r-')

% Plot White Response Weighted Average
Wrwa = [LSrta(1) LSrta(2)] * whitemat;
plot(Wrwa(1)*scale,Wrwa(2)*scale,'r*')
xlim([min(WLcc)*1.05 max(WLcc)*1.05])
ylim([min(WMcc)*1.05 max(WMcc)*1.05])
set(gca,'XTick',linspace(min(WLcc),max(WLcc),7))
set(gca,'XTickLabel',rndofferr(linspace(min(WLcc),max(WLcc),7),2))
set(gca,'YTick',linspace(min(WMcc),max(WMcc),7))
set(gca,'YTickLabel',rndofferr(linspace(min(WMcc),max(WMcc),7),2))

% Formatting
set(gcf,'PaperPositionMode','auto')
print('-depsc','Fig 2c')
disp('Fig 2c done.')

figure(2);
temppd = whitepd * whitemat;
plot([0 temppd(1)*10],[0 temppd(2)*10],'g--')


%% Equivalence of stimulus representations
% Here, I want to build two models: one that represents the stimulus as cone 
% excitations, and one that represents the stimulus as cone contrasts. My
% hypothesis is that, given the right weighting vector, these models can be
% made to make identical predictions.

% Background cone excitations
Lbce = 50;
Mbce = Lbce;

% Stimulus cone excitations
Lsce = linspace(0,2*Lbce,21); % This should sample 100% contrast across all conetypes
Msce = Lsce;
[Lsce0,Msce0] = meshgrid(Lsce,Msce);
Lsce = Lsce0(:);
Msce = Msce0(:);
LMsce = [Lsce Msce];

% Stimulus cone contrasts
Lscc0 = (Lsce0 - Lbce)./Lbce;
Mscc0 = (Msce0 - Mbce)./Mbce;
Lscc = (Lsce - Lbce)./Lbce;
Mscc = (Msce - Mbce)./Mbce;
LMscc = [Lscc Mscc];
Lbcc = 0;
Mbcc = 0;

% Cone contrast weighting vector
Lcew = 1;
Mcew = 2;

% Cone excitation weighting vector
Lccw = (Lcew * Lsce * Lbce) ./ (Lsce - Lbce);
Mccw = (Mcew * Msce * Mbce) ./ (Msce - Mbce);

% Responses in cone excitation
respce = LMsce * [Lcew Mcew]';
surfce = reshape(respce,size(Lsce0));

% Responses in cone contrast
respcc = nan(numel(Lscc),1);
for n = 1:numel(Lscc)
    respcc(n) = LMscc(n,:) * [Lccw(n) Mccw(n)]';
end
surfcc = reshape(respcc,size(Lscc0));

% Plot
figure(1); clf; hold on; grid on;
surfc(Lscc0,Mscc0,surfcc)
alpha(.3)
axis equal square
xlabel('L Contrast'); ylabel('M Contrasts'); title('Cone Contrast Space')

figure(2); clf; hold on; grid on;
surfc(Lsce0,Msce0,surfce)
alpha(.3)
axis equal square
xlabel('L Excitations'); ylabel('M Excitations'); title('Cone Excitation Space')


