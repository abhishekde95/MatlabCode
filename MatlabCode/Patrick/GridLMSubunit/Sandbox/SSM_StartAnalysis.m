

function SSM_StartAnalysis(~,~)
global gl

% Load variables
modelfig = get(gcf,'UserData');
conpanel = get(modelfig.conpanel,'UserData');
surfpanel = get(modelfig.surfpanel,'UserData');
fitspanel = get(modelfig.fitspanel,'UserData');

% GL struct stores experimental variables and tracks the progress of the
% model.
gl.nPres = str2num(get(conpanel.uicontrols.npresstim,'string'));
gl.nSamps = str2num(get(conpanel.uicontrols.nsamps,'string'));
gl.allAngs = linspace(-pi,0,17);
GOFs = nan(gl.nSamps,numel(gl.allAngs),3);

for rot = 1:numel(gl.allAngs)
    gl.currentAng = rot;
    
    for rnd = 1:3
        gl.currentRnd = rnd;
        
        for sampn = 1:gl.nSamps
            
            gl.currentSamp = sampn;
            
            disp(['Angle ' num2str(gl.currentAng) ' of ' num2str(numel(gl.allAngs))])
            disp(['# of Rounds = ' num2str(gl.currentRnd)])
            disp(['Sample ' num2str(gl.currentSamp) ' of ' num2str(gl.nSamps)])
            
            disp('Choosing Stim...')
            ChooseLMStimuli;
            
            disp('Modeling Responses...')
            CreateModelSurface;
            
            disp('Fitting Model Data...')
            Fit1Ax
                        
            % Some variables
            fitspanel = get(modelfig.fitspanel,'UserData');
            GOFs(sampn,rot,rnd) = fitspanel.nrnds(rnd).angle(rot).samp(sampn).error;
            meanerror = nanmean(GOFs(:,:,rnd),1);
            
            %plot GOF
            axes(fitspanel.axes.GOF(gl.currentRnd)); cla; hold on; grid on;
            set(fitspanel.axes.GOF(gl.currentRnd),'xlim',[gl.allAngs(1) gl.allAngs(end)]);
            h = plot(gl.allAngs,GOFs(:,:,rnd),'ko');
            h1 = plot(gl.allAngs,meanerror,'*m');
            set(h,'buttondownfcn',@PlotFits)
            set(h1,'buttondownfcn',@PlotFits)
            xlabel('Real Angle')
            ylabel('Radians Difference')
            title(['Fit Using ' num2str(gl.currentRnd) ' Rounds'])
                        
        end
        
    end
    
end


end



function ChooseLMStimuli()
global gl

% Load variables
modelfig = get(gcf,'UserData');
conpanel = get(modelfig.conpanel,'UserData');
surfpanel = get(modelfig.surfpanel,'UserData');
fitspanel = get(modelfig.fitspanel,'UserData');

% Load global variables
nRnds = gl.currentRnd;
nPres = gl.nPres;

% Other variables
thetaspace = pi/4;
rhospace = .5;
maxcc = .09;

% Construct Polar Grid
if nRnds == 2
    
    rhospace = rhospace/2;
    
elseif nRnds == 3
    
    rhospace = rhospace/2;
    thetaspace = thetaspace/2;
    
end

thetas = shiftdim(0:thetaspace:2*pi-thetaspace);
rhos = shiftdim(rhospace:rhospace:1);


% TOTAL HACK TO ADD HIGH CONTRAST LUM STIMULI
%         if gl.rnd == 1
%             rhos = linspace(1,gl.lumcc/gl.colcc,6)';
%             temprho = repmat(rhos(2:end),gl.nPres,1);
%             temptheta1 = ones(size(temprho)) * pi/4;
%             temptheta2 = ones(size(temprho)) * -3*pi/4;
%             highlumstim = cat(2,[temptheta1;temptheta2],[temprho;temprho]);
%             rhothetalist = cat(1,rhothetalist,highlumstim);
%         end
%

% Totally random order
PolRhoIdx = fullfact([numel(thetas) numel(rhos)]);
rhothetalist = repmat([thetas(PolRhoIdx(:,1)) rhos(PolRhoIdx(:,2))],nPres,1);
thetas = rhothetalist(:,1);
rhos = rhothetalist(:,2);

% Transform from polar to cartesian coordinates
[tempLcc,tempMcc] = pol2cart(thetas,rhos);

% Scale Cone Contrast Units for Monitor
surfpanel.nrnds(gl.currentRnd).angle(gl.currentAng).samp(gl.currentSamp).Lcc = tempLcc * maxcc;
surfpanel.nrnds(gl.currentRnd).angle(gl.currentAng).samp(gl.currentSamp).Mcc = tempMcc * maxcc;

% Save variables
set(modelfig.conpanel,'UserData',conpanel);
set(modelfig.surfpanel,'UserData',surfpanel);
set(modelfig.fitspanel,'UserData',fitspanel);
set(gcf,'UserData',modelfig);

%gl.Lcc = tempLcc * maxcc;
%gl.Mcc = tempMcc * maxcc;

end


function CreateModelSurface()
global gl

% Load Figure Variables
modelfig = get(gcf,'UserData');
conpanel = get(modelfig.conpanel,'UserData');
surfpanel = get(modelfig.surfpanel,'UserData');
fitspanel = get(modelfig.fitspanel,'UserData');

% Retrieve stimuli and parameters
ang = gl.allAngs(gl.currentAng);
L = surfpanel.nrnds(gl.currentRnd).angle(gl.currentAng).samp(gl.currentSamp).Lcc;
M = surfpanel.nrnds(gl.currentRnd).angle(gl.currentAng).samp(gl.currentSamp).Mcc;


% Set up some variables
params = [30 0 .045 0 3 0 0];
%x = linspace(-.09,.09,51);
%y = x;
%[xx yy] = meshgrid(x,y);


% Rotate L and M coordinates to make the 'fitting axis' the x-axis
rotMat = [cos(ang) -sin(ang); sin(ang) cos(ang)];
%rotSurfPts = [rotMat * [xx(:) yy(:)]']';
%surfprojs = rotSurfPts(:,1);
rotLMPts = [rotMat * [L(:),M(:)]']';
LMprojs = rotLMPts(:,1);


% Surface Points
%fitPts = ComputeNakaRushtonJPW(params,surfprojs,'asymmetric');
%zz = reshape(fitPts,size(xx));

% Responses to stimuli
respPts = ComputeNakaRushtonJPW(params,LMprojs,'asymmetric');
respPts = poissrnd(respPts);

surfpanel.nrnds(gl.currentRnd).angle(gl.currentAng).samp(gl.currentSamp).Lcc = L;
surfpanel.nrnds(gl.currentRnd).angle(gl.currentAng).samp(gl.currentSamp).Mcc = M;
surfpanel.nrnds(gl.currentRnd).angle(gl.currentAng).samp(gl.currentSamp).responses = respPts;


% Save Figure Variables
set(modelfig.conpanel,'UserData',conpanel);
set(modelfig.surfpanel,'UserData',surfpanel);
set(modelfig.fitspanel,'UserData',fitspanel);
set(gcf,'UserData',modelfig);


end



function Fit1Ax
global gl

% Load Figure Variables
modelfig = get(gcf,'UserData');
conpanel = get(modelfig.conpanel,'UserData');
surfpanel = get(modelfig.surfpanel,'UserData');
fitspanel = get(modelfig.fitspanel,'UserData');

% Set up some variables
angs = linspace(0,pi,181)';
nrots = numel(angs);
GOF = nan(nrots,1);

options = optimset('MaxFunEvals',800,'MaxIter',200,'TolFun',10^-6,'TolX',10^-6,'Display','off','Algorithm','active-set');
params = nan(nrots,7);
Lcc = surfpanel.nrnds(gl.currentRnd).angle(gl.currentAng).samp(gl.currentSamp).Lcc;
Mcc = surfpanel.nrnds(gl.currentRnd).angle(gl.currentAng).samp(gl.currentSamp).Mcc;
responses = surfpanel.nrnds(gl.currentRnd).angle(gl.currentAng).samp(gl.currentSamp).responses;


% Generating an initial guess
maxx = [];
topfr = max(responses);
sigmaguess = 0.045;
baselineguess = 0;
params0 = [topfr, topfr, sigmaguess, sigmaguess, 50, 50, baselineguess];  % need to constrain B across color directions
vlb = [5 5 0.01 0.01 1 1 0];
vub = [topfr*2 topfr*2 1 1 20 20 5];

% Rotate through angles
for rot = 1:nrots
    
    % Rotate L and M coordinates to make the 'fitting axis' the x-axis
    rotMat = [cos(angs(rot)) -sin(angs(rot)); sin(angs(rot)) cos(angs(rot))];
    tempRotPts = [rotMat * [Lcc Mcc]']';
    projs = tempRotPts(:,1);
    

    
    % Use previous parameters
    %if rot>1
        %params1 = mean([params(rot-1,:);params0]);
    %    [~,idx] = max(GOF);
    %    params1 = mean([params(rot-1,:);params(idx,:)]);
    %else
        params1 = params0;
    %end
    
    % Fit all rotated points
    [f1,fval] = fmincon('FitNakaRushtonFunJPW',params1,[],[],[],[],vlb,vub,[],options,projs,responses,'asymmetric');
    params(rot,:) = f1;
    GOF(rot) = -fval;
    
    %Variables for plotting
    maxx = cat(1,maxx,max(projs));
    pts = min(projs):.001:max(projs);
    fitPts = ComputeNakaRushtonJPW(params(rot,:),pts,'asymmetric');
    axlim = max(maxx);
    projline = [rotMat\[-axlim 0; 0 0; axlim 0]']';
    [x,y] = meshgrid(linspace(-axlim,axlim,50));
    rotxy = [rotMat * [x(:) y(:)]']';
    surface = ComputeNakaRushtonJPW(params(rot,:),rotxy(:,1),'asymmetric');
    surface = reshape(surface,size(x));
    
    % Plot figure
%     axes(surfpanel.axes.surffig); cla; hold on; grid on;
%     set(gca,'XLim',[-axlim axlim],'XTickMode','auto','YLim',[-axlim axlim],'YTickMode','auto')
%     surf(x,y,surface)
%     plot3(Lcc,Mcc,responses,'k*')
%     plot3(projline(:,1),projline(:,2),[topfr 0 topfr]','m')
%     xlabel('L Cone Contrast')
%     ylabel('M Cone Contrast')
%     zlabel('Number of Spikes')
%     title('Surface Fit')
  
    % Plot Projections
%     axes(surfpanel.axes.projfig);
%     cla; hold on; grid on;
%     xlim([-max(maxx) max(maxx)])
%     set(gca,'XTickMode','auto','YTickMode','auto')
%     plot(projs,responses,'k*');
%     plot(pts,fitPts,'r--')
%     xlabel('Cone Contrast')
%     ylabel('Number of Spikes')
%     title('Fit to Projections')
    
    % Plot Goodness of Fit
%     axes(fitspanel.axes.GOF(4))
%     if rot == 1
%         cla; hold on; grid on;
%         xlim([min(angs)/pi*180 max(angs)/pi*180])
%         set(gca,'XTick',[0 45 90 135 180],'YTickMode','auto')
%         xlabel('Rotation of Fitting Axis')
%         ylabel('-LL')
%         title('Fits to Sample')
%     end
%     h = plot(angs(1:rot)/pi*180,GOF(1:rot),'ko-');
%     set(h,'buttondownfcn',@PlotSurfaces)
    
end

%Variables for best fitting surface
[~,bestIdx] = max(GOF);
rotMat = [cos(angs(bestIdx)) -sin(angs(bestIdx)); sin(angs(bestIdx)) cos(angs(bestIdx))];
projline = [inv(rotMat) * [-axlim 0; 0 0; axlim 0]']';
[x,y] = meshgrid(linspace(-axlim,axlim,50));
rotxy = [rotMat * [x(:) y(:)]']';
surface = ComputeNakaRushtonJPW(params(bestIdx,:),rotxy(:,1),'asymmetric');
surface = reshape(surface,size(x));
fitPts = ComputeNakaRushtonJPW(params(bestIdx,:),pts,'asymmetric');
tempRotPts = [rotMat * [Lcc Mcc]']';
projs = tempRotPts(:,1);

% % Plot best 3d surface
% axes(surfpanel.axes.surffig); cla; hold on; grid on;
% surf(x,y,surface)
% plot3(Lcc,Mcc,responses,'k*')
% plot3(projline(:,1),projline(:,2),[topfr 0 topfr]','m')
% set(gca,'XLim',[-axlim axlim],'XTickMode','auto','YLim',[-axlim axlim],'YTickMode','auto')
% xlabel('L Cone Contrast')
% ylabel('M Cone Contrast')
% zlabel('Number of Spikes')
% title('Surface Fit')
% 
% % Plot Projections
% axes(surfpanel.axes.projfig)
% cla; hold on; grid on;
% xlim([-max(maxx) max(maxx)])
% plot(projs,responses,'k*');
% plot(pts,fitPts,'r--')
% set(gca,'XTickMode','auto','YTickMode','auto');
% xlabel('Cone Contrast')
% ylabel('Number of Spikes')
% title('Fit to Projections')

% Save results
fitspanel.nrnds(gl.currentRnd).angle(gl.currentAng).samp(gl.currentSamp).testedangs = angs;
fitspanel.nrnds(gl.currentRnd).angle(gl.currentAng).samp(gl.currentSamp).GOF = GOF;
fitspanel.nrnds(gl.currentRnd).angle(gl.currentAng).samp(gl.currentSamp).fitparams = params;
fitspanel.nrnds(gl.currentRnd).angle(gl.currentAng).samp(gl.currentSamp).bestIdx = bestIdx;
fitspanel.nrnds(gl.currentRnd).angle(gl.currentAng).samp(gl.currentSamp).bestAng = angs(bestIdx);
fitspanel.nrnds(gl.currentRnd).angle(gl.currentAng).samp(gl.currentSamp).error = rem(gl.allAngs(gl.currentAng) - angs(bestIdx),pi);

%Save variables
set(modelfig.conpanel,'UserData',conpanel);
set(modelfig.surfpanel,'UserData',surfpanel);
set(modelfig.fitspanel,'UserData',fitspanel);
set(gcf,'UserData',modelfig);

end


function PlotSample(~,~)
global gl

% Load variables
modelfig = get(gcf,'UserData');
conpanel = get(modelfig.conpanel,'UserData');
surfpanel = get(modelfig.surfpanel,'UserData');
fitspanel = get(modelfig.fitspanel,'UserData');

selected = get(fitspanel.axes.GOF(4),'UserData');
nrnds = selected.nrnds;
whichang = selected.ang;
selected.samp = get(conpanel.uicontrols.whichsamp,'value');
set(fitspanel.axes.GOF(4),'UserData',selected)

sampGOF = fitspanel.nrnds(nrnds).angle(whichang).samp(selected.samp).GOF;
sampangs = fitspanel.nrnds(nrnds).angle(whichang).samp(selected.samp).testedangs;
rnderrors = cat(1,fitspanel.nrnds(nrnds).angle(whichang).samp(:).error);
meanrnderror = mean(rnderrors);

% Plot rnd-type GOF and highlight first sample
axes(fitspanel.axes.GOF(gl.currentRnd)); cla; hold on; grid on;
set(fitspanel.axes.GOF(gl.currentRnd),'xlim',[gl.allAngs(1) gl.allAngs(end)],'XTickMode','auto','YTickMode','auto');
h = plot(repmat(gl.allAngs(gl.currentAng),numel(rnderrors),1),rnderrors,'ko');
h1 = plot(gl.allAngs(gl.currentRnd),meanrnderror,'*m');
h2 = plot(gl.allAngs(whichang),rnderrors(selected.samp),'c*');
set(h,'buttondownfcn',@PlotFits)
set(h1,'buttondownfcn',@PlotFits)
set(h2,'buttondownfcn',@PlotFits)
xlabel('Real Angle')
ylabel('Radians Difference')
title(['Fit Using ' num2str(gl.currentRnd) ' Rounds'])

% Plot fits of sample 1 (default)
axes(fitspanel.axes.GOF(4)); cla; hold on; grid on;
set(fitspanel.axes.GOF(4),'xlim',[sampangs(1) sampangs(end)]);
h = plot(sampangs,sampGOF,'ko-');
h1 = plot(sampangs(1),sampGOF(1),'r*');
set(h,'buttondownfcn',@PlotSurfaces);
set(h1,'buttondownfcn',@PlotSurfaces);
xlabel('Rotation of Fitting Axis')
ylabel('-LL')
title('Fits to Sample')

% Variables for surface plotting
Lcc = surfpanel.nrnds(nrnds).angle(whichang).samp(selected.samp).Lcc;
Mcc = surfpanel.nrnds(nrnds).angle(whichang).samp(selected.samp).Mcc;
responses = surfpanel.nrnds(nrnds).angle(whichang).samp(selected.samp).responses;
params = fitspanel.nrnds(nrnds).angle(whichang).samp(selected.samp).fitparams(1,:);
ang = sampangs(1);
rotMat = [cos(ang) -sin(ang); sin(ang) cos(ang)];
tempRotPts = [rotMat * [Lcc Mcc]']';
projs = tempRotPts(:,1);
axlim = max(projs);
projline = [inv(rotMat) * [-axlim 0; 0 0; axlim 0]']';
topfr = max(responses);
[x,y] = meshgrid(linspace(-axlim,axlim,50));
rotxy = [rotMat * [x(:) y(:)]']';
pts = min(projs):.001:max(projs);
surface = ComputeNakaRushtonJPW(params(whichang,:),rotxy(:,1),'asymmetric');
surface = reshape(surface,size(x));
fitPts = ComputeNakaRushtonJPW(params(whichang,:),pts,'asymmetric');

% Plot surface #1 (default)
axes(surfpanel.axes.surffig); cla; hold on; grid on;
set(gca,'XLim',[-axlim axlim],'XTickMode','auto','YLim',[-axlim axlim],'YTickMode','auto')
surf(x,y,surface)
plot3(Lcc,Mcc,responses,'k*')
plot3(projline(:,1),projline(:,2),[topfr 0 topfr]','m')
xlabel('L Cone Contrast')
ylabel('M Cone Contrast')
zlabel('Number of Spikes')
title('Surface Fit')

% Plot Projections of sample 1 at angle 1 (default)
axes(surfpanel.axes.projfig);
cla; hold on; grid on;
xlim([-axlim axlim])
set(gca,'XTickMode','auto','YTickMode','auto')
plot(projs,responses,'k*');
plot(pts,fitPts,'r--')
xlabel('Cone Contrast')
ylabel('Number of Spikes')
title('Fit to Projections')
 
% Allow selection of other samples
% set(conpanel.uicontrols.whichsamp,'string',1:str2double(get(conpanel.uicontrols.nsamps,'string')),'callback',@PlotSample);

% Save variables
set(modelfig.conpanel,'UserData',conpanel);
set(modelfig.surfpanel,'UserData',surfpanel);
set(modelfig.fitspanel,'UserData',fitspanel);
set(gcf,'UserData',modelfig);

end


function PlotFits(~,~)
global gl

h = gca;
whichpt = get(h,'CurrentPoint');
whichpt = whichpt(1,[1 2]);

% Load variables
modelfig = get(gcf,'UserData');
conpanel = get(modelfig.conpanel,'UserData');
surfpanel = get(modelfig.surfpanel,'UserData');
fitspanel = get(modelfig.fitspanel,'UserData');

if h == fitspanel.axes.GOF(1)
    nrnds = 1;
elseif h == fitspanel.axes.GOF(2);
    nrnds = 2;
elseif h == fitspanel.axes.GOF(3);
    nrnds = 3;
else
    keyboard
end    

[~,whichang] = min(diff([repmat(whichpt(1),1,numel(gl.allAngs));gl.allAngs]));
selected = get(fitspanel.axes.GOF(4),'UserData');
selected.nrnds = nrnds;
selected.ang = whichang;
selected.samp = get(conpanel.uicontrols.whichsamp,'value');
set(fitspanel.axes.GOF(4),'UserData',selected)

sampGOF = fitspanel.nrnds(nrnds).angle(whichang).samp(selected.samp).GOF;
sampangs = fitspanel.nrnds(nrnds).angle(whichang).samp(selected.samp).testedangs;
rnderrors = cat(1,fitspanel.nrnds(nrnds).angle(whichang).samp(:).error);
meanrnderror = mean(rnderrors);

% Plot rnd-type GOF and highlight first sample
axes(fitspanel.axes.GOF(gl.currentRnd)); cla; hold on; grid on;
set(fitspanel.axes.GOF(gl.currentRnd),'xlim',[gl.allAngs(1) gl.allAngs(end)],'XTickMode','auto','YTickMode','auto');
h = plot(repmat(gl.allAngs(gl.currentAng),numel(rnderrors),1),rnderrors,'ko');
h1 = plot(gl.allAngs(gl.currentRnd),meanrnderror,'*m');
h2 = plot(gl.allAngs(whichang),rnderrors(1),'c*');
set(h,'buttondownfcn',@PlotFits)
set(h1,'buttondownfcn',@PlotFits)
set(h2,'buttondownfcn',@PlotFits)
xlabel('Real Angle')
ylabel('Radians Difference')
title(['Fit Using ' num2str(gl.currentRnd) ' Rounds'])

% Plot fits of sample 1 (default)
axes(fitspanel.axes.GOF(4)); cla; hold on; grid on;
set(fitspanel.axes.GOF(4),'xlim',[sampangs(1) sampangs(end)]);
h = plot(sampangs,sampGOF,'ko-');
h1 = plot(sampangs(1),sampGOF(1),'r*');
set(h,'buttondownfcn',@PlotSurfaces);
set(h1,'buttondownfcn',@PlotSurfaces);
xlabel('Rotation of Fitting Axis')
ylabel('-LL')
title('Fits to Sample')

% Variables for surface plotting
Lcc = surfpanel.nrnds(nrnds).angle(whichang).samp(selected.samp).Lcc;
Mcc = surfpanel.nrnds(nrnds).angle(whichang).samp(selected.samp).Mcc;
responses = surfpanel.nrnds(nrnds).angle(whichang).samp(selected.samp).responses;
params = fitspanel.nrnds(nrnds).angle(whichang).samp(selected.samp).fitparams(1,:);
ang = sampangs(1);
rotMat = [cos(ang) -sin(ang); sin(ang) cos(ang)];
tempRotPts = [rotMat * [Lcc Mcc]']';
projs = tempRotPts(:,1);
axlim = max(projs);
projline = [inv(rotMat) * [-axlim 0; 0 0; axlim 0]']';
topfr = max(responses);
[x,y] = meshgrid(linspace(-axlim,axlim,50));
rotxy = [rotMat * [x(:) y(:)]']';
pts = min(projs):.001:max(projs);
surface = ComputeNakaRushtonJPW(params(whichang,:),rotxy(:,1),'asymmetric');
surface = reshape(surface,size(x));
fitPts = ComputeNakaRushtonJPW(params(whichang,:),pts,'asymmetric');

% Plot surface #1 (default)
axes(surfpanel.axes.surffig); cla; hold on; grid on;
set(gca,'XLim',[-axlim axlim],'XTickMode','auto','YLim',[-axlim axlim],'YTickMode','auto')
surf(x,y,surface)
plot3(Lcc,Mcc,responses,'k*')
plot3(projline(:,1),projline(:,2),[topfr 0 topfr]','m')
xlabel('L Cone Contrast')
ylabel('M Cone Contrast')
zlabel('Number of Spikes')
title('Surface Fit')

% Plot Projections of sample 1 at angle 1 (default)
axes(surfpanel.axes.projfig);
cla; hold on; grid on;
xlim([-axlim axlim])
set(gca,'XTickMode','auto','YTickMode','auto')
plot(projs,responses,'k*');
plot(pts,fitPts,'r--')
xlabel('Cone Contrast')
ylabel('Number of Spikes')
title('Fit to Projections')

% Allow selection of other samples
set(conpanel.uicontrols.whichsamp,'string',1:str2double(get(conpanel.uicontrols.nsamps,'string')),'callback',@PlotSample);

% Save variables
set(modelfig.conpanel,'UserData',conpanel);
set(modelfig.surfpanel,'UserData',surfpanel);
set(modelfig.fitspanel,'UserData',fitspanel);
set(gcf,'UserData',modelfig);

end

function PlotSurfaces(~,~)

h = gca;
whichpt = get(h,'CurrentPoint');
whichpt = whichpt(1,[1 2]);

% Load variables
modelfig = get(gcf,'UserData');
conpanel = get(modelfig.conpanel,'UserData');
surfpanel = get(modelfig.surfpanel,'UserData');
fitspanel = get(modelfig.fitspanel,'UserData');

% Load variables for selected surface
selected = get(fitspanel.axes.GOF(4),'UserData');
nrnds = selected.nrnds;
ang = selected.ang;
samp = selected.samp;

% Variables for surface plotting
sampangs = fitspanel.nrnds(nrnds).angle(ang).samp(samp).testedangs;
[~,selectedang] = min(abs(diff([repmat(whichpt(1),1,numel(sampangs));sampangs'])));
sampGOF = fitspanel.nrnds(nrnds).angle(ang).samp(samp).GOF;
sampangs = fitspanel.nrnds(nrnds).angle(ang).samp(samp).testedangs;
Lcc = surfpanel.nrnds(nrnds).angle(ang).samp(samp).Lcc;
Mcc = surfpanel.nrnds(nrnds).angle(ang).samp(samp).Mcc;
responses = surfpanel.nrnds(nrnds).angle(ang).samp(samp).responses;
params = fitspanel.nrnds(nrnds).angle(ang).samp(samp).fitparams(selectedang,:);
rotMat = [cos(sampangs(selectedang)) -sin(sampangs(selectedang)); sin(sampangs(selectedang)) cos(sampangs(selectedang))];
tempRotPts = [rotMat * [Lcc Mcc]']';
projs = tempRotPts(:,1);
axlim = max(projs);
projline = [inv(rotMat) * [-axlim 0; 0 0; axlim 0]']';
topfr = max(responses);
[x,y] = meshgrid(linspace(-axlim,axlim,50));
rotxy = [rotMat * [x(:) y(:)]']';
pts = min(projs):.001:max(projs);
surface = ComputeNakaRushtonJPW(params,rotxy(:,1),'asymmetric');
surface = reshape(surface,size(x));
fitPts = ComputeNakaRushtonJPW(params,pts,'asymmetric');

% Plot surface
axes(surfpanel.axes.surffig); cla; hold on; grid on;
set(gca,'XLim',[-axlim axlim],'XTickMode','auto','YLim',[-axlim axlim],'YTickMode','auto')
surf(x,y,surface)
plot3(Lcc,Mcc,responses,'k*')
plot3(projline(:,1),projline(:,2),[topfr 0 topfr]','m')
xlabel('L Cone Contrast')
ylabel('M Cone Contrast')
zlabel('Number of Spikes')

% Plot Projections
axes(surfpanel.axes.projfig);
cla; hold on; grid on;
xlim([-axlim axlim])
set(gca,'XTickMode','auto','YTickMode','auto')
plot(projs,responses,'k*');
plot(pts,fitPts,'r--')
xlabel('Cone Contrast')
ylabel('Number of Spikes')

% Plot fits of sample 1
axes(fitspanel.axes.GOF(4)); cla; hold on; grid on;
set(fitspanel.axes.GOF(4),'xlim',[sampangs(1) sampangs(end)]);
h = plot(sampangs,sampGOF,'ko-');
h1 = plot(sampangs(selectedang),sampGOF(selectedang),'r*');
set(h,'buttondownfcn',@PlotSurfaces);
set(h1,'buttondownfcn',@PlotSurfaces);

%Save variables
set(modelfig.conpanel,'UserData',conpanel);
set(modelfig.surfpanel,'UserData',surfpanel);
set(modelfig.fitspanel,'UserData',fitspanel);
set(gcf,'UserData',modelfig);

end

