% Ideal Observer Analysis for GLMS Datasets
function params = GLMSDGUI_ROC(glmp,sub)
global GLMP
GLMP = glmp;

if ismac
    library = '/Users/jpatrickweller/Dropbox/Patrick/GLMS Data/';
elseif ispc
    library = 'C:\Users\jpweller\Dropbox\Patrick\GLMS Data\';
end

% Organize ROC data
AUC = nan(numel(GLMP.subunit{sub}.uniqueLcc),1);

for s = 1:numel(GLMP.subunit{sub}.uniqueLcc) 
    
    idx = GLMP.subunit{sub}.uniqueIdx{s};
    spHist = GLMP.subunit{sub}.fr(idx)';
    blHist = GLMP.subunit{sub}.blfr';

    thresh = unique([spHist blHist]);
    TPR = nan(1,length(thresh));
    FPR = TPR;

    for n=1:length(thresh)
        TPR(n) = sum(spHist >= thresh(n))/length(spHist);
        FPR(n) = sum(blHist >= thresh(n))/length(blHist);
    end
    
    AUC(s) = trapz(fliplr(FPR),fliplr(TPR));
    
end

x = GLMP.subunit{sub}.uniqueLcc;
y = GLMP.subunit{sub}.uniqueMcc;
z =  AUC;
F = TriScatteredInterp(x,y,z);
[qx,qy] = meshgrid(min(x):.01:max(x),min(y):.01:max(y));


figure(4); clf; hold on;
plot3(x,y,AUC,'k*')
surfc(qx,qy,F(qx,qy))

% Fit data
OneDNeuro(x,y,z)

end


%% 1D fit

% Set up some variables
% vlb = [.0001  .005  .005    2   .5  -pi 0];
% vub = [1       50    50     5   .5   pi 0];
% options = optimset('MaxFunEvals',300,'MaxIter',300,'TolFun',10^-6,'TolX',10^-6,'Display','off','Algorithm','active-set');
% topfr = 1;
% sigmaguess = 0.1;
% expguess = 2;
% baselineguess = .5;
% angguess = linspace(-pi,pi,9);
% 
% resp = AUC;
% LMpts = [GLMP.subunit{sub}.uniqueLcc GLMP.subunit{sub}.uniqueMcc];
% f1 = nan(numel(angguess),numel(vlb));
% fval = nan(numel(angguess),1);
% 
% for n = 1:numel(angguess)
%        
%     paramsGuess = [topfr, sigmaguess, sigmaguess,...
%         expguess,baselineguess,angguess(n) 0];
%     
%     [f1(n,:),fval(n)] = fmincon('FitNakaRushtonFunJPW',paramsGuess,[],[],[],[],vlb,vub,[],options,...
%         LMpts,resp,'surface7','bernoulli');
%     
% end
% 
% [~,best] = max(fval);
% params = f1(best,:);
% fval = fval(best);



function OneDNeuro(xvlas,yvals,zvals)

disp('Fitting a 1D function to Neurometric Data...')


% Load Figure Variables
SurfFig = get(gcf,'UserData');
conpanel = get(SurfFig.conpanel,'UserData');
neuropanel = get(SurfFig.disp.neuropanel,'UserData');
neurostats = get(SurfFig.stats.neuro,'UserData');
    
% Set up some variables
angs = linspace(0,pi,9);
angs(end) = [];
vlb = [1   .001  .001  .001  .5  -pi];
vub = [1     10    10    10  .5   pi];
options = optimset('MaxFunEvals',5000,'MaxIter',5000,'TolFun',10^-6,'TolX',10^-6,'Display','off','Algorithm','active-set');
Aguess = 1;
sigguess = 0.2;
expguess = 2;
blguess = .5;
guessIdx = fullfact([numel(sigguess) numel(sigguess) numel(angs)]);
GOF = nan(size(guessIdx,1),1);
params = nan(numel(guessIdx),numel(vlb));

% Rotate through angles
for rot = 1:size(guessIdx,1)
    
    % Generating an initial guess
    %params0 = [upperA, sigmaguess, sigmaguess, expguess, baselineguess, angs(rot)];
    params0 = [Aguess sigguess(guessIdx(rot,1)) sigguess(guessIdx(rot,2))...
        expguess blguess angs(guessIdx(rot,3))];
    
    % Fit all rotated points
    [f1,fval] = fmincon('FitNakaRushtonFunJPW',params0,[],[],[],[],vlb,vub,[],options,[xvals yvals],zvals,'surface7','bernoulli');
    params(rot,:) = f1;
    GOF(rot) = fval;
    
end

% Find the best fitting surface
[~,bestIdx] = min(GOF);
params1 = params(bestIdx,:);

% Determine confidence intervals
thanks2greg = @(a) FitNakaRushtonFunJPW(a,[xvals yvals],zvals,'surface7','bernoulli');
[hessval,~] = hessian(thanks2greg,params1);
parvar = 1/hessval(end,end);
confint = (2* sqrt(parvar))/pi*180;
neurostats.conf.oneD = confint;

% Display surface and pts
x = linspace(-max(xvals),max(xvals),50);
[xx yy] = meshgrid(x,x);
surface = ComputeNakaRushtonJPW(params1,[xx(:) yy(:)],'surface7');
surface = reshape(surface,size(xx));
axes(neuropanel.axes.oneD); cla; hold on;
neuropanel.pts.oneD = plot3(xvals,yvals,zvals,'k*');
%neuropanel.surf.oneD = surfc(xx,yy,surface);
neuropanel.surf.oneD = contour(xx,yy,surface);
alpha(.4); drawnow

keyboard
end