function [nextX, totData] = computeNextStim_ALalgorithm_NewModel(x, r, stimDim, support,  numInitData, maxsupport, minsupport)

%% to compute next stimulus using active learning algorithm under a new model for overdispersed data

% --------- input -----------
% 1. x: all the stimuli presented
% 2. r: all the responses measured
% 3. support: grid of points where we estimate tuning curve
% 4. maxsupport: maximum value of the input support 
% 5. minsupport: minimum value of the input support

% --------- output -----------
% 1. nextX: selected stimulus to present next
% 2. totData: datastructure that has important quantities
%    totData.x: all the stimuli presented
%    totData.r: all the responses measured
%    totData.support: grid of points where we estimate tuning curve
%    totData.ndim: stimulus dimension
%    totData.norm_mat_support: matrix with squared distance of support
%    totData.norm_mat: matrix with squared distance of x
%    totData.norm_mat_Kstar: matrix with squared distance between x and support
%    totData.hinit: initial values of hmap
%    totData.g: nonlinearity (exponential)
%    totData.ginv: inverse nonlinearity (log)
%    totData.nstim: number of stimuli presented
%    totData.prs: estimated hyperparameters
%    totData.lambFinal: estimated tuning curve 

%%

global datastruct;

thTrial = length(r);

r(r==0) = 0.1;

%% construct data structure

if thTrial ==numInitData
    datastruct.x = x;
    datastruct.r = r;
    datastruct.support = support;
    datastruct.ndim = stimDim;
    datastruct.norm_mat_support = form_normMat(support, support);  % squared distance
    
    g = @(t)exp(t); 
    ginv = @(t) log(t);
    datastruct.g = g;
    datastruct.ginv = ginv;
    datastruct.hinit = datastruct.ginv(datastruct.r+0.1);
    
    datastruct.norm_mat = form_normMat(datastruct.x, datastruct.x);  % squared distance
    datastruct.norm_mat_Kstar = form_normMat(datastruct.support, datastruct.x);
    
else
    datastruct.x = x;
    datastruct.r = r;
    
    try
        sqrDist_new = form_normMat(datastruct.x(end,:), datastruct.x);
        datastruct.norm_mat = [[datastruct.norm_mat; sqrDist_new(1:end-1)] sqrDist_new'];
        
        normMat_Kstar_new = form_normMat(datastruct.support, datastruct.x(end,:));
        datastruct.norm_mat_Kstar = [datastruct.norm_mat_Kstar normMat_Kstar_new];
        datastruct.hinit = [datastruct.hinit; datastruct.ginv(r(end)+0.1)];
    catch
        disp('wow');
        datastruct.norm_mat = form_normMat(datastruct.x, datastruct.x);  % squared distance
        datastruct.norm_mat_Kstar = form_normMat(datastruct.support, datastruct.x);
        datastruct.hinit = datastruct.ginv(datastruct.r+0.1);
    end     
        
end

% update number of stimulus 
datastruct.nstim = length(datastruct.r);

%% hyperparameter estimation & update fmap

maxNumTrial = stimDim*30; % maximum number of trials for updating hyperparameters

if ((rem(thTrial, 10)==0)&&(thTrial<=maxNumTrial))
    
    % optimize hyperparameters with analytic form
    
    nsevar_init = 0.2;
    ovrscl_1 = datastruct.ginv(mean(datastruct.r)); % overall scale
    lngthscl_1 = (maxsupport-minsupport)/2; % variance
    prs0 = [datastruct.ginv(mean(datastruct.r)); ovrscl_1; lngthscl_1; nsevar_init];
    K = prs0(2)*exp(-.5/prs0(3).*datastruct.norm_mat);
    datastruct.H = K + prs0(end)*eye(size(datastruct.norm_mat));
    datastruct.muf = prs0(1);
    
    [prs, hmapFinal, aFinal, WFinal, sqrtLFinal, neglogev0, detH] = computeFmapAndUpdateTheta_nse(prs0, datastruct);
    datastruct.detH=detH;
    datastruct.prs = prs;
    datastruct.hinit = hmapFinal;
    datastruct.muf = prs(1);
    
else
    
    % update fmap only given hyperparameters 
    
    prs = datastruct.prs;
    K = abs(prs(2))*exp(-.5/abs(prs(3)).*datastruct.norm_mat);
    datastruct.H = K + abs(prs(end))*eye(size(datastruct.norm_mat));
    datastruct.muf = abs(prs(1));
    [neglogev, hmapFinal, aFinal, WFinal, sqrtLFinal]  = computeFmap_nse(datastruct);
    
    datastruct.neglogev = neglogev;
    datastruct.hinit = hmapFinal;
end

%% select next stimulus

datastruct.Kstar = prs(2)*exp(-.5/prs(3).*datastruct.norm_mat_Kstar);
[mean_lambda, nextX, idxNext] = makePrediction_nse(prs, datastruct, aFinal, WFinal, sqrtLFinal);
datastruct.idxNext = idxNext;

datastruct.lambFinal = mean_lambda;
datastruct.aFinal = aFinal;
datastruct.WFinal = WFinal;
datastruct.sqrtLFinal = sqrtLFinal;
totData = datastruct;