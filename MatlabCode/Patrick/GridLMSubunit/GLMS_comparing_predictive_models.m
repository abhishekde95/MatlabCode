% For making GLMS figures
global GLMSPopData

% Grab saved population data
if ismac
    library = '/Users/jpatrickweller/Dropbox/Patrick/GLMS Data/';
elseif ispc
    library = 'C:\Users\jpweller\Dropbox\Patrick\GLMS Data\';
end
load([library 'GLMSPopData.mat'])
datatypes = GLMSPopData(1,:);

GLMSPopGUI_Params;
ParamsFig = get(552,'userdata');
poppanel = get(ParamsFig.poppanel,'userdata');
includeidx = find(poppanel.oneDL | poppanel.twoDL);

glmps = GLMSPopData(includeidx+1,strcmp(datatypes,'GLMP'));
subs = GLMSPopData(includeidx+1,strcmp(datatypes,'Subunit'));

%% 2D model sampling col and lum dirs

surftype = 'conicsection_xy';
errortype = 'negativebinomial';
options = optimset('Algorithm','interior-point','MaxFunEvals',5000,...
    'MaxIter',5000,'FinDiffType','central','Hessian','bfgs','display','off',...
    'TolFun',0,'TolCon',0,'FunValCheck','on','AlwaysHonorConstraints','bounds');
PAs = [-3*pi/4 -pi/4 pi/4 3*pi/4];
a(1,:) = [0 -1 1 0 0 0 0 0];
a(2,:) = [-1 0 0 0 0 1 0 0];
b = zeros(size(a,1),1);
cardirs.params = nan(numel(includeidx),8);
cardirs.nll = nan(numel(includeidx),1);
cardirs.sse = nan(numel(includeidx),1);

for n = 1:numel(glmps)
    
    % Pull out datafile
    GLMP = glmps{n};
    sub = subs{n};
    n
    
    % Only look at col and lum axes
    L = any(GLMP.subunit{sub}.theta == PAs,2);
    lcc = cat(1,GLMP.subunit{sub}.Lcc(L),zeros(sum(L),1));
    mcc = cat(1,GLMP.subunit{sub}.Mcc(L),zeros(sum(L),1));
    nsp = cat(1,GLMP.subunit{sub}.nspikes(L),GLMP.subunit{sub}.blnspikes(L));
    
    % set up variables
    ub = [max(nsp)             300  300  300 10 max(nsp)  pi 5];
    lb = [min(nsp) 1/max(GLMP.rho)    0    0  1     .001 -pi 0];
    Aguess = max(GLMP.subunit{sub}.nspikes)*.8;
    sigguess = [0 1/.05];
    expguess = 3;
    blguess = mean(GLMP.subunit{sub}.blnspikes);
    angguess = linspace(-pi,pi,5);
    kappaguess = regress(GLMP.subunit{sub}.varnspikes-GLMP.subunit{sub}.meannspikes,...
        GLMP.subunit{sub}.meannspikes.^2);% Variance for negative binomial fit
    if kappaguess < lb(end)
        kappaguess = lb(end);
    end
    guessIdx = fullfact([numel(sigguess) numel(sigguess) numel(angguess)]);

    % preassign variables
    tparams = nan(size(guessIdx,1),8);
    tnll = nan(size(guessIdx,1),1);
        
    % Test each direction as the preferred direction
    for i = 1:size(guessIdx,1)

        
        paramsGuess = [Aguess sigguess(2) sigguess(guessIdx(i,1)) sigguess(guessIdx(i,2))  expguess blguess angguess(guessIdx(i,3)) kappaguess];
        [tparams(i,:),tnll(i)] = fmincon('FitNakaRushtonFunJPW',paramsGuess,...
            a,b,[],[],lb,ub,[],options,[lcc mcc],nsp,surftype,errortype);
                
    end
    
    % Use least nll
    [minnll,which] = min(tnll);
    params = tparams(which,:);
    
    % Calculate sse
    [predresp] = ComputeNakaRushtonJPW(params,[GLMP.subunit{sub}.Lcc GLMP.subunit{sub}.Mcc],surftype);
    
    % save data
    cardirs.params(n,:) = params;
    cardirs.nll(n) = minnll;
    cardirs.sse(n) = sum((predresp-GLMP.subunit{sub}.nspikes).^2);
    %save([library '/predmod/cardirs'],'cardirs')

    % Plot resulting fit
    figure(1); clf; hold on;
    set(gcf,'units','normalized','pos',[.2 .2 .6 .6],'numbertitle','off',...
        'renderer','paint')
    
    % Bubbble plot surface (this is annoying, but bubbles, contours, and
    % surface must all be saved sepeartely for Illustrator to recognize it.
    meannsp = GLMP.subunit{sub}.meannspikes;
    ticks = -.8:.2:.8; % regular axes
    cmap = repmat(linspace(.7,0,32)',1,3);
    colormap(cmap)
    
    % Define thetas and rhos for the whole grid
    thetas = linspace(-pi,pi,361)';
    majorax = max(GLMP.subunit{sub}.uniquerho(GLMP.subunit{sub}.uniquetheta==pi/4))*1.1;
    minorax = max(GLMP.subunit{sub}.uniquerho(GLMP.subunit{sub}.uniquetheta==-pi/4))*1.1;
    nom = majorax * minorax;
    denom = (majorax*sin(thetas-pi/4)).^2 + (minorax*cos(thetas-pi/4)).^2;
    rhos = nom ./ sqrt(denom);
    scalars = 0:.01:1;
    rhosgrid = rhos * scalars;
    thetasgrid = repmat(thetas,1,numel(scalars));
    [x,y] = pol2cart(thetasgrid,rhosgrid);
    surface = ComputeNakaRushtonJPW(params,[x(:) y(:)],surftype);
    surface = reshape(surface,size(rhosgrid));
    axlim = max(x(:));
    
    %%% (1) Surface %%%
    h = surf(x,y,surface); hold on;
    set(h,'edgecolor','none');
    alpha(.3);
    
    %%% (2) contours and pd %%%
    rotmat = [sin(params(end-1)) cos(params(end-1)); -cos(params(end-1)) sin(params(end-1))];
    PD = rotmat * [0 .5]';
    contour(x,y,surface,5,'linewidth',2); hold on; % masked grid
    plot([0 PD(1)],[0 PD(2)],'k-')
    set(gca,'tickdir','out','XTick',ticks,'YTick',ticks,'xlim',[-axlim axlim],'ylim',[-axlim axlim]);
    
    %%% (3) Bubbles %%%
    maxnsp = max(max(GLMP.subunit{sub}.meannspikes));
    polar(0,axlim,'k*'); hold on; % this is stupid, just to get axis limits right
    scalar = 30;
    for i = 1:numel(GLMP.subunit{sub}.uniqueLcc)
        mn = meannsp(i)/maxnsp*scalar+scalar/10;
        h = polar(GLMP.subunit{sub}.uniquetheta(i),GLMP.subunit{sub}.uniquerho(i),'ko');
        if any(GLMP.subunit{sub}.uniquetheta(i) == PAs)
            set(h,'MarkerFaceColor','b','MarkerSize',mn)
        else
            set(h,'MarkerFaceColor','none','MarkerSize',mn)
        end
    end

    box on; grid off;
    axis equal square tight
    set(gca,'xlim',[-axlim axlim],'ylim',[-axlim axlim])

end
    

%% 2D model fixed to cardinal dirs

surftype = 'conicsection_xy';
errortype = 'negativebinomial';
options = optimset('Algorithm','interior-point','MaxFunEvals',5000,...
    'MaxIter',5000,'FinDiffType','central','Hessian','bfgs','display','off',...
    'TolFun',0,'TolCon',0,'FunValCheck','on','AlwaysHonorConstraints','bounds');
PAs = [-3*pi/4 -pi/4 pi/4 3*pi/4];
a(1,:) = [0 -1 1 0 0 0 0 0];
a(2,:) = [-1 0 0 0 0 1 0 0];
b = zeros(size(a,1),1);
cardfix.params = nan(numel(includeidx),8);
cardfix.nll = nan(numel(includeidx),1);
cardfix.sse = nan(numel(includeidx),1);

for n = 1:numel(glmps)
    
    n
    % Pull out datafile
    GLMP = glmps{n};
    sub = subs{n};
    
    % Only look at col and lum axes
    L = any(GLMP.subunit{sub}.theta == PAs,2);
    lcc = cat(1,GLMP.subunit{sub}.Lcc(L),zeros(sum(L),1));
    mcc = cat(1,GLMP.subunit{sub}.Mcc(L),zeros(sum(L),1));
    nsp = cat(1,GLMP.subunit{sub}.nspikes(L),GLMP.subunit{sub}.blnspikes(L));
    
    % set up variables
    ub = [max(nsp)             300  300  300 10 max(nsp)  pi 5];
    lb = [min(nsp) 1/max(GLMP.rho)    0    0  1     .001 -pi 0];
    Aguess = max(GLMP.subunit{sub}.nspikes)*.8;
    sigguess = [0 1/.05];
    expguess = 3;
    blguess = mean(GLMP.subunit{sub}.blnspikes);
    angguess = PAs;
    kappaguess = regress(GLMP.subunit{sub}.varnspikes-GLMP.subunit{sub}.meannspikes,...
        GLMP.subunit{sub}.meannspikes.^2);% Variance for negative binomial fit
    if kappaguess < lb(end)
        kappaguess = lb(end);
    end
    guessIdx = fullfact([numel(sigguess) numel(sigguess) numel(angguess)]);
    
    % preassign variables
    tparams = nan(size(guessIdx,1),8);
    tnll = nan(size(guessIdx,1),1);
        
    % Test each direction as the preferred direction
    for i = 1:size(guessIdx,1)

        ub = [max(nsp)             300  300  300 10 max(nsp) angguess(guessIdx(i,3)) 5];
        lb = [min(nsp) 1/max(GLMP.rho)    0    0  1     .001 angguess(guessIdx(i,3)) 0];
        
        paramsGuess = [Aguess sigguess(2) sigguess(guessIdx(i,1)) sigguess(guessIdx(i,2))  expguess blguess angguess(guessIdx(i,3)) kappaguess];
        [parvals,nll] = fmincon('FitNakaRushtonFunJPW',paramsGuess,...
            a,b,[],[],lb,ub,[],options,[lcc mcc],nsp,surftype,errortype);
        if ~isempty(nll)
            tparams(i,:) = parvals;
            tnll(i) = nll;
        end
                
    end
    
    % Use least nll
    [minnll,which] = min(tnll);
    params = tparams(which,:);
    
    % Calculate sse
    [predresp] = ComputeNakaRushtonJPW(params,[GLMP.subunit{sub}.Lcc GLMP.subunit{sub}.Mcc],surftype);
    
    % save data
    cardfix.params(n,:) = params;
    cardfix.nll(n) = minnll;
    cardfix.sse(n) = sum((predresp-GLMP.subunit{sub}.nspikes).^2);
    %save([library '/predmod/cardfix'],'cardfix')

    % Plot resulting fit
    figure(1); clf; hold on;
    set(gcf,'units','normalized','pos',[.2 .2 .6 .6],'numbertitle','off',...
        'renderer','paint')
    
    % Bubbble plot surface 
    meannsp = GLMP.subunit{sub}.meannspikes;
    ticks = -.8:.2:.8; % regular axes
    cmap = repmat(linspace(.7,0,32)',1,3);
    colormap(cmap)
    
    % Define thetas and rhos for the whole grid
    thetas = linspace(-pi,pi,361)';
    majorax = max(GLMP.subunit{sub}.uniquerho(GLMP.subunit{sub}.uniquetheta==pi/4))*1.1;
    minorax = max(GLMP.subunit{sub}.uniquerho(GLMP.subunit{sub}.uniquetheta==-pi/4))*1.1;
    nom = majorax * minorax;
    denom = (majorax*sin(thetas-pi/4)).^2 + (minorax*cos(thetas-pi/4)).^2;
    rhos = nom ./ sqrt(denom);
    scalars = 0:.01:1;
    rhosgrid = rhos * scalars;
    thetasgrid = repmat(thetas,1,numel(scalars));
    [x,y] = pol2cart(thetasgrid,rhosgrid);
    surface = ComputeNakaRushtonJPW(params,[x(:) y(:)],surftype);
    surface = reshape(surface,size(rhosgrid));
    axlim = max(x(:));
    
    %%% (1) Surface %%%
    h = surf(x,y,surface); hold on;
    set(h,'edgecolor','none');
    alpha(.3);
    
    %%% (2) contours and pd %%%
    rotmat = [sin(params(end-1)) cos(params(end-1)); -cos(params(end-1)) sin(params(end-1))];
    PD = rotmat * [0 .5]';
    contour(x,y,surface,5,'linewidth',2); hold on; % masked grid
    plot([0 PD(1)],[0 PD(2)],'k-')
    set(gca,'tickdir','out','XTick',ticks,'YTick',ticks,'xlim',[-axlim axlim],'ylim',[-axlim axlim]);
    
    %%% (3) Bubbles %%%
    maxnsp = max(max(GLMP.subunit{sub}.meannspikes));
    polar(0,axlim,'k*'); hold on; % this is stupid, just to get axis limits right
    scalar = 30;
    for i = 1:numel(GLMP.subunit{sub}.uniqueLcc)
        mn = meannsp(i)/maxnsp*scalar+scalar/10;
        h = polar(GLMP.subunit{sub}.uniquetheta(i),GLMP.subunit{sub}.uniquerho(i),'ko');
        if any(GLMP.subunit{sub}.uniquetheta(i) == PAs)
            set(h,'MarkerFaceColor','b','MarkerSize',mn)
        else
            set(h,'MarkerFaceColor','none','MarkerSize',mn)
        end
    end

    box on; grid off;
    axis equal square tight
    set(gca,'xlim',[-axlim axlim],'ylim',[-axlim axlim])

end


%% 2D model sampling cone iso directions

surftype = 'conicsection_xy';
errortype = 'negativebinomial';
options = optimset('Algorithm','interior-point','MaxFunEvals',5000,...
    'MaxIter',5000,'FinDiffType','central','Hessian','bfgs','display','off',...
    'TolFun',0,'TolCon',0,'FunValCheck','on','AlwaysHonorConstraints','bounds');
PAs = [-pi/2 0 pi/2 pi];
a(1,:) = [0 -1 1 0 0 0 0 0];
a(2,:) = [-1 0 0 0 0 1 0 0];
b = zeros(size(a,1),1);
coneiso.params = nan(numel(includeidx),8);
coneiso.nll = nan(numel(includeidx),1);
coneiso.sse = nan(numel(includeidx),1);

for n = 87%1:numel(glmps)
    
    n
    % Pull out datafile
    GLMP = glmps{n};
    sub = subs{n};
    
    % Only look at col and lum axes
    L = any(GLMP.subunit{sub}.theta == PAs,2);
    lcc = cat(1,GLMP.subunit{sub}.Lcc(L),zeros(sum(L),1));
    mcc = cat(1,GLMP.subunit{sub}.Mcc(L),zeros(sum(L),1));
    nsp = cat(1,GLMP.subunit{sub}.nspikes(L),GLMP.subunit{sub}.blnspikes(L));
    
    % set up variables
    ub = [(max(nsp)+1)             300  300  300 10 max(nsp)+1  pi 5];
    lb = [min(nsp)     1/max(GLMP.rho)    0    0  1       .001 -pi 0];
    Aguess = (max(GLMP.subunit{sub}.nspikes(L))+.1)*.8;
    sigguess = [0 1/.05];
    expguess = 3;
    blguess = mean(GLMP.subunit{sub}.blnspikes)+.1;
    angguess = linspace(-pi,pi,5);
    kappaguess = regress(GLMP.subunit{sub}.varnspikes-GLMP.subunit{sub}.meannspikes,...
        GLMP.subunit{sub}.meannspikes.^2);% Variance for negative binomial fit
    if kappaguess < lb(end)
        kappaguess = lb(end);
    end
    guessIdx = fullfact([numel(sigguess) numel(sigguess) numel(angguess)]);

    % preassign variables
    tparams = nan(size(guessIdx,1),8);
    tnll = nan(size(guessIdx,1),1);
        
    % Test each direction as the preferred direction
    for i = 1:size(guessIdx,1)

        paramsGuess = [Aguess sigguess(2) sigguess(guessIdx(i,1)) sigguess(guessIdx(i,2))  expguess blguess angguess(guessIdx(i,3)) kappaguess];
        [parvals,nll] = fmincon('FitNakaRushtonFunJPW',paramsGuess,...
            a,b,[],[],lb,ub,[],options,[lcc mcc],nsp,surftype,errortype);
        if ~isempty(nll)
            tparams(i,:) = parvals;
            tnll(i) = nll;
        end
            
                
    end
    
    % Use least nll
    [minnll,which] = min(tnll);
    params = tparams(which,:);
    
    % Calculate sse
    [predresp] = ComputeNakaRushtonJPW(params,[GLMP.subunit{sub}.Lcc GLMP.subunit{sub}.Mcc],surftype);
    
    % save data
    coneiso.params(n,:) = params;
    coneiso.nll(n) = minnll;
    coneiso.sse(n) = sum((predresp-GLMP.subunit{sub}.nspikes).^2);
   

%     % Plot resulting fit
%     figure(1); clf; hold on;
%     set(gcf,'units','normalized','pos',[.2 .2 .6 .6],'numbertitle','off',...
%         'renderer','paint')
%     
%     % Bubbble plot surface (this is annoying, but bubbles, contours, and
%     % surface must all be saved sepeartely for Illustrator to recognize it.
%     meannsp = GLMP.subunit{sub}.meannspikes;
%     ticks = -.8:.2:.8; % regular axes
%     cmap = repmat(linspace(.7,0,32)',1,3);
%     colormap(cmap)
%     
%     % Define thetas and rhos for the whole grid
%     thetas = linspace(-pi,pi,361)';
%     majorax = max(GLMP.subunit{sub}.uniquerho(GLMP.subunit{sub}.uniquetheta==pi/4))*1.1;
%     minorax = max(GLMP.subunit{sub}.uniquerho(GLMP.subunit{sub}.uniquetheta==-pi/4))*1.1;
%     nom = majorax * minorax;
%     denom = (majorax*sin(thetas-pi/4)).^2 + (minorax*cos(thetas-pi/4)).^2;
%     rhos = nom ./ sqrt(denom);
%     scalars = 0:.01:1;
%     rhosgrid = rhos * scalars;
%     thetasgrid = repmat(thetas,1,numel(scalars));
%     [x,y] = pol2cart(thetasgrid,rhosgrid);
%     surface = ComputeNakaRushtonJPW(params,[x(:) y(:)],surftype);
%     surface = reshape(surface,size(rhosgrid));
%     axlim = max(x(:));
%     
%     %%% (1) Surface %%%
%     h = surf(x,y,surface); hold on;
%     set(h,'edgecolor','none');
%     alpha(.3);
%     
%     %%% (2) contours and pd %%%
%     rotmat = [sin(params(end-1)) cos(params(end-1)); -cos(params(end-1)) sin(params(end-1))];
%     PD = rotmat * [0 .5]';
%     contour(x,y,surface,5,'linewidth',2); hold on; % masked grid
%     plot([0 PD(1)],[0 PD(2)],'k-')
%     set(gca,'tickdir','out','XTick',ticks,'YTick',ticks,'xlim',[-axlim axlim],'ylim',[-axlim axlim]);
%     
%     %%% (3) Bubbles %%%
%     maxnsp = max(max(GLMP.subunit{sub}.meannspikes));
%     polar(0,axlim,'k*'); hold on; % this is stupid, just to get axis limits right
%     scalar = 30;
%     for i = 1:numel(GLMP.subunit{sub}.uniqueLcc)
%         mn = meannsp(i)/maxnsp*scalar+scalar/10;
%         h = polar(GLMP.subunit{sub}.uniquetheta(i),GLMP.subunit{sub}.uniquerho(i),'ko');
%         if any(GLMP.subunit{sub}.uniquetheta(i) == PAs)
%             set(h,'MarkerFaceColor','b','MarkerSize',mn)
%         else
%             set(h,'MarkerFaceColor','none','MarkerSize',mn)
%         end
%     end
% 
%     box on; grid off;
%     axis equal square tight
%     set(gca,'xlim',[-axlim axlim],'ylim',[-axlim axlim])

end

 save([library '/predmod/coneiso'],'coneiso')

%% 1D model sampling sparse coneiso stim

surftype = 'conicsection_xy';
errortype = 'negativebinomial';
options = optimset('Algorithm','interior-point','MaxFunEvals',5000,...
    'MaxIter',5000,'FinDiffType','central','Hessian','bfgs','display','off',...
    'TolFun',0,'TolCon',0,'FunValCheck','on','AlwaysHonorConstraints','bounds');
PAs = [-pi/2 0 pi/2 pi];
a(1,:) = [0 -1 1 0 0 0 0 0];
a(2,:) = [-1 0 0 0 0 1 0 0];
b = zeros(size(a,1),1);
coneweights.params = nan(numel(includeidx),8);
coneweights.nll = nan(numel(includeidx),1);
coneweights.sse = nan(numel(includeidx),1);

for n = 1:numel(glmps)
    
    n
    % Pull out datafile
    GLMP = glmps{n};
    sub = subs{n};
    
    % Only look at col and lum axes
    maxcc = max(GLMP.subunit{sub}.rho(GLMP.subunit{sub}.theta == PAs(1)));
    L = any(GLMP.subunit{sub}.theta == PAs,2) & GLMP.subunit{sub}.rho == maxcc;
    uniqL = any(GLMP.subunit{sub}.uniquetheta == PAs,2) & GLMP.subunit{sub}.uniquerho == maxcc;
    lcc = cat(1,GLMP.subunit{sub}.Lcc(L),zeros(sum(L),1));
    mcc = cat(1,GLMP.subunit{sub}.Mcc(L),zeros(sum(L),1));
    nsp = cat(1,GLMP.subunit{sub}.nspikes(L),GLMP.subunit{sub}.blnspikes(L));
    
    % set up variables
    ub = [max(nsp)             300  300  0 10 max(nsp)  pi 5];
    lb = [min(nsp) 1/max(GLMP.rho)    0  0  1     .001 -pi 0];
    Aguess = max(GLMP.subunit{sub}.nspikes)*.8;
    sigguess = [0 1/.05];
    expguess = 3;
    blguess = mean(GLMP.subunit{sub}.blnspikes);
    angguess = linspace(-pi,pi,5);
    kappaguess = regress(GLMP.subunit{sub}.varnspikes-GLMP.subunit{sub}.meannspikes,...
        GLMP.subunit{sub}.meannspikes.^2);% Variance for negative binomial fit
    if kappaguess < lb(end)
        kappaguess = lb(end);
    end
    guessIdx = fullfact([numel(sigguess) numel(angguess)]);

    % preassign variables
    tparams = nan(size(guessIdx,1),8);
    tnll = nan(size(guessIdx,1),1);
        
    % Test each direction as the preferred direction
    for i = 1:size(guessIdx,1)

        paramsGuess = [Aguess sigguess(2) sigguess(guessIdx(i,1)) 0 expguess blguess angguess(guessIdx(i,2)) kappaguess];
        [parvals,nll] = fmincon('FitNakaRushtonFunJPW',paramsGuess,...
            a,b,[],[],lb,ub,[],options,[lcc mcc],nsp,surftype,errortype);
        if ~isempty(nll)
            tparams(i,:) = parvals;
            tnll(i) = nll;
        end
                
    end
    
    % Use least nll
    [minnll,which] = min(tnll);
    params = tparams(which,:);
    
    % Calculate sse
    [predresp] = ComputeNakaRushtonJPW(params,[GLMP.subunit{sub}.Lcc GLMP.subunit{sub}.Mcc],surftype);
    
    % save data
    coneweights.params(n,:) = params;
    coneweights.nll(n) = minnll;
    coneweights.sse(n) = sum((predresp-GLMP.subunit{sub}.nspikes).^2);
    %save([library '/predmod/coneweights'],'coneweights')

    % Plot resulting fit
    figure(1); clf; hold on;
    set(gcf,'units','normalized','pos',[.2 .2 .6 .6],'numbertitle','off',...
        'renderer','paint')
    
    %Bubbble plot surface
    meannsp = GLMP.subunit{sub}.meannspikes;
    ticks = -.8:.2:.8; % regular axes
    cmap = repmat(linspace(.7,0,32)',1,3);
    colormap(cmap)
    
    % Define thetas and rhos for the whole grid
    thetas = linspace(-pi,pi,361)';
    majorax = max(GLMP.subunit{sub}.uniquerho(GLMP.subunit{sub}.uniquetheta==pi/4))*1.1;
    minorax = max(GLMP.subunit{sub}.uniquerho(GLMP.subunit{sub}.uniquetheta==-pi/4))*1.1;
    nom = majorax * minorax;
    denom = (majorax*sin(thetas-pi/4)).^2 + (minorax*cos(thetas-pi/4)).^2;
    rhos = nom ./ sqrt(denom);
    scalars = 0:.01:1;
    rhosgrid = rhos * scalars;
    thetasgrid = repmat(thetas,1,numel(scalars));
    [x,y] = pol2cart(thetasgrid,rhosgrid);
    surface = ComputeNakaRushtonJPW(params,[x(:) y(:)],surftype);
    surface = reshape(surface,size(rhosgrid));
    axlim = max(x(:));
    
    %%% (1) Surface %%%
    h = surf(x,y,surface); hold on;
    set(h,'edgecolor','none');
    alpha(.3);
    
    %%% (2) contours and pd %%%
    rotmat = [sin(params(end-1)) cos(params(end-1)); -cos(params(end-1)) sin(params(end-1))];
    PD = rotmat * [0 .5]';
    contour(x,y,surface,5,'linewidth',2); hold on; % masked grid
    plot([0 PD(1)],[0 PD(2)],'k-')
    set(gca,'tickdir','out','XTick',ticks,'YTick',ticks,'xlim',[-axlim axlim],'ylim',[-axlim axlim]);
    
    %%% (3) Bubbles %%%
    maxnsp = max(max(GLMP.subunit{sub}.meannspikes));
    polar(0,axlim,'k*'); hold on; % this is stupid, just to get axis limits right
    scalar = 30;
    for i = 1:numel(GLMP.subunit{sub}.uniqueLcc)
        mn = meannsp(i)/maxnsp*scalar+scalar/10;
        h = polar(GLMP.subunit{sub}.uniquetheta(i),GLMP.subunit{sub}.uniquerho(i),'ko');
        if uniqL(i)
            set(h,'MarkerFaceColor','b','MarkerSize',mn)
        else
            set(h,'MarkerFaceColor','none','MarkerSize',mn)
        end
    end

    box on; grid off;
    axis equal square tight
    set(gca,'xlim',[-axlim axlim],'ylim',[-axlim axlim])

end


%% Compare models
maxsse = max(cat(1,cardirs.sse,coneiso.sse,cardfix.sse));

figure(1); clf; hold on;
plot(cardirs.sse,coneiso.sse,'ko')
set(gca,'xscale','log','yscale','log')
xlabel('cardinal free fit sse')
ylabel('cone iso free fit sse')
xlim([0 maxsse])
ylim([0 maxsse])
yyaxis right
plot(cardirs.sse,cardfix.sse,'ro')
set(gca,'yscale','log')
ylabel('cardinal fixed fit sse')
ylim([0 maxsse])

[t,p] = ttest(coneiso.sse,coneweights.sse)
nanmean(coneiso.sse-coneweights.sse);

[t,p] = ttest(cardirs.sse,coneiso.sse)
nanmean(cardirs.sse-coneiso.sse);

[t,p] = ttest(cardirs.sse,cardfix.sse)
nanmean(cardirs.sse-cardfix.sse);

% did all sse's come from the same dist? obvi no.
allsse = cat(1,coneweights.sse,coneiso.sse,cardirs.sse,cardfix.sse);
gr = repmat({'M1'},numel(coneweights.sse),1);
gr = cat(1,gr,repmat({'M2'},numel(coneweights.sse),1));
gr = cat(1,gr,repmat({'M3'},numel(coneweights.sse),1));
gr = cat(1,gr,repmat({'M4'},numel(coneweights.sse),1));
p = kruskalwallis(allsse,gr,'off')

%%
figure(2); clf; hold on; box on;
plot([0 1],[nanmean(coneweights.sse) nanmean(coneweights.sse)],'r-')
plot([0 1],[nanmean(coneiso.sse) nanmean(coneiso.sse)],'g-')
plot([0 1],[nanmean(cardirs.sse) nanmean(cardirs.sse)],'b-')
plot([0 1],[nanmean(cardfix.sse) nanmean(cardfix.sse)],'k-')
legend('Method 1','Method 2','Method 3','Method 4')
ylim([0 nanmean(coneweights.sse)*1.1])
set(gca,'xtick',[])
ylabel('mean SSE')

%% Plot of sse values across the 4 models

allsse = cat(1,coneweights.sse,coneiso.sse,cardirs.sse,cardfix.sse);
bins = linspace(min(allsse),max(allsse),50);
binsoffset = mean(diff(bins))/10;
figure(3); clf; hold on; box on;

h = histogram(coneweights.sse,bins,'DisplayStyle','stairs','edgecolor','r');
i = histogram(coneiso.sse,bins+binsoffset,'DisplayStyle','stairs','edgecolor','g');
j = histogram(cardirs.sse,bins+binsoffset*2,'DisplayStyle','stairs','edgecolor','b');
k = histogram(cardfix.sse,bins+binsoffset*3,'DisplayStyle','stairs','edgecolor','c');

maxval = max(cat(2,h.Values,i.Values,j.Values,k.Values));

plot([nanmean(coneweights.sse) nanmean(coneweights.sse)],[0 maxval],'r','linewidth',2)
plot([nanmean(coneiso.sse) nanmean(coneiso.sse)],[0 maxval],'g','linewidth',2)
plot([nanmean(cardirs.sse) nanmean(cardirs.sse)],[0 maxval],'b','linewidth',2)
plot([nanmean(cardfix.sse) nanmean(cardfix.sse)],[0 maxval],'c','linewidth',2)

xlim([bins(1) bins(end)+binsoffset*3])
ylim([0 maxval*1.01])

%% Make example figures for each distribution
idx = 135; % pan color * **

GLMP = GLMSPopData{idx,strcmp(datatypes,'GLMP')};
sub = GLMSPopData{idx,strcmp(datatypes,'Subunit')};

meannsp = GLMP.subunit{sub}.meannspikes;
%ticks = -.08:.04:.08; % RS axes
ticks = -.8:.2:.8; % regular axes
cmap = repmat(linspace(.7,0,32)',1,3);
popidx = find(poppanel.oneDL | poppanel.twoDL);
analidx = find(popidx==idx)-1;
params = coneweights.params(analidx,:);

% Define thetas and rhos for the whole grid
thetas = linspace(-pi,pi,361)';
majorax = max(GLMP.subunit{sub}.uniquerho(GLMP.subunit{sub}.uniquetheta==pi/4))*1.1;
minorax = max(GLMP.subunit{sub}.uniquerho(GLMP.subunit{sub}.uniquetheta==-pi/4))*1.1;
nom = majorax * minorax;
denom = (majorax*sin(thetas-pi/4)).^2 + (minorax*cos(thetas-pi/4)).^2;
rhos = nom ./ sqrt(denom);
scalars = 0:.01:1;
rhosgrid = rhos * scalars;
thetasgrid = repmat(thetas,1,numel(scalars));
[x,y] = pol2cart(thetasgrid,rhosgrid);
surface = ComputeNakaRushtonJPW(params,[x(:) y(:)],poppanel.surftype);
surface = reshape(surface,size(rhosgrid));
axlim = max(x(:));


%%% Bubbles %%%
figure(32); clf;
set(gcf,'units','normalized','pos',[.2 .2 .6 .6],'numbertitle','off',...
    'renderer','paint')
maxnsp = max(max(GLMP.subunit{sub}.meannspikes));

polar(0,axlim,'k*'); hold on; % this is stupid, just to get axis limits right
scalar = 30;
for i = 1:numel(GLMP.subunit{sub}.uniqueLcc)
    mn = meannsp(i)/maxnsp*scalar+scalar/10;
    h = polar(GLMP.subunit{sub}.uniquetheta(i),GLMP.subunit{sub}.uniquerho(i),'ko'); 
    set(h,'MarkerFaceColor','none','MarkerSize',mn)
end
set(gca,'XTick',[],'YTick',[],'xlim',[-axlim axlim],'ylim',[-axlim axlim]);

name = 'bubbles';
box on; grid off;
axis equal square tight
set(gca,'CameraPosition',[0 0 275.3526],'xlim',[-axlim axlim],'ylim',[-axlim axlim])
print(gcf,['/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/' name],'-depsc')


%%% Contours and PD %%%
figure(31); clf;  colormap(cmap);
set(gcf,'units','normalized','pos',[.2 .2 .6 .6],'numbertitle','off',...
    'renderer','paint')

rotmat = [sin(params(end-1)) cos(params(end-1)); -cos(params(end-1)) sin(params(end-1))];
PD = rotmat * [0 .5]';
contour(x,y,surface,5,'linewidth',2); hold on; % masked grid
plot([0 PD(1)],[0 PD(2)],'k-')
set(gca,'tickdir','out','XTick',ticks,'YTick',ticks,'xlim',[-axlim axlim],'ylim',[-axlim axlim]);

name = 'contours';
box on; grid off;
axis equal square tight
set(gca,'CameraPosition',[0 0 275.3526],'xlim',[-axlim axlim],'ylim',[-axlim axlim])
print(gcf,['/Users/jpatrickweller/Dropbox/GLMS Paper/Figures/' name],'-depsc')

