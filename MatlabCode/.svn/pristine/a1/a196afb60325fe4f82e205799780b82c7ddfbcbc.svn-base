% This code is a rhudamentary analysis of GridLMPlane datafiles
% 10/24/11  Created.    JPW
% 12/11     Modified.   JPW
%       Added radial grid.
% 1/12     Modified.   JPW
%       Added plateaus.
% 2/12      Modified.   JPW
%       Stimulus moves in two different directions.
% 3/2/12    Modified.   JPW
%       Added a new transMat to accomodate new stimulus gamut.
% 5/15/14   Modified.    JPW
%       Updated GLMP data structures to accomodate GLMS analyses.


clear all
close all

% GLMP v1.0
% Cartesian Grid
%datafile = 'S102111002.nex'; % Ellipse! (#1)
%datafile = 'S102111003.nex'; % Ellipse (#2)
%datafile = 'S102811002.nex'; % +/- Luminance (#3)
%datafile = 'S103111002.nex'; % - Luminance (#4)
%datafile = 'S103111003.nex'; % -L-M, with some +L+M (#5)
%datafile = 'S103111008.nex'; % +/- Luminance (#6)
%datafile = 'S110211007.nex'; % +/- Luminance  (Most Data) (#7)
%datafile = 'S110411012.nex'; % ellipse (check iso) (#8)
%datafile = 'S110911003.nex'; % +/- Luminance (#9) *
%datafile = 'S110911007.nex'; % +/- Luminance (#10) *
%datafile = 'S111011003.nex'; % -L-M, with some +L+M (Second to Most Data) (#11) *
%datafile = 'S112311006.nex'; % -L-M, with some +L+M (#12) *
%datafile = 'S112811003.nex'; % -lum (check iso) (#13)
%datafile = 'S113011002.nex'; % ellipse (check iso) (#14) *
%datafile = 'S120211008.nex';% -L+M (maybe just +L?) (#15) 
%datafile = 'S120511008.nex'; % +/- Luminance (#16)

% GLMP v2.0
% Radial Grid
%datafile = 'S120911003.nex';% ellipse? (check iso) (#17)
%datafile = 'S121211003.nex';% +/- Luminance (#18)
%datafile = 'S121211006.nex';% +/- Luminance (#19)
%datafile = 'S121511004.nex';% +/- Luminance (#20) *

% GLMP v2.1
% Plateaus: S-Cone
%datafile = 'S010212003.nex';% Elipse (check iso) (plat 1 is #21)
%datafile = 'S011312002.nex';% +/- Luminance (ellipse?) (plat 1 is #22)
%datafile = 'S020112002.nex';% +/- Luminance (plat 1 is #23)
%datafile = 'S020212007.nex';% -L-M, with some +L+M (plat 1 is #24) *
%datafile = 'S020812003.nex';% (not included in data list, check iso)
%datafile = 'S020812004.nex';% (not included in data list, check iso)

% GLMP v2.2
% Plateaus: 2 Directions of Motion
%datafile = 'S021012007.nex';% +L+M, with some -L-M (plat 1 #25)
%datafile = 'S021712004.nex';% (not included in data list, check iso)
%datafile = 'S021712006.nex';% -luminance (check iso) (plat 1 #26)
%datafile = 'S022112003.nex';% (not included in data list, check iso)
%datafile = 'S030112011.nex';% -luminance (check iso) (plat 1 #27)

% GLMP v3.2
% Plateaus: 2 Directions of Motion
%datafile = 'S032312002.nex';% ellipse? (check iso) (plat 1 #28)
%datafile = 'S032812004.nex';% Dropped 2 headers - must trim orignal file
%datafile = 'S032812006.nex';% Funny raw data... lookup in notebook
%datafile = 'S032912004.nex';% +/- Luminance (plat 2 is #29)
%datafile = 'S040412002.nex';% R-G Chromatic (plat 2 is #30) *
%datafile = 'S041012002.nex';% Ellipsoidal Cell!! (plat 2 #31)

% GLMP v3.3 
% Logarithmic spacing
%datafile = 'S041912004.nex'; % Luminance (#32)
%datafile = 'S042012002.nex'; % -L+M (#33)

% GLMP v3.0
% Symmetric, linear sampling
%datafile = 'A060112003.nex';% Ellispoidal Cell! (Apollo's First Cell!) (#34)
%datafile = 'A062612003.nex';% Ellispoidal Cell! (#35)
%datafile = 'A070312003.nex';% +/- Luminance (maybe ellipse) (#36)
%datafile = 'A070312005.nex';% Ellipsoidal Cell (#37)
%datafile = 'A071012004.nex';% Ellipsoidal (10 repeats) (#38)
%datafile = 'A071212005.nex';% Ellipsoidal (#39)
%datafile = 'A071612002.nex';% Ellipsoidal Cell (#40)
%datafile = 'A071612005.nex';% Ellipsoidal Cell (#41)

% GLMP v4.0 (Adaptive + Nonadaptive)
datafile = 'S101812003.nex'; % Luminance! (plat 1 is #42)
%datafile = 'S102312007.nex'; % +L-M Chromatic (plat 1 is #43)
%datafile = 'S102412003.nex'; % Ellipse! (plat 1 is #44)

% GLMP v3.4 (3 bar widths)
%datafile = 'S102412004.nex'; % Ellipse (plat 3 is #45)

% GLMP v4.0 (Adaptive + Nonadaptive)
%datafile = 'S110212005.nex'; % -L+M Chromatic (plat 1 is #46)



library = '/Users/jpatrickweller/Documents/Matlab_Working_Directory/Patrick/GridLMPlane/GLMPDatafiles/';
%library = 'Users/jpatrickweller/Dropbox/Patrick/DatafilesAll/';
%library = 'C:\Users\jpweller\Dropbox\Patrick\DatafilesAll\';

try
    rawdata = nex2stro([char(library) datafile]);
catch
    rawdata = nex2stro;
end


[trial par plat] = OrganizeRawGLMPData(rawdata);

% Plot subunits
for s = 1:numel(plat)
    
    if s == 1
        mfc = 'ro';
    elseif s == 2
        mfc = 'go';
    elseif s == 3
        mfc = 'ro';
    end
    
    markSize = (plat{s}.par.meanfr / max(plat{s}.par.meanfr)+.1)*20;
    %markSize = trial.fr;
    
    figure(279+s); clf;
    grid on; hold on;
    for i = 1:numel(plat{s}.par.Lcc_orig)
        h = plot(plat{s}.par.Lcc_orig(i),plat{s}.par.Mcc_orig(i),mfc,'Markersize',markSize(i));
    end
    axis equal tight;
    xlabel('Lcc')
    ylabel('Mcc')
    title('Mean Firing Rate (sp/s)')
    set(gcf,'Name',[num2str(datafile)],'NumberTitle','off')    

end

%%
set(gcf,'PaperPositionMode','auto')
filename = datafile(1:end-4);
print('-depsc',filename)


%% Figure(1): Raw Data and Fit Surfaces

% Removing 2nd and 3rd order polynomial from plot. 4/10/12 JPW

frMat = cat(1,par(:).meanfr);
tPtMat = cat(1,par(:).stim);

% For plotting
figure(1); clf; hold on;
xfit = -1:.1:1;
yfit = -1:.1:1;
[XFIT,YFIT] = meshgrid(xfit,yfit);

% Preallocate Space
%fitvals2OP = nan(numel(platvals),sum(platL));
%fitvals3OP = fitvals2OP;
%ssqErr2OP = nan(1,numel(platvals));
%ssqErr3OP = ssqErr2OP;

% Plot raw data
for n = 1:nplateaus
    
    % Plot each trial seperately
    if numel(platvals) > 1
        platL = platCatTr == platvals(n);
    else
        platL = logical(ones(size(trial.stim,1),1));
    end
        
    subplot(nplateaus,2,n*2-1); hold on; grid on; axis([-1.1 1.1 -1.1 1.1]);
    plot3(trial.stim(platL,1),trial.stim(platL,2),trial.fr(platL),'k*')
    title('All Raw Data')
    xlabel('L - M')
    ylabel('L + M')
    zlabel('Firing Rate')
    set(gca,'view',[-75 38]);
    
    % Plot averages of repeated stimuli
    if numel(platvals) > 1
        platL = platCatPar == platvals(n);
    else
        platL = logical(ones(size(tPtMat,1),1));
    end

    % Lollypop Plot
    subplot(nplateaus,2,n*2); hold on; grid on; axis([-1.1 1.1 -1.1 1.1]);
    stem3(tPtMat(platL,1),tPtMat(platL,2),frMat(platL),'fill','MarkerFaceColor','g')
    title('Mean Firing Rate by Stimulus')
    xlabel('L-M')
    ylabel('L+M')
    zlabel('Firing Rate')
    set(gca,'view',[-75 38]);
    
    % Second Order Polynomial
    expr1 = [ones(size(tPtMat(platL),1),1) tPtMat(platL,1:2) tPtMat(platL,1).^2 tPtMat(platL,2).^2 tPtMat(platL,1).*tPtMat(platL,2)];
    %weights = 1./cat(1,par(platL).varfrs);
    %weights(isinf(weights)) = median(weights);
    weights = ones(numel(frMat),1);
    b1 = lscov(expr1,frMat(platL),weights(platL));
    %[b1, dev1, stats1] = glmfit(expr1(:,2:end), frMat(platL), 'poisson', 'link', 'identity');
    ZFIT1 = b1(1) + b1(2)*XFIT + b1(3)*YFIT + b1(4)*XFIT.^2 + b1(5)*YFIT.^2 + b1(6)*XFIT.*YFIT;
    
    
    % Third Order Polynomial
    expr2 = [ones(size(tPtMat(platL),1),1) tPtMat(platL,1:2) tPtMat(platL,1).^2 tPtMat(platL,2).^2 tPtMat(platL,1).*tPtMat(platL,2) ...
        tPtMat(platL,1).^3 tPtMat(platL,2).^3 tPtMat(platL,1).^2.*tPtMat(platL,2) tPtMat(platL,1).*tPtMat(platL,2).^2];
    %weights = 1./cat(1,par(platL).varfrs);
    %weights(isinf(weights)) = median(weights);
    weights = ones(numel(frMat),1);
    %b2 = glmfit(expr2(:,2:end), frMat(platL), 'poisson', 'link', 'identity');
    b2 = lscov(expr2,frMat(platL),weights(platL));
    ZFIT2 = b2(1) + b2(2)*XFIT + b2(3)*YFIT + b2(4)*XFIT.^2 + b2(5)*YFIT.^2 + b2(6)*XFIT.*YFIT...
        + b2(7)*XFIT.^3 + b2(8)*YFIT.^3 + b2(9)*XFIT.^2.*YFIT + b2(10)*XFIT.*YFIT.^2;
    
%     % Plot
%     subplot(nplateaus,4,n*4-1); hold on; grid on; axis([-1.1 1.1 -1.1 1.1]);
%     plot3(tPtMat(platL,1),tPtMat(platL,2),frMat(platL),'k*')
%     mesh(XFIT,YFIT,ZFIT1)
%     xlabel('L - M')
%     ylabel('L + M')
%     zlabel('Firing Rate')
%     title('Second Order Polynomial Surface Fit to Data')
%     set(gca,'view',[-75 38]);
%     
%     subplot(nplateaus,4,n*4); hold on; grid on; axis([-1.1 1.1 -1.1 1.1]);
%     plot3(tPtMat(platL,1),tPtMat(platL,2),frMat(platL),'k*')
%     mesh(XFIT,YFIT,ZFIT2)
%     xlabel('L - M')
%     ylabel('L + M')
%     zlabel('Firing Rate')
%     title('Third Order Polynomial Surface Fit to Data')
%     set(gca,'view',[-75 38]);

    % Compute residuals
    fitvals2OP{n} = b1(1) + b1(2)*tPtMat(platL,1) + b1(3)*tPtMat(platL,2) + b1(4)*tPtMat(platL,1).^2 + b1(5)*tPtMat(platL,2).^2 + b1(6)*tPtMat(platL,1).*tPtMat(platL,2);
    ssqErr2OP(n) = sum((fitvals2OP{n} - frMat(platL)).^2);
    fitvals3OP{n} = b2(1) + b2(2)*tPtMat(platL,1) + b2(3)*tPtMat(platL,2) + b2(4)*tPtMat(platL,1).^2 + b2(5)*tPtMat(platL,2).^2 + b2(6)*tPtMat(platL,1).*tPtMat(platL,2)...
        + b2(7)*tPtMat(platL,1).^3 + b2(8)*tPtMat(platL,2).^3 + b2(9)*tPtMat(platL,1).^2.*tPtMat(platL,2) + b2(10)*tPtMat(platL,1).*tPtMat(platL,2).^2;
    ssqErr3OP(n) = sum((fitvals3OP{n} - frMat(platL)).^2);
    
end

%figure; hold on; grid on;
%plot3(tPtMat(platL,1),tPtMat(platL,2),fitvals2OP,'bo')
%plot3(tPtMat(platL,1),tPtMat(platL,2),fitvals3OP,'ro')
%plot3(tPtMat(:,1),tPtMat(:,2),frMat,'k*')


%% Figure(3): Firing Rate Over Time

%Mean FR over time
allfrs = cat(2,par.frs);
meanfrs = nanmean(allfrs,2);

figure(3); clf; 
subplot(2,1,1); hold on; grid on;
plot(1:size(trial.fr,1),trial.fr);
xlabel('Trial Number')
ylabel('Firing Rate')
title('Firing Rate Per Trial')

%Variance of FR over time
[orderedmeanfrs, orderedidx] = sort(cat(1,par.meanfr));
orderedvarfrs = cat(1,par(orderedidx).varfrs);

figure(3); subplot(2,1,2); hold on; grid on;
plot(orderedmeanfrs, orderedvarfrs,'o')
xlabel('Mean Firing Rates')
ylabel('Variance in Firing Rates')
title('Dependence of Variance on Mean Firing Rate')


%% Figure(4): Plot Rasters ***Under Construction***

%maxt = max(cat(1,trial.normspiketimes{:}));

% for n = 1:size(trial.stim,1)
%     spiketimes = round((trial.tspikes{n} - trial.stimon(n)) * 1000);
%     spiketimes = spiketimes(spiketimes > 0 & spiketimes < trial.duration(n)*1000);
%     deltafunc = zeros(1,round(trial.duration(n)*1000));
%     deltafunc(spiketimes) = 1;
%     %trial.normspiketimes{n} = spiketimes; %Units go from s to ms
% end  
    
% ***** Get this figure to plot all rasters stacked vertically *****   
    
%lengthSpike = 1/(size(trial.stim,1)*2);


% figure(5); clf; hold on; grid on;% axis([0 1 0 max(
%     plot(1:round(trial.duration(n)*1000),0,'k-')
%     for k = 1:numel(spiketimes)
%        plot([spiketimes(k) spiketimes(k)], [-.2 .2],'k-')
%     end
% end


%% Figure(5): PSTH

if nplateaus == 1
    
    % Segregate Points into Quadrants    
    pLpML = sign(trial.stim(abs(trial.stim(:,1)) > abs(trial.stim(:,2)),1))==1;
    mLmML = sign(trial.stim(abs(trial.stim(:,1)) > abs(trial.stim(:,2)),1))==-1;
    pLmML = sign(trial.stim(abs(trial.stim(:,1)) < abs(trial.stim(:,2)),2))==1;
    mLpML = sign(trial.stim(abs(trial.stim(:,1)) < abs(trial.stim(:,2)),2))==-1;
    
    allspikes = cat(1,trial.normspiketimes{:});
    
    bincenters = 0:5:max(allspikes);
    
    [allStimHist, binsloc] = hist(allspikes,bincenters);
    mLpMhist = hist(cat(1,trial.normspiketimes{mLpML}),bincenters);
    mLmMhist = hist(cat(1,trial.normspiketimes{mLmML}),bincenters);
    pLmMhist = hist(cat(1,trial.normspiketimes{pLmML}),bincenters);
    pLpMhist = hist(cat(1,trial.normspiketimes{pLpML}),bincenters);
    
    histmax = max([mLpMhist mLmMhist pLmMhist pLpMhist]);

    % Plot Figure
    figure(5);clf;
    subplot(5,1,1); hold on; grid on; axis([0 max(bincenters) 0 1]); axis 'auto y';
    title('All Responses')
    bar(binsloc,allStimHist,'k')
    subplot(5,1,2); hold on; grid on; axis([0 max(bincenters) 0 histmax]);
    title('Responses to -L+M Stimuli')
    bar(binsloc,mLpMhist,'k')
    subplot(5,1,3); hold on; grid on; axis([0 max(bincenters) 0 histmax]);
    title('Responses to -L-M Stimuli')
    ylabel('Number of Spikes')
    bar(binsloc,mLmMhist,'k')
    subplot(5,1,4); hold on; grid on; axis([0 max(bincenters) 0 histmax]);
    title('Responses to +L-M Stimuli')
    bar(binsloc,pLmMhist,'k')
    subplot(5,1,5); hold on; grid on; axis([0 max(bincenters) 0 histmax]);
    title('Responses to +L+M Stimuli')
    bar(binsloc,pLpMhist,'k')
    xlabel('Miliseconds from Stimulus Onset')
    
else
    
    for n = 1:nplateaus
        
        % Segregate Points into Quadrants
        platL = platCatTr == platvals(n);

        %pLpML = sign(trial.stim(abs(trial.stim(platL,1)) > abs(trial.stim(platL,2)),1))==1;
        %mLmML = sign(trial.stim(abs(trial.stim(platL,1)) > abs(trial.stim(platL,2)),1))==-1;
        %pLmML = sign(trial.stim(abs(trial.stim(platL,1)) < abs(trial.stim(platL,2)),2))==1;
        %mLpML = sign(trial.stim(abs(trial.stim(platL,1)) < abs(trial.stim(platL,2)),2))==-1;
        % This should correct for an asymmetric original space
        pLpML = sign(trial.stim(abs(trial.stim(platL,1)) > abs(trial.stim(platL,2)),1))==1;
        mLmML = sign(trial.stim(abs(trial.stim(platL,1)) > abs(trial.stim(platL,2)),1))==-1;
        pLmML = sign(trial.stim(abs(trial.stim(platL,1)) < abs(trial.stim(platL,2)),2))==1;
        mLpML = sign(trial.stim(abs(trial.stim(platL,1)) < abs(trial.stim(platL,2)),2))==-1;
        
        allspikes = cat(1,trial.normspiketimes{platL});
        
        bincenters = 0:5:max(allspikes);
        
        [allStimHist, binsloc] = hist(allspikes,bincenters);
        mLpMhist = hist(cat(1,trial.normspiketimes{mLpML}),bincenters);
        mLmMhist = hist(cat(1,trial.normspiketimes{mLmML}),bincenters);
        pLmMhist = hist(cat(1,trial.normspiketimes{pLmML}),bincenters);
        pLpMhist = hist(cat(1,trial.normspiketimes{pLpML}),bincenters);
        
        histmax = max([mLpMhist mLmMhist pLmMhist pLpMhist]);
        
        % Plot Figure
        figure(5); if n==1; clf; end;
        subplot(5,nplateaus,n); hold on; grid on; axis([0 max(bincenters) 0 1]); axis 'auto y';
        title('All Responses')
        bar(binsloc,allStimHist,'k')
        subplot(5,nplateaus,n+nplateaus); hold on; grid on; axis([0 max(bincenters) 0 histmax]);
        title('Responses to -L+M Stimuli')
        bar(binsloc,mLpMhist,'k')
        subplot(5,nplateaus,n+2*nplateaus); hold on; grid on; axis([0 max(bincenters) 0 histmax]);
        title('Responses to -L-M Stimuli')
        ylabel('Number of Spikes')
        bar(binsloc,mLmMhist,'k')
        subplot(5,nplateaus,n+3*nplateaus); hold on; grid on; axis([0 max(bincenters) 0 histmax]);
        title('Responses to +L-M Stimuli')
        bar(binsloc,pLmMhist,'k')
        subplot(5,nplateaus,n+4*nplateaus); hold on; grid on; axis([0 max(bincenters) 0 histmax]);
        title('Responses to +L+M Stimuli')
        bar(binsloc,pLpMhist,'k')
        xlabel('Miliseconds from Stimulus Onset')
        
    end
    
end


%% Figure(8): Lollipop-style plots, but counting spikes in narrow time windows that go from stimon to stimoff.

for q = 1:nplateaus
    
    
    allcastim = cat(1,par.stim);
    %platL = allcastim(:,3) == platvals(q);
    platL = platCatPar == platvals(q);
    platidx = find(platL);
    
    tbins = linspace(0,max(cat(1,trial.normspiketimes{:})),9);
    nsubplots = numel(tbins)-1;
    
    for k = 2:numel(tbins)
        %for n = 1:size(par,2)
        for n = 1:numel(platidx)   
            %par(n).timehist(k-1) = numel(par(n).catspks(par(n).catspks > tbins(k-1) & par(n).catspks <= tbins(k)));
            par(platidx(n)).timehist(k-1) = numel(par(platidx(n)).catspks(par(platidx(n)).catspks > tbins(k-1) & par(platidx(n)).catspks <= tbins(k)));
        end
    end
    tHistMat = cat(1,par(platL).timehist);
    
    % Plot figure
    if q == 1; figure(8);clf; hold on; end
    for n = 1:nsubplots
        subplot(nplateaus,nsubplots,(q-1)*nsubplots+n); hold on; grid on;
        %subplot(1,nsubplots,n); hold on; grid on;
        stem3(tPtMat(platL,1),tPtMat(platL,2),tHistMat(:,n),'fill','MarkerFaceColor','g')
        set(gca,'view',[-75 38]); axis([-1 1 -1 1 0 max(tHistMat(:))]);
        xlabel('L-M')
        ylabel('L+M')
        zlabel('Firing Rate')
        title(['t = ', num2str(round(tbins(n))),'ms to ',num2str(round(tbins(n+1))),'ms'])
    end
    
end


%% Figure(9): Wedge Analysis

% Choice angle between 0 and 180
choiceAng = 180;
% Theta is added to either side of the choiceAng
theta = 5;

angMat = rnddec(((atan2(tPtMat(:,2),tPtMat(:,1))./(pi))*180),4);
angMat(sign(angMat)==-1) = (180 - abs(angMat(sign(angMat)==-1)))+180;

for q = 1:nplateaus
    
    angMatPlat = angMat(platCatPar == platvals(q));    
    if sign(choiceAng-theta)==1 && sign(choiceAng+theta)==1
        angL = angMatPlat >= (choiceAng - theta) & angMatPlat <= (choiceAng + theta);
    elseif sign(choiceAng-theta)==-1 && sign(choiceAng+theta)==1
        angL = angMatPlat >= (choiceAng - theta + 360) | angMatPlat <= (choiceAng + theta);
    elseif sign(choiceAng-theta)==1 && choiceAng+theta>360
        angL = angMatPlat >= (choiceAng - theta) & angMatPlat <= (choiceAng + theta - 360);
    elseif sign(choiceAng-theta)==-1 && sign(choiceAng+theta)==-1
        disp('Wedge is impossible')
    end
    tPtMatPlat = tPtMat(platCatPar == platvals(q),:);
    frMatPlat = frMat(platCatPar == platvals(q));
    frWedgePlat = frMat(angL);
    tPtWedgePlat = tPtMat(angL,:);
    angWedgePlat = angMat(angL);
    
    % Radial Contrast Response Function
    radCon = nan(numel(frWedgePlat),1);
    for n = 1:size(tPtWedgePlat,1)
        radCon(n) = vectlength([0 0 0],tPtWedgePlat(n,:));
    end
    
    % Projection Contrast Response Function
    thetas = angWedgePlat - choiceAng;
    vectProjs = radCon .* cosd(thetas);
    
    %Plot Results
    figure(9); if q==1;clf;end; hold on;
    subplot(nplateaus,3,q*3-2); hold on; grid on; axis([-1 1 -1 1]);
    stem3(tPtMatPlat(~angL,1),tPtMatPlat(~angL,2),frMatPlat(~angL))
    stem3(tPtMatPlat(angL,1),tPtMatPlat(angL,2),frWedgePlat,'fill','MarkerFaceColor','g','MarkerEdgeColor','g','Color','g')
    plot([0 cosd(choiceAng)],[0 sind(choiceAng)],'g')
    plot([0 cosd(choiceAng-theta)],[0 sind(choiceAng-theta)],'g--')
    plot([0 cosd(choiceAng+theta)],[0 sind(choiceAng+theta)],'g--')
    set(gca,'view',[-75 38]);
    title('Firing Rate by Stimulus')
    xlabel('L-M')
    ylabel('L+M')
    zlabel('Firing Rate')
    subplot(nplateaus,3,q*3-1); hold on; grid on;
    plot(radCon,frWedgePlat,'o')
    title('Radial Contrast')
    xlabel('Radial Contrast')
    ylabel('Firing Rate (sp/s)')
    subplot(nplateaus,3,q*3); hold on; grid on;
    plot(vectProjs,frWedgePlat,'o')
    title('Projection Contrast')
    xlabel('Projected Contrast')
    ylabel('Firing Rate (sp/s)')
    
end


%% Gregs PSTH

% 
% % A couple of analyses for Patrick's experiment
% % switching over the B-splines
% 
% % Greg's labeling
% %stro = nex2stro(findfile('S121211003'));
% %stro = nex2stro('Users/jpatrickweller/MATLAB/GridLMPlane/Datafiles/S121211003.nex');
% stro = rawdata;
% Lcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'Lcc'));
% Mcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'Mcc'));
% stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimon_t'));
% stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fpoff_t'));
% ntrials = size(stro.trial,1);
% LM = Lcc+Mcc;
% LvM = 10*(Lcc-Mcc);
% 
% 
% 
% bins = linspace(0,2,50);
% PSTH = zeros(1,length(bins));
% % First PSTHs to find good values for "offset"
% for i = 1:ntrials
%     spiketimes = stro.ras{i,1}-stimon_t(i);
%     PSTH = PSTH +histc(spiketimes',bins);
% end
% 
% figure; axes; hold on;
% bar(bins,PSTH);
% gauss = @(x,params)(params(2)+params(1)*exp((-(x-params(3)).^2)/params(4).^2));
% bguess = mean(PSTH(bins<.25));
% aguess = max(PSTH)-bguess;
% muguess = bins(find(PSTH == max(PSTH),1,'first'));
% sigmaguess = .1; % terrible. Need something better.
% %plot(bins,gauss(bins,[aguess,bguess,muguess,sigmaguess]),'m.')
% fittedparams = fminsearch(@(params)sum((gauss(bins,params)-PSTH).^2),[aguess,bguess,muguess,sigmaguess]);
% plot(bins,gauss(bins,fittedparams),'b-.','LineWidth',3)
% offset = [fittedparams(3)-fittedparams(4) fittedparams(3)+fittedparams(4)];
% plot(offset,[fittedparams(1) fittedparams(1)],'k-','linewidth',3);
% title('Gregs PSTH (Check with mine)')
% 
% % Now counting up spikes in a window
% ntrials = size(stro.trial,1);
% nspikes = zeros(ntrials,1);
% for i = 1:ntrials
%     nspikes(i) = sum(stro.ras{i,1} > stimon_t(i)+offset(1) & stro.ras{i,1}<stimon_t(i)+offset(2));
% end


%% Figure(10): 1-D Spline Fit

% Projecting onto a few unit vectors in the LM plane and
% Fitting the 1-D function (response vs. projection) with
% regular old 1-D cubic smoothing spline.

nthetas = 1000;
thetas = linspace(0,pi-pi/nthetas,nthetas);
data = nan*ones(nthetas,1);
smoothingparam = .9;

figure(10); clf; hold on;

for q = 1:nplateaus;
    
    platL = platCatTr == platvals(q);
    subplot(nplateaus,2,nplateaus*(q-1)+1); hold on;
    if q == 1; title('Spline fit to Projections'); end;
    ylabel('Spikes');
    if q == nplateaus; xlabel('Contrast'); end
    
    for i = 1:length(thetas);
        unitvect = [cos(thetas(i)) sin(thetas(i))]';
        %projs = [LM LvM]*unitvect;
        projs = [trial.stim(platL,1) trial.stim(platL,2)] * unitvect;
        lastwarn('');
        %pp  = csaps(projs,nspikes,smoothingparam);
        pp  = csaps(projs,trial.nspikes(platL),smoothingparam);
        %plot(projs,nspikes,'k.');
        figure(10); subplot(nplateaus,2,nplateaus*(q-1)+1); hold on;
        fnplt(pp);
        plot(projs,trial.nspikes(platL),'k.');
        drawnow;
        %data(i) = sum((fnval(pp,projs)-nspikes).^2);
        data(i) = sum((fnval(pp,projs)-trial.nspikes(platL)).^2);
        if ~isempty(lastwarn)
            data(i) = nan;
        end
        cla
    end
    
    bp = find(data==min(data));
    preftheta = thetas(bp);
    unitvect = [cos(thetas(bp)) sin(thetas(bp))]';
    prefproj = [trial.stim(platL,1) trial.stim(platL,2)] * unitvect;
    prefpp = csaps(prefproj,trial.nspikes(platL),smoothingparam);
    fitvals1Dsp{q} = fnval(prefpp,unique(prefproj));
    fnplt(prefpp);
    plot(prefproj,trial.nspikes(platL),'k.');
    
    % Plot "goodness of fits"
    subplot(nplateaus,2,nplateaus*(q-1)+2); axis([0 180 0 max(data)*1.1]) ; hold on; grid on;
    plot(linspace(0,180-180/nthetas,nthetas),data);
    xlabel('Theta (Deg)'); ylabel('SSE from Spline Fit');
    title('Best Projection Vector')
    plot(preftheta/pi*180,min(data),'*r')
    
end


%% Figure(11): Fitted Value vs Residual

for n = 1:numel(platvals)
    
    platL = platCatPar == platvals(n);
    
    %[fitvals2OP_sorted, fitvals2OP_sortidx] = sort(fitvals2OP)
    resid2OP{n} = frMat(platL) - fitvals2OP{n};
    resid3OP{n} = frMat(platL) - fitvals3OP{n};
    resid1Dsp{n} = frMat(platL) - fitvals1Dsp{n};
    
    SSEresid2OP(n) = sum(resid2OP{n}.^2);
    SSEresid3OP(n) = sum(resid3OP{n}.^2);
    SSEresid1Dsp(n) = sum(resid1Dsp{n}.^2);
    
    maxresid(n) = max([resid2OP{n}' resid3OP{n}' resid1Dsp{n}']);
    minresid(n) = min([resid2OP{n}' resid3OP{n}' resid1Dsp{n}']);
    maxfit(n) = max([fitvals2OP{n}' fitvals3OP{n}' fitvals1Dsp{n}']);
    
    % Plot figure
    figure(11);
    %title('Model Goodness of Fit')
    ylabel('Residuals')
    subplot(3,nplateaus,n); hold on; grid on; axis square equal;% axis([0 maxfit(n)*1.1 minresid(n)*1.1 maxresid(n)*1.1])
    title('Second Order Polynomial')
    ylabel('Residuals')
    plot(fitvals2OP{n},resid2OP{n},'*k')
    subplot(3,nplateaus,n+2); hold on; grid on; axis square equal;% axis([0 maxfit(n)*1.1 minresid(n)*1.1 maxresid(n)*1.1])
    title('Third Order Polynomial')
    ylabel('Residuals')
    plot(fitvals3OP{n},resid3OP{n},'*b')
    subplot(3,nplateaus,n+4); hold on; grid on;axis square equal;% axis([0 maxfit(n)*1.1 minresid(n)*1.1 maxresid(n)*1.1])
    title('1-D Spline')
    plot(fitvals1Dsp{n}, resid1Dsp{n},'*r')
    xlabel('Fitted Value')
    ylabel('Residuals')

end


%% Figure(12): Using L+M and L-M to Predict Test Points
figure(12); clf; cnt = 1;
for n = 1:nplateaus
    
    platL = platCatPar == platvals(n); 
    
    [theta rho] = cart2pol(tPtMat(platL,1),tPtMat(platL,2));
    thetas = theta ./ pi * 180 + 180;
    rotAngs = unique(thetas(thetas>0 & thetas<=90));
    platpts = tPtMat(platL,1:2);
    platfrs = frMat(platL);
        
    for q = 1:numel(rotAngs)
    
        % Define Principle axes
        [LpMx LpMy] = pol2cart(rotAngs(q)/180*pi,1)
        LpMx = rndofferr(LpMx,4);
        LpMy = rndofferr(LpMy,4);
        [LmMx LmMy] = pol2cart((rotAngs(q)+90)/180*pi,1)
        LmMx = rndofferr(LmMx,4);
        LmMy = rndofferr(LmMy,4);
        
        % Contrasts along axes
%         LpML = angMat == rotAngs(q) | angMat == rotAngs(q)+180 & platL;
%         LmML = angMat == rotAngs(q)+90 | angMat == rotAngs(q)+270 & platL;
%         LpMpts = tPtMat(LpML,1:2);
%         LmMpts = tPtMat(LmML,1:2);
%         LpMcont = nan(sum(LpML),1);
%         LmMcont = nan(sum(LmML),1);
%         for s = 1:sum(LpML)  
%             LpMcont(s) = norm(LpMpts(s,:));
%         end
%         for s = 1:sum(LmML)
%             LmMcont(s) = norm(LmMpts(s,:));
%         end
        
        LpMcont = rho(thetas==rotAngs(q) | thetas==rotAngs(q)+180)
        LmMcont = rho(thetas==rotAngs(q)+90 | thetas==rotAngs(q)+270)
        
        % Responses at Contrasts
        %LpMresp = frMat(LpML);
        %LmMresp = frMat(LmML);
        LpMresp = platfrs(thetas==rotAngs(q) | thetas==rotAngs(q)+180)
        LmMresp = platfrs(thetas==rotAngs(q)+90 | thetas==rotAngs(q)+270)

        
        % Fit spline
        [LpMpp smoothLpM] = csaps(LpMcont,LpMresp,1);
        [LmMpp smoothLmM] = csaps(LmMcont,LmMresp,1);
        
        % Project test points onto principle axes
        %LpMax = [0 1];
        %LmMax = [1 0];
        LpMax = [LpMx LpMy];
        LmMax = [LmMx LmMy];
        %projLpMax = LpMax * tPtMat(:,1:2)';
        %projLmMax = LmMax * tPtMat(:,1:2)';
        projLpMax = LpMax * platpts'
        projLmMax = LmMax * platpts'
        
        predLpMax = fnval(LpMpp,projLpMax);
        predLmMax = fnval(LmMpp,projLmMax);
        
        % Predict firing rates
        %predictions = cosd(angMat).^2 .* predLmMax' + sind(angMat).^2 .* predLpMax';
        predictions = cosd(thetas).^2 .* predLmMax' + sind(thetas).^2 .* predLpMax';
        
        % Plot predictions
        subplot(numel(rotAngs),nplateaus,cnt);hold on;cnt=cnt+1;
        plot(platfrs,predictions,'o')
        %plot(platfrs,predLpMax,'go')
        xlabel('Measured Mean Firing Rate (sp/s)')
        ylabel('Predicted Firing Rate (sp/s)')
        %legend('Two Principle Axes','One Principle Axis')
        legend(['Theta = ' num2str(rotAngs(q))])
        axis equal 
        
    end
end

%%

% % % Taking a look at the residuals from a 1-D fit
% % unitvect = [cos(preftheta) sin(preftheta)]';
% % projs = [LM LvM]*unitvect;
% % pp  = csaps(projs,nspikes,smoothingparam);
% % predfr1D = fnval(pp,projs);
% % figure; 
% % plot(predfr1D,nspikes-predfr1D,'k.');
% % xlabel('fitted value'); ylabel('residual');
% % [x,y]= meshgrid(linspace(-1,1,50),linspace(-1,1,50));
% % projs = [x(:) y(:)]*unitvect;
% % predfr2D = fnval(pp,projs);
% % figure; axes; hold on; axis square;
% % surf(x,y,reshape(predfr2D,size(x)));
% % plot3(LM,LvM,nspikes,'k.');
% % LMconeweights = 1/sqrt(2)*[1 1;1 -1]*unitvect;
% % title(['L weight = ',num2str(LMconeweights(1),1),' M weight = ',num2str(LMconeweights(2),1)]);
% % 
% 


%%
% Trying a thinplate spline on the sqrt-transformed spikecounts.
stro = rawdata;
Lcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'Lcc'));
Mcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'Mcc'));
stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimon_t'));
stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'fpoff_t'));
ntrials = size(stro.trial,1);
LM = Lcc+Mcc;
LvM = 10*(Lcc-Mcc);
nspikes = trial.nspikes;


figure; axes; hold on;
out2d = tpaps([LM, LvM]', nspikes');
plot3(LM,LvM,nspikes,'k.');
h = fnplt(out2d);
surf(h{1},h{2},h{3})

%out2d = tpaps([LM, LvM]', sqrt(nspikes'));
%plot3(LM,LvM,sqrt(nspikes),'k.');
%h = fnplt(out2d);
%surf(h{1},h{2},h{3}.^2)
axis square;

out2d_noxform = tpaps([LM, LvM]', nspikes');
h = fnplt(out2d_noxform);
h =surf(h{1},h{2},h{3},repmat(5,size(h{1})))
set(h,'FaceAlpha',.7)

%The square root transformation looks like it's pretty subtle.

%Now some 1-D splines
outx = tpaps([LM, repmat(mean(LvM),size(LvM,1),1)]', sqrt(nspikes)',1);
outx = csaps(LM, sqrt(nspikes));
outy = csaps(LvM, sqrt(nspikes));
figure;
subplot(2,1,1);
fnplt(outx);
subplot(2,1,2);
fnplt(outy);


% %%
% % % Bootstrapping test to determine whether the 2-D fit is actually helping at all.
% % nboot = 2000;
% % PARAMETRIC = 2; % 0=Nonparametric, 1=Poisson, 2=Gaussian
% % unitvect = [cos(preftheta) sin(preftheta)]';
% % projs = [LM LvM]*unitvect;
% % %projs = [LM LvM]*[unitvect(2) unitvect(1)]';
% % pp1D  = csaps(projs,nspikes,smoothingparam);
% % predfr1D = fnval(pp1D,projs);
% % resid = predfr1D-nspikes;
% % SS1D = sum(resid.^2);
% % SS2D = sum((fnval(out2d,[LM LvM]')-nspikes').^2);
% % [u,a,b] = unique([LM LvM],'rows');
% % SS = zeros(nboot,2);
% % for i = 1:nboot
% %      i
% %     tmpdata = zeros(size(u,1),1);
% %     for j = 1:max(b)    
% %         L = b == j;
% %         r = resid(L);
% %         if (PARAMETRIC == 0)
% %             tmpdata(L) = predfr1D(L) + r(unidrnd(length(r),length(r),1)); % Nonparametric
% %         elseif (PARAMETRIC == 1)
% %             tmpdata(L) = poissrnd(predfr1D(L)); % Parametric (Poisson)
% %         else
% %             tmpdata(L) = predfr1D(L) + normrnd(mean(r),std(r),size(r)); % Parametric (Gaussian)
% %         end
% %     end
% %     pp  = csaps(projs,tmpdata,smoothingparam);
% %     predfr = fnval(pp,projs);
% %     SS(i,1) = sum((predfr-tmpdata).^2);
% %     out = tpaps([LM, LvM]', tmpdata');
% %     SS(i,2) = sum((fnval(out,[LM LvM]')-tmpdata').^2);
% % end
% % SSratio = SS(:,1)./SS(:,2);  % Large values mean that TP spline fits way better
% % figure; axes; hold on;
% % hist(SSratio);
% % plot(SS1D./SS2D,0,'m*');
% % p = sum((SS1D./SS2D)<SSratio)/nboot;
% % title(['p = ',num2str(p)]);
% % 
% %

%%
% Fitting individual contrast-response functions by least-squares
% This only works for cells that were tested on a radial grid 
thetas = atan2(LvM, LM);
bins = linspace(-pi,pi,1000);
n = histc(thetas,bins);
uniquethetas = bins(n>(max(n)/2));
fittedparams = [];
figure; axes; hold on;
options = optimset;
options = optimset(options,'Diagnostics','off','Display','off');
options = optimset(options,'LargeScale','off');
vlb = [0  0.001 0.001 0];
vub = [1000 1000 100 100];
for i = 1:length(uniquethetas)-1
    unitvect = [cos(uniquethetas(i)) sin(uniquethetas(i))]';
    projs = [LM LvM]*unitvect;
    L = thetas > uniquethetas(i) & thetas < uniquethetas(i+1);
    plot(projs(L),nspikes(L),'k.')
    % Using the thin-plate spline to derive good initial guesses
    params0 = [max(nspikes(L)), prctile(abs(projs(L)),1), 1, 0];
    params1 = fmincon('MyFitNakaRushtonFun',params0,[],[],[],[],vlb,vub,[],options,abs(projs(L)),fnval(out2d,[LM(L) LvM(L)]')'); 
    f = fmincon('MyFitNakaRushtonFun',params1,[],[],[],[],vlb,vub,[],options,abs(projs(L)),nspikes(L)); 
    
    fittedparams(i,:) = f;
    x = linspace(0,max(projs(L)),100);
    y = MyComputeNakaRushton(f,x);
    plot(x,y,'b-');
    title(num2str(uniquethetas(i)));
    drawnow;    
    pause(.5);
    cla;
end

figure;
for i = 1:size(fittedparams,2)
    subplot(3,1,i);
    polar(uniquethetas(1:end-1),fittedparams(:,i)','k.-')
end


%%
% Now fitting a single naka-rushton to the full data set, assuming L+M tuning
options = optimset;
options = optimset(options,'Diagnostics','off','Display','off');
options = optimset(options,'LargeScale','off');
vlb = [0  0.001 0.001 0];
vub = [1000 1000 100 100];

params0 = fittedparams(find(uniquethetas<0,1,'last'),:); % Using the L+M spoke as a first guess

f = fmincon('MyFitNakaRushtonFun',params0,[],[],[],[],vlb,vub,[],options,abs(LM),nspikes);
% Note that in the above line we're fitting to the absolute values of LM

% Plotting the fit
figure; axes; hold on;
plot3(LM,LvM,nspikes,'k.');
x = linspace(min(LM),max(LM),100);
z = MyComputeNakaRushton(f,abs(x));
h = surf(repmat(x,length(x),1),repmat(x,length(x),1)',repmat(z,length(z),1))
set(h,'EdgeAlpha',0);

% %%
% Fitting a Naka-Rushton to the full data set, allowing different
% A, sigma, and n for different luminance polarities
params0 = [f(1) f(1) f(2) f(2) f(3) f(3) f(4)];
vlb = [0  0 0.001 0.001 0.001 0.001 0 0];
vub = [1000 1000 1000 1000 100 100 100 100];

ff = fmincon('MyFitNakaRushtonFun',params0,[],[],[],[],vlb,vub,[],options,LM,nspikes);
% Plotting the fit
figure; axes; hold on;
plot3(LM,LvM,nspikes,'k.');
x = linspace(min(LM),max(LM),100);
z = MyComputeNakaRushton(ff,x);
h = surf(repmat(x,length(x),1),repmat(x,length(x),1)',repmat(z,length(z),1))
set(h,'EdgeAlpha',0);

