function [plat] = GLMP_ContrastSensitivity(plat)

% This code is intended to compare CRFs in every sampled direction

% 1/16/13   Created.    JPW


filename = plat{1}.datafile;

for p = 1:numel(plat)
    
    rots = unique(plat{p}.par.poltheta_orig);
    maxcon = max(plat{p}.par.polrho_orig);
    fitrots = [];
    color = [];
    
    % Some parameters for fitting
    vub = [1000 6 6 100];
    vlb = [0 0 0 0];
    options = optimset('MaxFunEvals',5000,'MaxIter',5000,'TolFun',10^-6,'TolX',10^-6,'Display','off','Algorithm','active-set');
    sigmaguess = 0.01;
    baselineguess = nanmean(plat{p}.trial.baselinensp);
    expguess = 2;
    fitL = ones(numel(rots),1);
    fitparams = nan(numel(rots),4,numel(plat));
    
    for r = 1:numel(rots)
        
        %Brifly set up figure
        figure(773+p); if r==1; clf; end;
        plotTitle = [filename ' Plat # ' num2str(p)];
        set(gcf,'Name',plotTitle,'NumberTitle','off')
        subplot(3,2,[1:2]); hold on; grid on;
        title('Contrast Sensitivity')
        xlabel('Cone Contrast')
        ylabel('Response (# of Spikes)')

        
        L = plat{p}.trial.poltheta_orig == rots(r);
        contrasts = plat{p}.trial.polrho_orig(L);
        responses = plat{p}.trial.nspikes(L);

        if numel(unique(contrasts)) < 3
            fitL(r) = 0;
            continue
        end
        
        %fitrots = cat(1,fitrots,rots(r));
        
        Aguess = max(responses);
        params0 = [Aguess sigmaguess expguess baselineguess];
        
        [fitparams(r,:,p),f] = fmincon('FitNakaRushtonFunJPW',params0,[],[],[],[],vlb,vub,[],options,contrasts,responses,'symmetric');
       
        % Plot Figure!
        hold on;
        x = linspace(0,maxcon,100);
        y = ComputeNakaRushtonJPW(fitparams(r,:,p),x,'symmetric');
        
        color = cat(1,color,unifrnd(0,1,1,3));
        h = plot(x,y,'--');
        set(h,'Color',color(end,:));
    
    end
        
    % Weed out untested directions
    fitL = logical(fitL);
    fitrots = rots(fitL);
    fitparams = fitparams(fitL,:,:);
    
    % Plot Raw Responses
    legend([num2str(fitrots./pi*180)],'Location','EastOutside')
    for r = 1:numel(fitrots)
        L = plat{p}.trial.poltheta_orig == fitrots(r);
        contrasts = plat{p}.trial.polrho_orig(L);
        responses = plat{p}.trial.nspikes(L);
        h = plot(contrasts,responses,'*'); hold on;
        set(h,'Color',color(r,:));
    end
    
    % Display parameters
    subplot(3,2,3); hold on; grid on;
    polar([fitrots; fitrots(1)],[fitparams(:,1,p); fitparams(1,1,p)],'o--')
    polar(0,0,'k*')
    title('Parameter: A')
    
    subplot(3,2,4); hold on; grid on;
    polar([fitrots; fitrots(1)],[fitparams(:,2,p); fitparams(1,2,p)],'o--')
    polar(0,0,'k*')
    title('Parameter: sigma')
    
    subplot(3,2,5); hold on; grid on;
    polar([fitrots; fitrots(1)],[fitparams(:,3,p); fitparams(1,3,p)],'o--')
    polar(0,0,'k*')
    title('Parameter: exponent')
    
    subplot(3,2,6); hold on; grid on;
    polar([fitrots; fitrots(1)],[fitparams(:,4,p); fitparams(1,4,p)],'o--')
    polar(0,0,'k*')
    title('Parameter: baseline')

    
    
end
