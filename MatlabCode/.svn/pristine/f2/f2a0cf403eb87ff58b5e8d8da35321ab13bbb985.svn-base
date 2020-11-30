function [plat] = GLMP_TimingAnalysis(plat,ver)

% This code is intended to analyze the arrival times of signals to various
% L/M combinations.

% 1/14/13   Created.    JPW
% 1/22/13   Modified (Rosters separated by trial).     JPW
% 1/24/13   Modified (Added second version for contrast display).   JPW


filename = plat{1}.datafile;

if nargin < 2
    ver = 2;
end


%% Plot rasters by color direction

if ver == 1
    
   for p = 1:numel(plat)
        
        % Set up variables
        mindiffrots = min(diff(unique(plat{p}.trial.poltheta_orig)))/pi*180;
        %maxdur = max(plat{p}.trial.stimdur);
        maxdur = plat{p}.GaussianSpikingProfile.mu + plat{p}.GaussianSpikingProfile.sigma*2;
        maxsamp = max(plat{p}.par.nsamps);
        spacing = linspace(-mindiffrots/2,mindiffrots/2,maxsamp*2+1);
        spikelines = spacing(2:2:end);
        spikeheight = diff(spacing(1:2))*.45;
        angs4plot = [];
        
        % Set up figure
        figure(653+p); clf; hold on;
        plotTitle = [filename ' Plat # ' num2str(p)];
        set(gcf,'Name',plotTitle,'NumberTitle','off')
        title('Spike Rasters by Contrast (Low to High)')
        xlabel('Time (s) from Stimulus Onset')
        ylabel('Color Direction (deg)')
        xlim([0 maxdur])
        
        plot([0 maxdur],[-135 -135],'y--')
        plot([0 maxdur],[-90 -90],'g--')
        plot([0 maxdur],[-45 -45],'m--')
        plot([0 maxdur],[0 0],'r--')
        legend('L+M','M','L-M','L')
        
        % Minor correction for error... this SHOULD be done in OrgRawData...
        plat{p}.trial.polrho_orig = rndofferr(plat{p}.trial.polrho_orig,3);
        
        allrhos = sort(unique(plat{p}.trial.polrho_orig),'descend');
        
        % Loop through contrasts
        for r = 1:2%:numel(allrhos)
            
            Lrho = plat{p}.trial.polrho_orig == allrhos(r);
            rots = unique(plat{p}.trial.poltheta_orig(Lrho));
            
            % Loop through angles
            for th = 1:numel(rots)
                
                Lrot = plat{p}.trial.polrho_orig == allrhos(r) & plat{p}.trial.poltheta_orig == rots(th);
                nsamp = sum(Lrot);
                trnum = find(Lrot);
                
                % For plotting
                tempang = rots(th)/pi*180 + 360*(r-1);
                
                % Loop through samples
                for s = 1:nsamp
                    
                    if s == ceil(nsamp/2)
                        if rots(th) == 0 || rots(th) == pi
                            plot([0 maxdur],[tempang+spikelines(s) tempang+spikelines(s)],'r--')
                            %set(gca,'LineWidth',.9)
                        elseif rots(th) == -pi/2 || rots(th) == pi/2
                            plot([0 maxdur],[tempang+spikelines(s) tempang+spikelines(s)],'g--')
                        elseif rots(th) == -3*pi/4 || rots(th) == pi/4
                            plot([0 maxdur],[tempang+spikelines(s) tempang+spikelines(s)],'y--')
                        elseif rots(th) == 3*pi/4 || rots(th) == -pi/4
                            plot([0 maxdur],[tempang+spikelines(s) tempang+spikelines(s)],'m--')
                        end
                    end
                    
                    spikes = plat{p}.trial.normtspikes{trnum(s)};
                    plot([spikes spikes]',repmat([spikelines(s)+tempang+spikeheight spikelines(s)+tempang-spikeheight],numel(spikes),1)','k')
                    %plot([0 maxdur],[spikelines(s)+tempang spikelines(s)+tempang],'b--')
                    
                end
                
                plot([0 maxdur],[spacing(1)+tempang spacing(1)+tempang],'k')
                plot([0 maxdur],[spacing(end)+tempang spacing(end)+tempang],'k')
                
                angs4plot = cat(1,angs4plot,tempang);
                
            end
            
        end
        
        % Plotting stuff
        ylim([min(angs4plot)-mindiffrots/2 max(angs4plot)+mindiffrots/2])
        set(gca,'YTick',-135:45:360*numel(allrhos)-180)
        set(gca,'YTickLabel',-135:45:180)
        set(gca,'xlim',[maxdur - plat{p}.GaussianSpikingProfile.sigma*4 maxdur]);
        
   end
    
    
elseif ver == 2
    
    for p = 1:numel(plat)
        
        % Set up variables
        mindiffrots = min(diff(unique(plat{p}.trial.poltheta_orig)))/pi*180;
        %maxdur = max(plat{p}.trial.stimdur);
        maxdur = plat{p}.GaussianSpikingProfile.mu + plat{p}.GaussianSpikingProfile.sigma*2;
        maxsamp = max(plat{p}.par.nsamps);
        spacing = linspace(-mindiffrots/2,mindiffrots/2,maxsamp*2+1);
        spikelines = spacing(2:2:end);
        spikeheight = diff(spacing(1:2))*.45;
        
        % Set up figure
        figure(653+p); clf; hold on;

        % Minor correction for error... this SHOULD be done in OrgRawData...
        plat{p}.trial.polrho_orig = rndofferr(plat{p}.trial.polrho_orig,3);
        
        allrhos = sort(unique(plat{p}.trial.polrho_orig),'descend');
        
        % Loop through contrasts
        for r = 1%:numel(allrhos)
            
            Lrho = plat{p}.trial.polrho_orig == allrhos(r);
            rots = unique(plat{p}.trial.poltheta_orig(Lrho));
            subplot(numel(allrhos),1,r); hold on;
            title(['Contrast = ' num2str(allrhos(r))])
            if r ==1
                plot([0 maxdur],[-135 -135],'y--')
                plot([0 maxdur],[-90 -90],'g--')
                plot([0 maxdur],[-45 -45],'m--')
                plot([0 maxdur],[0 0],'r--')
                legend('L+M','M','L-M','L')
            end
            
            % Loop through angles
            for th = 1:numel(rots)
                
                Lrot = plat{p}.trial.polrho_orig == allrhos(r) & plat{p}.trial.poltheta_orig == rots(th);
                nsamp = sum(Lrot);
                trnum = find(Lrot);
                tempang = rots(th)/pi*180;
                
                % Loop through samples
                for s = 1:nsamp
                    
                     if s == ceil(nsamp/2)
                        if rots(th) == 0 || rots(th) == pi
                            plot([0 maxdur],[tempang+spikelines(s) tempang+spikelines(s)],'r--')
                            %set(gca,'LineWidth',.9)
                        elseif rots(th) == -pi/2 || rots(th) == pi/2
                            plot([0 maxdur],[tempang+spikelines(s) tempang+spikelines(s)],'g--')
                        elseif rots(th) == -3*pi/4 || rots(th) == pi/4
                            plot([0 maxdur],[tempang+spikelines(s) tempang+spikelines(s)],'y--')
                        elseif rots(th) == 3*pi/4 || rots(th) == -pi/4
                            plot([0 maxdur],[tempang+spikelines(s) tempang+spikelines(s)],'m--')
                        end
                    end
                    
                    spikes = plat{p}.trial.normtspikes{trnum(s)};
                    plot([spikes spikes]',repmat([spikelines(s)+tempang+spikeheight spikelines(s)+tempang-spikeheight],numel(spikes),1)','k','LineWidth',2)
                    %plot([0 maxdur],[spikelines(s)+tempang spikelines(s)+tempang],'b--')
                    
                end
                
                plot([0 maxdur],[spacing(1)+tempang spacing(1)+tempang],'k')
                plot([0 maxdur],[spacing(end)+tempang spacing(end)+tempang],'k')
                
                                
            end
            
            % Plotting stuff
            xlim([0 maxdur])
            ylim([min(rots)/pi*180-mindiffrots/2 max(rots)/pi*180+mindiffrots/2])
            if r == ceil(numel(allrhos)/2)
                ylabel('Color Direction (deg)')
            end
            set(gca,'YTick',rots(2:2:end)/pi*180)
            set(gca,'xlim',[plat{p}.GaussianSpikingProfile.mu - plat{p}.GaussianSpikingProfile.sigma*2 maxdur]);

        end
        
        % Plotting stuff
        plotTitle = [filename ' Plat # ' num2str(p)];
        set(gcf,'Name',plotTitle,'NumberTitle','off')
        xlabel('Time (s) from Stimulus Onset')

        
    end
    
    
end


% Create histogram for highest contrast stimuli (the same as above)
%
% for p = 1:numel(plat)
%
%     % Set up variables
%     bins = linspace(0,max(plat{p}.trial.stimdur),50);
%     PSTH = zeros(size(bins));
%
%     for i = 1:numel(rots)
%
%         L = plat{p}.par.poltheta_orig == rots(r);
%
%         rhos = plat{p}.par.polrho_orig(L);
%
%         LL = plat{p}.par.poltheta_orig == rots(r) & plat{p}.par.polrho_orig == max(rhos);
%
%         PSTH = PSTH + histc(plat{p}.par.catntspikes{LL}',bins);
%
%     end
%
%     %Gaussian profile
%     gauss = @(x,params)(params(2)+params(1)*exp((-(x-params(3)).^2)/params(4).^2));
%     bguess = plat{p}.GaussianSpikingProfile.b;
%     aguess = plat{p}.GaussianSpikingProfile.a;
%     muguess = plat{p}.GaussianSpikingProfile.mu;
%     sigmaguess = plat{p}.GaussianSpikingProfile.sigma;
%     fittedparams = fminsearch(@(params)sum((gauss(bins,params)-PSTH).^2),[aguess,bguess,muguess,sigmaguess]);
%     offset = [fittedparams(3)-fittedparams(4) fittedparams(3)+fittedparams(4)];
%
%
%     % Plot spiking profile
%     subplot(3,1,1); hold on; grid on;
%     set(gcf,'Name',filename,'NumberTitle','off')
%     xlim([0 maxdur])
%     bar(bins,PSTH);
%     plot(bins,gauss(bins,fittedparams),'b-.','LineWidth',1)
%     plot(offset,[fittedparams(1) fittedparams(1)],'k-','linewidth',1);
%     title('Histogram of Spikes')
%     xlabel('Time (ms) from Stimulus Onset')
%     ylabel('# of Spikes')
%
%
% end
%
