%% This code is intended to calculate population data...

% Created   4/11/12     JPW
% Making a variant for cell FR variability...   7/26/12     JPW
% Changing code to look at likelihood fits
%    across the population (and saving as new script)  8/12     JPW
% Changing code to look at variance of
%    residuals vs predicted values (new script)     9/12    JPW


clear all
close all

library = 'Users/jpatrickweller/Documents/Matlab_Working_Directory/Patrick/GridLMPlane/GLMPDatafiles/';
n=1;

% Bring in raw data
% Cartesian Grid
%datafile{n} = [char(library) 'S110211007.nex']; n = n+1;% 1 Symmetric Trench (Most Data)
%datafile{n} = [char(library) 'S110411012_2.nex']; n = n+1;% 2
%datafile{n} = [char(library) 'S110911003.nex']; n = n+1;% 3 Symmetric Trench
%datafile{n} = [char(library) 'S110911007.nex']; n = n+1;% 4 Symmetric Trench
%datafile{n} = [char(library) 'S111011003.nex']; n = n+1;% 5 Asymmetric Trench (Second to Most Data)
%datafile{n} = [char(library) 'S112311006.nex']; n = n+1;% 6 Asymmetric Trench
%datafile{n} = [char(library) 'S113011002.nex']; n = n+1;% 7 Complex Cell
%datafile{n} = [char(library) 'S120511008.nex']; n = n+1;% 8 Symmetric Trench

% Radial Grid
% (Mostly) Luminance
datafile{n} = 'S121211003.nex'; n = n+1;% 9 Symmetric Trench
datafile{n} = 'S121211006.nex'; n = n+1;% 10 Symmetric L-M Trench
datafile{n} = 'S121511004.nex'; n = n+1;% 11 Symmetric L-M Trench (Small Dataset)

% (Mostly) Chromatic
datafile{n} = 'S041912004.nex'; n = n+1;% 12 Chromatic Cell? (Bug: Theta Unknown)

% Ellipsoidal
%datafile{n} = 'S120811003.nex'; n = n+1;% 13 Ellipsoidal Cell? Extremely
        %noisy. NOTE: NOT INCLUDING B/C FIT SURFACE HAS COMPLEX NUMBERS
datafile{n} = 'S120911003.nex'; n = n+1;% 14 Ellipsoidal Cell? Extremely noisy.
datafile{n} = 'S042012002.nex'; n = n+1;% 15 Ellipsoidal Cell! (Bug: Theta Unknown)
datafile{n} = 'A060112003.nex'; n = n+1;% 16 Ellispoidal Cell! (Apollo's First Cell!)
datafile{n} = 'A062612003.nex'; n = n+1;% 17 Ellispoidal Cell!
%datafile{n} = 'A062612004.nex'; n = n+1;% 18 Ellipsoidal (Continuation of 062612003. Not getting all waveforms?)
datafile{n} = 'A070312003.nex'; n = n+1;% 19 Ellipsoidal Cell (Tiny Dataset)
datafile{n} = 'A070312005.nex'; n = n+1;% 20 Ellipsoidal Cell
datafile{n} = 'A071012004.nex'; n = n+1;% 21 Ellipsoidal Cell (10 repeats)
datafile{n} = 'A071212005.nex'; n = n+1;% 22 Ellipsoidal Cell
datafile{n} = 'A071612002.nex'; n = n+1;% 23 Ellipsoidal Cell
datafile{n} = 'A071612005.nex'; n = n+1;% 24 Ellipsoidal Cell (Same cell as A071612002)

% Plateaus: S-Cone
datafile{n} = 'S011312002.nex'; n = n+1;% 25 Ellipse!
datafile{n} = 'S020112002.nex'; n = n+1;% 26 Symmetric L-M Trench (check iso)
datafile{n} = 'S020212007.nex'; n = n+1;% 27 Asymmetric L-M Trench

% Plateaus: 2 Directions of Motion
datafile{n} = 'S021012007.nex'; n = n+1;% 28 Symmetric L-M Trench
datafile{n} = 'S032312002.nex'; n = n+1;% 29 Luminance Cell
datafile{n} = 'S032912004.nex'; n = n+1;% 30 Symmetric L-M Trench
datafile{n} = 'S040412002.nex'; n = n+1;% 31 Chromatic Cell!!
datafile{n} = 'S041012002.nex'; n = n+1;% 32 Ellipsoidal Cell!!
%datafile{n} = 'S040612002.nex'; n = n+1;% 33 Great iso, but not much structure...

% Hybrid Paradigm
datafile{n} = 'S101812003.nex'; n = n+1;% 34 
datafile{n} = 'S102312007.nex'; n = n+1;% 35 


% Bugs...
%rawdata = 'S032812004.nex';% Dropped 2 headers - must trim orignal file
%rawdata = 'S032812006.nex';% Grating File!

% Setting up some variables
fig = 4;
totpl = 1;
sig2gauss = 8;

ndatafiles = size(datafile,2);

% Begin grand loop...
for loop = 1%1:ndatafiles
    
    %Organize Raw Data
    rawdata = nex2stro([char(library) datafile{loop}]);
    disp(['Now processing ' datafile{loop}])
    [trial par plat] = OrganizeRawGLMPData(rawdata);
    
    
    %% 2-D Fit Surfaces To All Data
    
    nplat = numel(plat);
    for p = 1:nplat
        
        disp('Fitting 2D surface to all of the data...')
        
        % Fit 2-D Surfaces to All Data (Testing only cardinal axes)
        vlb = [0    0 0.001 0.001 1 1 0];
        vub = [200 200    1     1 6 6 100];
        options = optimset('MaxFunEvals',5000,'MaxIter',5000,'TolFun',10^-6,'TolX',10^-6,'Display','off','Algorithm','active-set');
        sigmaguess = 0.1;
        maxpred = []; f = [];
        titles = {'L+M','L-M'};
        for i = 1:2
            if (i == 1)
                L = softEq(plat{p}.trial.poltheta_norm,0) | softEq(plat{p}.trial.poltheta_norm,pi);
            else
                L = softEq(plat{p}.trial.poltheta_norm,pi/2) | softEq(plat{p}.trial.poltheta_norm,-pi/2);
            end
            projs = plat{p}.trial.polrho_norm(L);
            topfr = max(plat{p}.trial.nspikes(L));
            if (topfr == 0)
                f(i,:) = [0 0 1 1 2 2 0];
            elseif isempty(topfr)
                f(i,:) = [0 0 1 1 2 2 0];
            else
                params0 = [topfr, topfr, sigmaguess, sigmaguess, 2, 2, 0];  % need to constrain B across color directions
                f(i,:) = fmincon('MyFitNakaRushtonFun',params0,[],[],[],[],vlb,vub,[],options,projs,plat{p}.trial.nspikes(L),'asymmetric');
            end
            
            % May want to use all data to generate 1-d nakarushtons, which
            % in turn generate initial guesses...
            pred = MyComputeNakaRushton(f(i,:),projs,'asymmetric');
            
        end
        a = mean(f(1,[1 2]));
        b = mean(f(2,[1 2]));
        params0 = max(plat{p}.trial.nspikes);
        params0(2) = a./sum([a b]); %this doesn't work so well
        params0(2) = .5;
        params0(3) = .25;
        params0(4) = 2;
        params0(5) = (f(1,2)./f(1,4))./(f(1,1)./f(1,3)); if (isnan(params0(5))) params0(5) = 1; end;
        params0(6) = (f(2,2)./f(2,4))./(f(2,1)./f(2,3)); if (isnan(params0(6))) params0(6) = 1; end;
        params0(7) = mean(f(:,end));
        
        vlb = [10   0 0.0001 1  0  0  0];
        vub = [1000 1     1  6 10 10 20];
        options = optimset('MaxFunEvals',5000,'MaxIter',5000,'TolFun',10^-6,'TolX',10^-6,'Display','off','Algorithm','active-set');
        [f0,fval] = fmincon('MyFitNakaRushtonFun',params0,[],[],[],[],vlb,vub,[],options,plat{p}.trial.stim_norm,plat{p}.trial.nspikes,'surface1');
        [x,y] = meshgrid(linspace(min(plat{p}.trial.LpM_norm),max(plat{p}.trial.LpM_norm),50), linspace(min(plat{p}.trial.LmM_norm),max(plat{p}.trial.LmM_norm),50));
        surface = MyComputeNakaRushton(f0,[x(:) y(:)],'surface1');
        
        %         % Plot Surface
        %         figure(fig); fig=fig+1; hold on; grid on;
        %         surf(x,y,reshape(surface,size(x)))
        %         axis([-1 1 -1 1])
        %         tempRotPts = plat{p}.trial.stim_norm*rotMat;
        %         plot3(tempRotPts(:,1),tempRotPts(:,2),plat{p}.trial.nspikes,'ko')
        %         xlabel('L-M')
        %         ylabel('L+M')
        %         zlabel('Number of Spikes')
        %         title(['Fit to all the data (n = ',num2str(size(plat{p}.trial.LpM_norm,1)),')']);
        
        
        %         % Plot Non-Rotated Surface
        %         figure(fig); fig=fig+1; hold on; grid on;
        %         surf(x,y,reshape(surface,size(x)))
        %         axis([-1 1 -1 1])
        %         plot3(tempRotPts(:,1),tempRotPts(:,2),plat{p}.trial.nspikes,'ko')
        %         xlabel('L-M')
        %         ylabel('L+M')
        %         zlabel('Number of Spikes')
        %         title('Non-Rotated Surface (Fit to All Pts)')
        
        
        % Calculate residuals
        residSurf = MyComputeNakaRushton(f0,plat{p}.par.stim_norm,'surface1');
        %resids = residSurf - plat{p}.trial.nspikes;
        resids = residSurf - plat{p}.par.meannspikes;
        resMAD(totpl) = regress(resids,residSurf);
        
        
        %         % Plot Residuals
        %         figure(fig); fig=fig+1; hold on; grid on;
        %         plot(residSurf,resids,'*')
        %         plot(residSurf,residSurf*resM(totpl),'k--')
        %         xlabel('Prediction (# of Spikes)')
        %         ylabel('Residual (Prediction - Actual Spikes)')
        %         title('Prediction vs. Residuals')
        
        % Calculate slope of absolute value of residuals
        resMabsAD(totpl) = regress(abs(resids),residSurf);
        
        % Plot absolute value of Residuals
        figure(fig); subplot(2,5,[4 5]); hold on; grid on;
        set(gcf,'Name',datafile{loop},'NumberTitle','off')
        plot(residSurf,abs(resids),'*')
        plot(residSurf,residSurf*resMabsAD(totpl),'k--')
        xlabel('Prediction (# of Spikes)')
        ylabel('Residual (Absolute Value of Mean # of Spikes)')
        title('Prediction vs. Residuals (Fit to All Data)')
        
        
        
        %% Building 2-D Surfaces with 2 Vectors
        
        disp('Using 2 orthagonal vectors to create 2D surfaces...')
        
        % Generating an initial guess
        vlb = [0    0 0.001 0.001 1 1 0];
        vub = [200 200    1     1 6 6 100];
        options = optimset('MaxFunEvals',50000,'MaxIter',50000,'TolFun',10^-6,'TolX',10^-6,'Display','off','Algorithm','active-set');
        sigmaguess = 0.1;
        maxpred = []; f = [];
        
        for i = 1:2
            if (i == 1)
                L = softEq(plat{p}.trial.poltheta_norm,0) | softEq(plat{p}.trial.poltheta_norm,pi);
            else
                L = softEq(plat{p}.trial.poltheta_norm,pi/2) | softEq(plat{p}.trial.poltheta_norm,-pi/2);
            end
            projs = plat{p}.trial.polrho_norm(L);
            topfr = max(plat{p}.trial.nspikes(L));
            if (topfr == 0)
                f(i,:) = [0 0 1 1 2 2 0];
            else
                params1 = [topfr, topfr, sigmaguess, sigmaguess, 2, 2, 0];  % need to constrain B across color directions
                f(i,:) = fmincon('MyFitNakaRushtonFun',params1,[],[],[],[],vlb,vub,[],options,projs,plat{p}.trial.nspikes(L),'asymmetric');
            end
            
            % May want to use all data to generate 1-d nakarushtons, which
            % in turn generate initial guesses...
            pred = MyComputeNakaRushton(f(i,:),projs,'asymmetric');
            
        end
        
        % Redefine the index as all 4 directions
        L = softEq(plat{p}.trial.poltheta_norm,-pi/2) | softEq(plat{p}.trial.poltheta_norm,0)...
            | softEq(plat{p}.trial.poltheta_norm,pi/2) | softEq(plat{p}.trial.poltheta_norm,pi);
        
        a = mean(f(1,[1 2]));
        b = mean(f(2,[1 2]));
        params1 = max(plat{p}.trial.nspikes);
        params1(2) = .5;
        params1(3) = .25;
        params1(4) = 2;
        params1(5) = (f(1,2)./f(1,4))./(f(1,1)./f(1,3)); if (isnan(params1(5))) params1(5) = 1; end;
        params1(6) = (f(2,2)./f(2,4))./(f(2,1)./f(2,3)); if (isnan(params1(6))) params1(6) = 1; end;
        params1(7) = mean(f(:,end));
        
        vlb = [10   0 0.0001 1  0  0  0];
        vub = [1000 1     1  6 10 10 20];
        options = optimset('MaxFunEvals',50000,'MaxIter',50000,'TolFun',10^-6,'TolX',10^-6,'Display','off','Algorithm','active-set');
        [f1] = fmincon('MyFitNakaRushtonFun',params1,[],[],[],[],vlb,vub,[],options,plat{p}.trial.stim_norm(L,:),plat{p}.trial.nspikes(L),'surface1');
        [x,y] = meshgrid(linspace(min(plat{p}.par.LpM_norm),max(plat{p}.par.LpM_norm),50), linspace(min(plat{p}.par.LmM_norm),max(plat{p}.par.LmM_norm),50));
        surface = MyComputeNakaRushton(f1,[x(:) y(:)],'surface1');

        %         % Plot Surface
        %         figure(fig); fig=fig+1; hold on; grid on;
        %         surf(x,y,reshape(surface,size(x)))
        %         axis([-1 1 -1 1])
        %         plot3(rotPts(L,1),rotPts(L,2),plat{p}.trial.nspikes(L),'k*')
        %         plot3(rotPts(~L,1),rotPts(~L,2),plat{p}.trial.nspikes(~L),'go');
        %         xlabel('L-M')
        %         ylabel('L+M')
        %         zlabel('Number of Spikes')
        %         title(['Fit to Only 2 Vectors: ',num2str(angs(rot)+pi/2),' and ',num2str(angs(rot)+pi)])
        
        
        % Plot Non-Rotated Surface
        figure(fig); subplot(2,5,[1 2 6 7]); hold on; grid on;
        surf(x,y,reshape(surface,size(x)))
        axis([-1 1 -1 1])
        plot3(plat{p}.par.LmM_norm,plat{p}.par.LpM_norm,plat{p}.par.meannspikes,'ko')
        xlabel('L-M')
        ylabel('L+M')
        zlabel('Mean # of Spikes')
        title('Surface Fit to 2-Axes')
        
        
        % Calculate residuals
        residSurf = MyComputeNakaRushton(f1,plat{p}.par.stim_norm,'surface1');
        %resids = residSurf - plat{p}.trial.nspikes;
        resids = residSurf - plat{p}.par.meannspikes;
        resM2A(totpl) = regress(resids,residSurf);
        
        
        %         % Plot Residuals
        %         figure(fig); fig=fig+1; hold on; grid on;
        %         plot(residSurf,resids,'g*')
        %         plot(residSurf,residSurf*resM,'k--')
        %         xlabel('Prediction (# of Spikes)')
        %         ylabel('Residual (Prediction - Actual Spikes)')
        %         title('Prediction vs. Residuals')
        
        % Calculate slope of absolute value of residuals
        resMabs2A(totpl) = regress(abs(resids),residSurf);
        
        % Plot absolute value of Residuals
        figure(fig); subplot(2,5,[9 10]); fig=fig+1; hold on; grid on;
        plot(residSurf,abs(resids),'g*')
        plot(residSurf,residSurf*resMabs2A(totpl),'k--')
        xlabel('Prediction (# of Spikes)')
        ylabel('Residual (Absolute Value of Mean)')
        title('Prediction vs. Residuals (Fit to 2 Axes)')
        
        totpl = totpl+1;
    end
    
    
end


%% Plot final analysis

figure(fig); fig=fig+1; hold on; grid on;
[n1, xout1] = hist(resMabsAD);
[n2, xout2] = hist(resMabs2A);
bar(xout2,n2,'g')
bar(xout1,n1)
legend('Fit to All Data','Fit Using 2 Axes')
xlabel('Slope of Regression')
ylabel('# of Cells')
title('Slope of Pred vs Res')



