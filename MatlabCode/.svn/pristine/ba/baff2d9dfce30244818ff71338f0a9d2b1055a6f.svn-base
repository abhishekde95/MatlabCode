% This code is intended as an analysis for GLMP datafiles. It takes a plat
% structure and fits a two dimensional surface to all the datapoints.  The
% two dimensions are orthogonal, and are systematically rotated to find the
% two vectors that fit the data the best.  For details on the fitting
% proceedure and function, see below.

% NOTE: This analysis works best with radial datasets which have equal
% L+M/L-M contrasts in the original space

% TO DO: nothing at the moment.

% 10/29/12      Created.    JPW
% 10/30/12      Modified. (Rotations now smarter, plotting better)   JPW
% 11/14/12      Modified. (Rotating through smaller changes in angle, even
%                   angles that do not have CRF's along them.)  JPW
% 11/15/12      Modified. (Fitting surfaces to the original space of the
%                   data instead of the normalized space.)  JPW
% 11/16/12      Modified. (Rotating 180 degrees just to be sure the inital
%                   guess is accurate.  Beginning with L+M, L-M, the guess 
%                   for each rotation is the fit from the previous 
%                   rotation.)   JPW
% 11/20/12      Created New Function.  (Modified TwoDFit_AllData to fit a
%                   1-D Naka Rushton.)  JPW


function plat = OneDFit_AllData(plat)

disp('Fitting A 1-D Function to All Data...')

filename = plat{1}.datafile;
fig = 200;
figure(fig); clf;

% Generating an initial guess
for p = 1:numel(plat)
    
    disp(['Evaluating Platform # ',num2str(p),' of ',num2str(numel(plat)),'.'])
    
    % Organize rotations
    angs = linspace(0,pi,181);
    nrots = numel(angs);
    maxx = [];
    
    % Preallocate some space
    if p == 1
        GOF = cell(nrots,numel(plat));
    end
    params = nan(nrots,7);
    
    % Some Variables
    vlb = [0     0   0.0001 0.0001  0  0  0];
    vub = [1000 1000  100    100   10 10 50];
    options = optimset('MaxFunEvals',5000,'MaxIter',5000,'TolFun',10^-6,'TolX',10^-6,'Display','off','Algorithm','active-set');
    sigmaguess = 0.1;
    dt = plat{p}.trial.nspikes ./ plat{p}.trial.fr;
    dt = nanmean(dt);
    baselineguess = mean(plat{p}.trial.baselinefr) * dt;
    topfr = max(plat{p}.trial.nspikes);
    if (topfr == 0)
        params0 = [0 0 1 1 2 2 0];
    elseif isempty(topfr)
        params0 = [0 0 1 1 2 2 0];
    else
        params0 = [topfr, topfr, sigmaguess, sigmaguess, 2, 2, baselineguess];  % need to constrain B across color directions
    end
    
    % Rotate through angles
    for rot = 1:nrots
        
        rotMat = [cos(angs(rot)) -sin(angs(rot)); sin(angs(rot)) cos(angs(rot))];
        tempRotPts = [rotMat *plat{p}.trial.stim_orig']';
        projs = tempRotPts(:,1);
        
        if rot>1
            %params1 = params(rot-1,:);
            params1 = mean([params(rot-1,:);params0]);
        else
            params1 = params0;
        end
        
        [f1,fval] = fmincon('FitNakaRushtonFunJPW',params1,[],[],[],[],vlb,vub,[],options,projs,plat{p}.trial.nspikes,'asymmetric');

        params(rot,:) = f1;
        GOF{rot,p} = -fval;
        
        maxx = cat(1,maxx,max(projs));
        currentGOF = cat(1,GOF{:,p});
        pts = min(projs):.001:max(projs);
        fitPts = ComputeNakaRushtonJPW(params(rot,:),pts,'asymmetric');
        
        
        %Variables for plotting
        
        %Identify Important Lines and Points in the Original Space
        axlim = max(maxx);
        Lline = [-axlim 0; axlim 0];
        Mline = [0 -axlim; 0 axlim];
        LpMline = [axlim axlim; -axlim -axlim];
        LmMline = [-axlim axlim; axlim -axlim];
        LccL = softEq(plat{p}.trial.poltheta_orig,0)...
            | softEq(plat{p}.trial.poltheta_orig,pi);
        MccL = softEq(plat{p}.trial.poltheta_orig,pi/2)...
            | softEq(plat{p}.trial.poltheta_orig,-pi/2);
        LpML = softEq(plat{p}.trial.poltheta_orig,pi/4)...
            | softEq(plat{p}.trial.poltheta_orig,-3*pi/4);
        LmML = softEq(plat{p}.trial.poltheta_orig,3*pi/4)...
            | softEq(plat{p}.trial.poltheta_orig,-pi/4);
        LlineNorm = [plat{p}.transMat*Lline']';
        MlineNorm = [plat{p}.transMat*Mline']';
        LpMlineNorm = [plat{p}.transMat*LpMline']';
        LmMlineNorm = [plat{p}.transMat*LmMline']';
        normalXAxis = [0 -axlim; 0 axlim];
        normalYAxis = [-axlim 0; axlim 0];
        
        
        %Identify Important Lines and Points in the Transformed Space
        %TLline = [rotMat*Lline']';
        %TMline = [rotMat*Mline']';
        %TLpMline = [rotMat*LpMline']';
        %TLmMline = [rotMat*LmMline']';
        [x,y] = meshgrid(linspace(-axlim,axlim,50));

%        surface = ComputeNakaRushtonJPW(params(rot,:),[x(:) y(:)],'surface1');
        %surface = ComputeNakaRushtonJPW([params(rot,1:6) zeros(1,6) params(rot,7)],[x(:) y(:)],'surface2');
       

        L = softEq(plat{p}.trial.poltheta_orig,angs(rot))...
            | softEq(plat{p}.trial.poltheta_orig,angs(rot)+pi)...
            | softEq(plat{p}.trial.poltheta_orig,angs(rot)-pi/2)...
            | softEq(plat{p}.trial.poltheta_orig,angs(rot)-pi);
        currentXAxis = [plat{p}.transMat*rotMat*normalXAxis']';
        currentYAxis = [plat{p}.transMat*rotMat*normalYAxis']';
        
        
        % Plot Surface (sanity checking)
        figure(fig); if rot==1 clf; end
        plotTitle = [filename ' Plat # ' num2str(p)];
        set(gcf,'Name',plotTitle,'NumberTitle','off')

        % Plot Porjections
        subplot(4,4,[1 2 5 6 9 10]); cla; hold on; grid on;
        title('Projections of Orignal Dataset (Rotated)')
        xlim([-max(maxx) max(maxx)])
        plot(projs,plat{p}.trial.nspikes,'k*');
        plot(pts,fitPts,'r--')
        drawnow;

        
        %Plot Normalized Space
        subplot(4,4,[3 4 7 8 11 12]); cla; hold on; grid on;
        xlabel('Lcc (Normalized)')
        ylabel('Mcc (Normalized)')
        title('Noramalized Dataset')
        axis([-1 1 -1 1])
        plot3(plat{p}.trial.stim_norm(LccL,1),plat{p}.trial.stim_norm(LccL,2),plat{p}.trial.nspikes(LccL),'r*')
        plot3(plat{p}.trial.stim_norm(MccL,1),plat{p}.trial.stim_norm(MccL,2),plat{p}.trial.nspikes(MccL),'g*')
        plot3(plat{p}.trial.stim_norm(LpML,1),plat{p}.trial.stim_norm(LpML,2),plat{p}.trial.nspikes(LpML),'m*')
        plot3(plat{p}.trial.stim_norm(LmML,1),plat{p}.trial.stim_norm(LmML,2),plat{p}.trial.nspikes(LmML),'c*')
        plot3(LlineNorm(:,1),LlineNorm(:,2),[0;0],'r--')
        plot3(MlineNorm(:,1),MlineNorm(:,2),[0;0],'g--')
        plot3(LpMlineNorm(:,1),LpMlineNorm(:,2),[0;0],'m--')
        plot3(LmMlineNorm(:,1),LmMlineNorm(:,2),[0;0],'c--')
        %legend('Lcc (Normalized)','Mcc (Normalized)','L+M (Normalized)','L-M (Normalized)','Location','SouthEastOutside')
        plot3(plat{p}.trial.stim_norm(:,1),plat{p}.trial.stim_norm(:,2),plat{p}.trial.nspikes,'ko')
        plot3(plat{p}.trial.stim_norm(L,1),plat{p}.trial.stim_norm(L,2),plat{p}.trial.nspikes(L),'*')
        plot3(currentXAxis(:,1),currentXAxis(:,2),[0;0],'k')
        plot3(currentYAxis(:,1),currentYAxis(:,2),[0;0],'k')
        
        % Plot Goodness of Fit
        subplot(4,4,[13 14 16]); hold on; grid on;
        xlim([min(angs)/pi*180 max(angs)/pi*180])
        plot(angs(1:rot)/pi*180,currentGOF(1:rot),'o-')
        
        % Sanity check surface
%         [x,y] = meshgrid(linspace(-axlim,axlim,50));
%         surface = ComputeNakaRushtonJPW([params(rot,:) zeros(1,6)],[x(:) y(:)],'surface2');
%         figure(1); cla; hold on; grid on;
%         surf(x,y,reshape(surface,size(x)))
%         plot3(tempRotPts(:,1),tempRotPts(:,2),plat{p}.trial.nspikes,'ko')
        
    end
    
    % Find best fit surface

    platGOF = cell2mat(GOF(:,p));
    bestIdx = find(platGOF==max(platGOF));
    bestL = softEq(plat{p}.trial.poltheta_orig,angs(bestIdx))...
        | softEq(plat{p}.trial.poltheta_orig,angs(bestIdx)+pi/2)...
        | softEq(plat{p}.trial.poltheta_orig,angs(bestIdx)+pi)...
        | softEq(plat{p}.trial.poltheta_orig,angs(bestIdx)-pi/2)...
        | softEq(plat{p}.trial.poltheta_orig,angs(bestIdx)-pi)... 
        | softEq(plat{p}.trial.poltheta_orig,angs(bestIdx)-3*pi/4);
    rotMat = [cos(angs(bestIdx)) -sin(angs(bestIdx)); sin(angs(bestIdx)) cos(angs(bestIdx))];
    tempRotPts = [rotMat*plat{p}.trial.stim_orig']';
    %tempRotPts = plat{p}.trial.stim_orig*rotMat';
    axlim = max(plat{p}.trial.polrho_orig);
    [x,y] = meshgrid(linspace(-axlim,axlim,50));
    surfparams = [params(bestIdx,1:6) zeros(1,6) params(bestIdx,7)];
    surface = ComputeNakaRushtonJPW(surfparams,[x(:) y(:)],'surface2');

    
    %Identify Important Lines and Points in the Original Space
    Lline = [-axlim 0; axlim 0];
    Mline = [0 -axlim; 0 axlim];
    LpMline = [axlim axlim; -axlim -axlim];
    LmMline = [-axlim axlim; axlim -axlim];
    LccL = softEq(plat{p}.trial.poltheta_orig,0)...
        | softEq(plat{p}.trial.poltheta_orig,pi);
    MccL = softEq(plat{p}.trial.poltheta_orig,pi/2)...
        | softEq(plat{p}.trial.poltheta_orig,-pi/2);
    LpML = softEq(plat{p}.trial.poltheta_orig,pi/4)...
        | softEq(plat{p}.trial.poltheta_orig,-3*pi/4);
    LmML = softEq(plat{p}.trial.poltheta_orig,3*pi/4)...
        | softEq(plat{p}.trial.poltheta_orig,-pi/4);
    LlineNorm = [plat{p}.transMat*Lline']';
    MlineNorm = [plat{p}.transMat*Mline']';
    LpMlineNorm = [plat{p}.transMat*LpMline']';
    LmMlineNorm = [plat{p}.transMat*LmMline']';
    normalXAxis = [0 -axlim; 0 axlim];
    normalYAxis = [-axlim 0; axlim 0];
    
    
    %Identify Important Lines and Points in the Transformed Space
    TLline = [rotMat*Lline']';
    TMline = [rotMat*Mline']'; 
    TLpMline = [rotMat*LpMline']';
    TLmMline = [rotMat*LmMline']';
    
    %Identify Important Lines and Points in the Normalized Space
    currentXAxis = [plat{p}.transMat*rotMat*normalXAxis']';
    currentYAxis = [plat{p}.transMat*rotMat*normalYAxis']';
    
    % Plot Goodness of Fit
    figure(fig); clf; 
    plotTitle = [filename ' Plat # ' num2str(p)];
    set(gcf,'Name',plotTitle,'NumberTitle','off')
    subplot(4,2,[7 8]); hold on; grid on;
    plot(angs/pi*180,platGOF,'o-')
    plot(angs(bestIdx)/pi*180,platGOF(bestIdx),'*r')
    title('1-D Fit Using All Data')
    xlabel('Rotation (Degrees)')
    ylabel('Log Likelihood')

    % Plot Original Surface
    subplot(4,2,[1 3 5]); cla; hold on; grid on;
    plot3(TLline(:,1),TLline(:,2),[0;0],'r--')
    plot3(TMline(:,1),TMline(:,2),[0;0],'g--')
    plot3(TLpMline(:,1),TLpMline(:,2),[0;0],'m--')
    plot3(TLmMline(:,1),TLmMline(:,2),[0;0],'c--')
    %legend('Lcc', 'Mcc', 'L+M', 'L-M')
    surf(x,y,reshape(surface,size(x)))
    axis([-axlim axlim -axlim axlim])
    plot3(tempRotPts(:,1),tempRotPts(:,2),plat{p}.trial.nspikes,'ko')
    plot3(tempRotPts(bestL,1),tempRotPts(bestL,2),plat{p}.trial.nspikes(bestL),'k*')
    plot3(tempRotPts(LccL,1),tempRotPts(LccL,2),plat{p}.trial.nspikes(LccL),'r*')
    plot3(tempRotPts(MccL,1),tempRotPts(MccL,2),plat{p}.trial.nspikes(MccL),'g*')
    plot3(tempRotPts(LpML,1),tempRotPts(LpML,2),plat{p}.trial.nspikes(LpML),'m*')
    plot3(tempRotPts(LmML,1),tempRotPts(LmML,2),plat{p}.trial.nspikes(LmML),'c*')
    zlabel('Number of Spikes')
    title('Original Dataset (Rotated)')
    

    % Plot Normalized Surface
    figure(fig)
    subplot(4,2,[2,4,6]);cla; hold on; grid on;
    plot3(plat{p}.trial.stim_norm(LccL,1),plat{p}.trial.stim_norm(LccL,2),plat{p}.trial.nspikes(LccL),'r*')
    plot3(plat{p}.trial.stim_norm(MccL,1),plat{p}.trial.stim_norm(MccL,2),plat{p}.trial.nspikes(MccL),'g*')
    plot3(plat{p}.trial.stim_norm(LpML,1),plat{p}.trial.stim_norm(LpML,2),plat{p}.trial.nspikes(LpML),'m*')
    plot3(plat{p}.trial.stim_norm(LmML,1),plat{p}.trial.stim_norm(LmML,2),plat{p}.trial.nspikes(LmML),'c*')
    plot3(LlineNorm(:,1),LlineNorm(:,2),[0;0],'r--')
    plot3(MlineNorm(:,1),MlineNorm(:,2),[0;0],'g--')
    plot3(LpMlineNorm(:,1),LpMlineNorm(:,2),[0;0],'m--')
    plot3(LmMlineNorm(:,1),LmMlineNorm(:,2),[0;0],'c--')
    %legend('Lcc (Normalized)','Mcc (Normalized)','L+M (Normalized)','L-M (Normalized)','Location','SouthEastOutside')
    plot3(plat{p}.trial.stim_norm(:,1),plat{p}.trial.stim_norm(:,2),plat{p}.trial.nspikes,'ko')
    plot3(plat{p}.trial.stim_norm(bestL,1),plat{p}.trial.stim_norm(bestL,2),plat{p}.trial.nspikes(bestL),'*')
    plot3(currentXAxis(:,1),currentXAxis(:,2),[0;0],'k')
    plot3(currentYAxis(:,1),currentYAxis(:,2),[0;0],'k')
    axis([-1 1 -1 1])
    zlabel('Number of Spikes')
    title('Normalized Dataset')
    fig=fig+1;

    keyboard
    
    % For grant application
    rotMat = [cos(angs(bestIdx)) -sin(angs(bestIdx)); sin(angs(bestIdx)) cos(angs(bestIdx))];
    tempRotPts = [(rotMat) *[plat{2}.trial.Lcc_orig'; plat{2}.trial.Mcc_orig']]';
    projs = tempRotPts(:,1);
    fitPts = ComputeNakaRushtonJPW(params(bestIdx,:),pts,'asymmetric');
    figure; hold on;
    title('Projections Onto Best Axis')
    xlim([-max(maxx) max(maxx)])
    plot(projs,plat{2}.trial.nspikes,'k*');
    plot(pts,fitPts,'r--','LineWidth',3)
    
    figure; hold on;
    plot(angs(find(angs==0):end)/pi*180,platGOF(find(angs==0):end),'o-')
    plot(angs(bestIdx)/pi*180,subunitGOF(bestIdx),'*r')
    title('1-D Fit Using All Data')
    xlabel('Rotation (Degrees)')
    ylabel('Log Likelihood')
    
    figure(3000); clf; hold on; grid on;
    surf(x,y,reshape(surface,size(x)))
    axis([-axlim axlim -axlim axlim])
%    plot3(tempRotPts(:,1),tempRotPts(:,2),plat{p}.trial.nspikes,'k*')
    %plot3(tempRotPts(bestL,1),tempRotPts(bestL,2),plat{p}.trial.nspikes(bestL),'k*')
   
    
    
    %Plot normalized stim rotated to "Natural Axes"?
%     rotNormPts = [rotNormMat*plat{p}.trial.stim_orig']';
%     TLline_norm = [rotNormMat*Lline']';
%     TMline_norm = [rotNormMat*Mline']';
%     TLpMline_norm = [rotNormMat*LpMline']';
%     TLmMline_norm = [rotNormMat*LmMline']';
%     plot3(TLline_norm(:,1),TLline_norm(:,2),[0;0],'r--')
%     plot3(TMline_norm(:,1),TMline_norm(:,2),[0;0],'g--')
%     plot3(TLpMline_norm(:,1),TLpMline_norm(:,2),[0;0],'m--')
%     plot3(TLmMline_norm(:,1),TLmMline_norm(:,2),[0;0],'c--')
%     plot3(rotNormPts(:,1),rotNormPts(:,2),plat{p}.trial.nspikes,'ko')
%     plot3(rotNormPts(bestL,1),rotNormPts(bestL,2),plat{p}.trial.nspikes(bestL),'k*')
%     plot3(rotNormPts(LccL,1),rotNormPts(LccL,2),plat{p}.trial.nspikes(LccL),'r*')
%     plot3(rotNormPts(MccL,1),rotNormPts(MccL,2),plat{p}.trial.nspikes(MccL),'g*')
%     plot3(rotNormPts(LpML,1),rotNormPts(LpML,2),plat{p}.trial.nspikes(LpML),'m*')
%     plot3(rotNormPts(LmML,1),rotNormPts(LmML,2),plat{p}.trial.nspikes(LmML),'c*')
%     plot3(currentXAxis(:,1),currentXAxis(:,2),[0;0],'k')
%     plot3(currentYAxis(:,1),currentYAxis(:,2),[0;0],'k')

    
    % Resturn analysis with plat data
    plat{p}.OneDFitAllData.All.RotationInDeg = angs'/pi*180;
    plat{p}.OneDFitAllData.All.LogLikelihood = platGOF;
    plat{p}.OneDFitAllData.All.Params = params;
    plat{p}.OneDFitAllData.Best.RotationInDeg = angs(bestIdx)/pi*180;
    plat{p}.OneDFitAllData.Best.LogLikelihood = platGOF(bestIdx);
    plat{p}.OneDFitAllData.Best.Params = params(bestIdx,:);

    
end


