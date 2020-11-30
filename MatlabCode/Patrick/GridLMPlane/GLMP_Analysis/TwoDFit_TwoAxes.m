
% This code is intended as an analysis for GLMP datafiles. It takes a plat
% structure and fits a two dimensional surface to two orthogonal axes.  The
% axes are systematically rotated to find the two vectors that best predict
% the data.  For details on the fitting proceedure and function, see below.

% NOTE: Currently, this code works best with radial datasets.

% TO DO: Nothing at the moment.
 
% 10/30/12      Created.    JPW
% 10/31/12      Modified (rotations now smarter, plotting better)   JPW


function plat = TwoDFit_TwoAxes(plat)


disp('Using 2 orthagonal vectors to create 2D surfaces...')


filename = plat{1}.datafile;
fig = 300;

% Generating an initial guess
for p = 1:numel(plat)
    
    disp(['Evaluating Platform # ',num2str(p),' of ',num2str(numel(plat)),'.'])
    
    % Organize rotations
    deltaAngs = min(diff(unique(plat{p}.trial.poltheta_norm)));
    angs = 0:deltaAngs:pi/2-deltaAngs;
    nrots = numel(angs);
    
    % Preallocate some space
    if p == 1
        GOF = nan(nrots,numel(plat));
    end
    params = nan(nrots,7);
    f = nan(2,7);

    % Rotate through angles
    for rot = 1:nrots
        
        disp(['Rotating Data. Current Rotation = ',num2str(angs(rot)/pi*180),' degrees.'])
        
        % Evaluate NakaRushton along single vectors to calculate guesses for 2-D fit.
        vlb = [0    0 0.001 0.001 1 1 0];
        vub = [200 200    1     1 6 6 100];
        options = optimset('MaxFunEvals',500,'MaxIter',500,'TolFun',10^-6,'TolX',10^-6,'Display','off','Algorithm','active-set');
        sigmaguess = 0.1;
        for i = 1:2
            if (i == 1)
                L = plat{p}.trial.poltheta_norm == angs(rot) ...
                    | softEq(plat{p}.trial.poltheta_norm,angs(rot)+pi)...
                    | softEq(plat{p}.trial.poltheta_norm,angs(rot)-pi);
            else
                L = softEq(plat{p}.trial.poltheta_norm,angs(rot)+pi/2)...
                    | softEq(plat{p}.trial.poltheta_norm,angs(rot)+3*pi/2)...
                    | softEq(plat{p}.trial.poltheta_norm,angs(rot)-pi/2)...
                    | softEq(plat{p}.trial.poltheta_norm,angs(rot)-3*pi/2);
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
            
        end
        
        L = plat{p}.trial.poltheta_norm == angs(rot)...
            | softEq(plat{p}.trial.poltheta_norm,angs(rot)+pi)...
            | softEq(plat{p}.trial.poltheta_norm,angs(rot)+pi/2)...
            | softEq(plat{p}.trial.poltheta_norm,angs(rot)-pi/2)...
            | softEq(plat{p}.trial.poltheta_norm,angs(rot)-pi);
        
        params0 = max(plat{p}.trial.nspikes);
        params0(2) = .5;
        params0(3) = .25;
        params0(4) = 2;
        params0(5) = (f(1,2)./f(1,4))./(f(1,1)./f(1,3)); if (isnan(params0(5))) params0(5) = 1; end;
        params0(6) = (f(2,2)./f(2,4))./(f(2,1)./f(2,3)); if (isnan(params0(6))) params0(6) = 1; end;
        params0(7) = mean(f(:,end));
        
        vlb = [10   0 0.0001 1  0  0  0];
        vub = [1000 1     1  6 10 10 20];
        options = optimset('MaxFunEvals',5000,'MaxIter',5000,'TolFun',10^-6,'TolX',10^-6,'Display','off','Algorithm','active-set');
        rotMat = [cos(angs(rot)) -sin(angs(rot)); sin(angs(rot)) cos(angs(rot))];
        tempRotPts = plat{p}.trial.stim_norm*rotMat;
        [f0] = fmincon('MyFitNakaRushtonFun',params0,[],[],[],[],vlb,vub,[],options,tempRotPts(L,:),plat{p}.trial.nspikes(L),'surface1');
        params(rot,:) = f0;
        
        GOF(rot,p) = -MyFitNakaRushtonFun(f0,tempRotPts,plat{p}.trial.nspikes,'surface1');
        
    end
    
    % Find best fit surface
    platGOF = GOF(:,p);
    bestIdx = find(platGOF==max(platGOF));
    rotMat = [cos(angs(bestIdx)) -sin(angs(bestIdx)); sin(angs(bestIdx)) cos(angs(bestIdx))];
    LmMline = [-1 0; 1 0];
    LpMline = [0 -1; 0 1];
    TLmMline = LmMline*rotMat;
    TLpMline = LpMline*rotMat;
    [x,y] = meshgrid(linspace(min(plat{p}.trial.LpM_norm),max(plat{p}.trial.LpM_norm),50), linspace(min(plat{p}.trial.LmM_norm),max(plat{p}.trial.LmM_norm),50));
    surface = MyComputeNakaRushton(params(bestIdx,:),[x(:) y(:)],'surface1');
    L = softEq(plat{p}.trial.poltheta_norm,angs(bestIdx))...
        | softEq(plat{p}.trial.poltheta_norm,angs(bestIdx)+pi)...
        | softEq(plat{p}.trial.poltheta_norm,angs(bestIdx)+pi/2)...
        | softEq(plat{p}.trial.poltheta_norm,angs(bestIdx)-pi/2)...
        | softEq(plat{p}.trial.poltheta_norm,angs(bestIdx)-pi); 
    LmML = softEq(plat{p}.trial.poltheta_norm,0)...
        | softEq(plat{p}.trial.poltheta_norm,pi);
    LpML = softEq(plat{p}.trial.poltheta_norm,pi/2)...
        | softEq(plat{p}.trial.poltheta_norm,-pi/2);
    tempRotPts = plat{p}.trial.stim_norm*rotMat;
    
    % Plot Log Likelihood
    figure(fig); clf; fig=fig+1;
    plotTitle = [filename ' Plat # ' num2str(p)];
    set(gcf,'Name',plotTitle,'NumberTitle','off')
    subplot(4,1,4); hold on; grid on;   
    plot(angs/pi*180,platGOF,'o-')
    plot(angs(bestIdx)/pi*180,platGOF(bestIdx),'*r')
    title('Fit Using Two Axes')
    xlabel('Rotation (Degrees)')
    ylabel('Log Likelihood')

    % Plot Surface
    subplot(4,1,1:3); hold on; grid on;
    plot3(TLmMline(:,1),TLmMline(:,2),[0;0],'b--')
    plot3(TLpMline(:,1),TLpMline(:,2),[0;0],'g--')
    legend('L - M (Normalized)', 'L + M (Normalized)')
    surf(x,y,reshape(surface,size(x)))
    axis([-1 1 -1 1])
    plot3(tempRotPts(L,1),tempRotPts(L,2),plat{p}.trial.nspikes(L),'k*')
    plot3(tempRotPts(:,1),tempRotPts(:,2),plat{p}.trial.nspikes,'ko')
    plot3(tempRotPts(LmML,1),tempRotPts(LmML,2),plat{p}.trial.nspikes(LmML),'b*')
    plot3(tempRotPts(LpML,1),tempRotPts(LpML,2),plat{p}.trial.nspikes(LpML),'g*')
    zlabel('Number of Spikes')
    title(['Fit to Two Best Axes (axes = ',num2str(angs(bestIdx)/pi*180),' and ',num2str((angs(bestIdx)+pi)/pi*180),')']);
    
        
    % Resturn analysis with plat data
    rotData = cat(1,cellstr('Rotation (deg)'),num2cell(angs/pi*180)');
    GOFData = cat(1,cellstr('Log Likelihood'),num2cell(platGOF));
    plat{p}.TwoDFitTwoAxes = cat(2,rotData,GOFData);
    
end


