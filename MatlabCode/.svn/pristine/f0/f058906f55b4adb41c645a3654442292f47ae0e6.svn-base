% This code is intended as an analysis for GLMP datafiles. It takes a plat
% structure and fits the best fit two dimensional surface to two orthogonal 
% axes.  It then finds the direction of maximal sensitivity.

% NOTE: Currently, this code works best with radial datasets.

% TO DO: 
 
% 11/9/12      Created.    JPW

function plat = ConeWeightAnalysis2(plat)

disp('Finding Most Sensitive Direction')

fig = 450;

for p = 1:numel(plat)
         
    disp(['Evaluating Platform # ',num2str(p),' of ',num2str(numel(plat)),'.'])
    
    if ~isfield(plat{p},'TwoDFitAllData')
        
        disp('Must First Fit Suface to Data. Executing TwoDFit_AllData...')
        
        plat = TwoDFit_AllData(plat);
    
    end   
        
    bestAng = plat{p}.TwoDFitAllData.Best.RotationInDeg/180*pi;
    rotMat = [cos(bestAng) -sin(bestAng); sin(bestAng) cos(bestAng)];
    rotStim = [rotMat*plat{p}.trial.stim_orig']';
    Lline = [-max(plat{p}.trial.polrho_orig) 0; max(plat{p}.trial.polrho_orig) 0];
    Mline = [0 -max(plat{p}.trial.polrho_orig); 0 max(plat{p}.trial.polrho_orig)];
    rotLline = [rotMat*Lline']';
    rotMline = [rotMat*Mline']';
    LL = softEq(plat{p}.trial.poltheta_orig,0)...
        | softEq(plat{p}.trial.poltheta_orig,pi);
    ML = softEq(plat{p}.trial.poltheta_orig,pi/2)...
        | softEq(plat{p}.trial.poltheta_orig,-pi/2);
    
    bestL = softEq(plat{p}.trial.poltheta_orig,bestAng)...
        | softEq(plat{p}.trial.poltheta_orig,bestAng+pi/2)...
        | softEq(plat{p}.trial.poltheta_orig,bestAng+pi)...
        | softEq(plat{p}.trial.poltheta_orig,bestAng-pi/2)...
        | softEq(plat{p}.trial.poltheta_orig,bestAng-pi);
    
    
    [x,y] = meshgrid(linspace(min(rotStim(:,1)),max(rotStim(:,1)),50), linspace(min(rotStim(:,2)),max(rotStim(:,2)),50));
    surface = MyComputeNakaRushton(plat{p}.TwoDFitAllData.Best.Params,[x(:) y(:)],'surface1');
    
    
    % Replot Results
    figure(fig); clf; subplot(1,2,1); hold on; grid on;
    plotTitle = [plat{1}.datafile ' Plat # ' num2str(p)];
    set(gcf,'Name',plotTitle,'NumberTitle','off')
    plot3(rotLline(:,1),rotLline(:,2),[0;0],'b--')
    plot3(rotMline(:,1),rotMline(:,2),[0;0],'g--')
    legend('L-Isolating Direction', 'M-Isolating Direction')
    surf(x,y,reshape(surface,size(x)))
    plot3(rotStim(:,1),rotStim(:,2),plat{p}.trial.nspikes,'ko')
    plot3(rotStim(LL,1),rotStim(LL,2),plat{p}.trial.nspikes(LL),'b*')
    plot3(rotStim(ML,1),rotStim(ML,2),plat{p}.trial.nspikes(ML),'g*')
    plot3(rotStim(bestL,1),rotStim(bestL,2),plat{p}.trial.nspikes(bestL),'k*')
    zlabel('# of Spikes')
    title('Rotated Data')
    axis([min(rotStim(:,1)) max(rotStim(:,1)) min(rotStim(:,2)) max(rotStim(:,2))])
    
    
    % Find most sensitive axis
    [cartX0,cartY0] = pol2cart(linspace(-pi,pi-pi/360,360),.01);
    responses0 = MyComputeNakaRushton(plat{p}.TwoDFitAllData.Best.Params,[cartX0' cartY0'],'surface1');
    [maxResp,maxIdx] = max(responses0);
    [theta,rho] = cart2pol(cartX0(maxIdx),cartY0(maxIdx));
    [cartX,cartY] = pol2cart(linspace(-pi,pi-pi/360,360),max(plat{p}.trial.polrho_orig));
    responses = MyComputeNakaRushton(plat{p}.TwoDFitAllData.Best.Params,[cartX' cartY'],'surface1');
    
    % Plot Results
    plot([0 cartX(maxIdx)],[0 cartY(maxIdx)]','r')
    plot3(cartX(maxIdx),cartY(maxIdx),responses(maxIdx),'r*')
    plot3([0 cartX(maxIdx)],[0 cartY(maxIdx)],[0 responses(maxIdx)],'r--')
    plot3([cartX(maxIdx) cartX(maxIdx)],[cartY(maxIdx) cartY(maxIdx)],[0 responses(maxIdx)],'r--')
    
    %Plot orignals (sanity check)
    figure(fig); fig=fig+1;
    subplot(1,2,2); hold on; grid on;
    plot3(plat{p}.trial.Lcc_orig,plat{p}.trial.Mcc_orig,plat{p}.trial.nspikes,'ko')
    temp = [inv(rotMat)*[cartX(maxIdx) cartY(maxIdx)]']';
    plot3(temp(1),temp(2),responses(maxIdx),'r*')
    plot3([0 temp(1)],[0 temp(2)],[0 responses(maxIdx)],'r--')
    plot3([0 temp(1)],[0 temp(2)],[0 0],'r')
    plot3([temp(1) temp(1)],[temp(2) temp(2)],[0 responses(maxIdx)],'r--')
    plot3(plat{p}.trial.Lcc_orig(LL),plat{p}.trial.Mcc_orig(LL),plat{p}.trial.nspikes(LL),'b*')
    plot3(plat{p}.trial.Lcc_orig(ML),plat{p}.trial.Mcc_orig(ML),plat{p}.trial.nspikes(ML),'g*')
    xlabel('L-Isolating Direction')
    ylabel('M-Isolating Direction')
    zlabel('# of Spikes')
    title('Sanity Check')
    
    mostSensitiveDir = bestAng/pi*180 + theta/pi*180;
    
    if sign(mostSensitiveDir) == -1
        mostSensitiveDir = mostSensitiveDir + 360;
    end        
    
    
    plat{p}.ConeWeightAnal2.MostSensitiveDirInDeg = mostSensitiveDir;

    
end    

