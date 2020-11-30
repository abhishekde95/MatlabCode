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
% 11/21/12      Created.  (Fitting to surface2 instead of surface1.)  JPW


function plat = TwoDFit_AllData2(plat,usenorm)

disp('Fitting A 2-D Function to All Data...')

filename = plat{1}.datafile;
fig = 200;

if nargin < 2
    usenorm = 0;
end

% Generating an initial guess
for p = 1:numel(plat)
    
    if isfield(plat{p},'Analyze')
        
        if plat{p}.Analyze == 1
            
            disp(['Evaluating Platform # ',num2str(p),' of ',num2str(numel(plat)),'.'])
            
            % Organize rotations
            angs = -pi/4:pi/128:7*pi/4;
            nrots = numel(angs);
            
            % Preallocate some space
            if p == 1
                GOF = cell(nrots,numel(plat));
            end
            params = nan(nrots,13);
            f = nan(2,7);
            
            % Rotate through angles
            for rot = 1:nrots
                
                disp(['Rotating Data. Current Rotation = ',num2str(angs(rot)/pi*180),' degrees.'])
                
                % Evaluate NakaRushton along single vectors to calculate guesses for 2-D fit.
                
                if rot == 1
                    vub = [1000 1000     6  6 2 2 50];
                    vlb = [0   0 0.0001 .0001 2 2  0];
                    options = optimset('MaxFunEvals',5000,'MaxIter',5000,'TolFun',10^-6,'TolX',10^-6,'Display','off','Algorithm','active-set');
                    sigmaguess = 0.01;
                    dt = plat{p}.trial.nspikes ./ plat{p}.trial.fr;
                    dt = nanmean(dt);
                    baselineguess = mean(plat{p}.trial.baselinefr) * dt;
                    for i = 1:2
                        if (i == 1)
                            L = softEq(plat{p}.trial.poltheta_orig,angs(rot)) ...
                                | softEq(plat{p}.trial.poltheta_orig,angs(rot)+pi)...
                                | softEq(plat{p}.trial.poltheta_orig,angs(rot)-pi);
                        else
                            L = softEq(plat{p}.trial.poltheta_orig,angs(rot)+pi/2)...
                                | softEq(plat{p}.trial.poltheta_orig,angs(rot)+3*pi/2)...
                                | softEq(plat{p}.trial.poltheta_orig,angs(rot)-pi/2)...
                                | softEq(plat{p}.trial.poltheta_orig,angs(rot)-3*pi/2);
                        end
                        if usenorm == 1
                            projs = plat{p}.trial.polrho_norm(L);
                        else
                            projs = plat{p}.trial.polrho_orig(L);
                        end
                        topfr = max(plat{p}.trial.nspikes(L));
                        if (topfr == 0)
                            f(i,:) = [0 0 1 1 2 2 0];
                        elseif isempty(topfr)
                            f(i,:) = [0 0 1 1 2 2 0];
                        else
                            params0 = [topfr topfr sigmaguess sigmaguess 2 2 baselineguess];  % need to constrain B across color directions
                            f(i,:) = fmincon('FitNakaRushtonFunJPW',params0,[],[],[],[],vlb,vub,[],options,projs,plat{p}.trial.nspikes(L),'asymmetric',1);
                        end
                    end
                    
                    
                    params0 = [f(1,1:6) f(2,1:6) mean(f(:,7))];
                    params01 = params0;
                    
                    
                    if usenorm == 1
                        axlim = max(plat{p}.trial.polrho_norm);
                    else
                        axlim = max(plat{p}.trial.polrho_orig);
                    end
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
                    
                else
                    
                    params0 = params(rot-1,:);
                    
                    if ~isreal(cell2mat(GOF(rot-1)))
                        disp('Warning: Imaginary Values in the Fit...')
                        
                        GOF{rot-1,:} = 0;
                        params0 = params01;
                        
                        
                        %                 reals = nan(1,rot-1);
                        %                 for n=1:rot-1
                        %                     reals(n)=isreal(cell2mat(GOF(n)));
                        %                 end
                        %                 lastreal = find(reals,1,'last');
                        %                 params0 = params(lastreal,:);
                        %
                    end
                    
                end
                
                %Fit Surface to Rotated Points
                rotMat = [cos(angs(rot)) -sin(angs(rot)); sin(angs(rot)) cos(angs(rot))];
                
                if usenorm == 1
                    tempRotPts = [rotMat*plat{p}.trial.stim_norm']';
                else
                    tempRotPts = [rotMat*plat{p}.trial.stim_orig']';
                end
                
                vub = [1000 1000 100 100 2 2 1000 1000 100 100 2 2 50];
                vlb = [0 0 .00001 .00001 2 2 0 0 .00001 .00001 2 2 0];
                [f0,fval] = fmincon('FitNakaRushtonFunJPW',params0,[],[],[],[],vlb,vub,[],options,tempRotPts,plat{p}.trial.nspikes,'surface2');
                if isnan(f0(1)) %This is a hack fix...
                    f0 = params0;
                end
                params(rot,:) = f0;
                GOF{rot,p} = -fval;
                currentGOF = cat(1,GOF{:,p});
                
                %Variables for plotting
                TLline = [rotMat*Lline']';
                TMline = [rotMat*Mline']';
                TLpMline = [rotMat*LpMline']';
                TLmMline = [rotMat*LmMline']';
                [x,y] = meshgrid(linspace(-axlim,axlim,50));
                surface = ComputeNakaRushtonJPW(params(rot,:),[x(:) y(:)],'surface2');
                surface(sign(surface) == -1) = 0;
                L = softEq(plat{p}.trial.poltheta_orig,angs(rot))...
                    | softEq(plat{p}.trial.poltheta_orig,angs(rot)+pi/2)...
                    | softEq(plat{p}.trial.poltheta_orig,angs(rot)+pi)...
                    | softEq(plat{p}.trial.poltheta_orig,angs(rot)+3*pi/2)...
                    | softEq(plat{p}.trial.poltheta_orig,angs(rot)-pi/2)...
                    | softEq(plat{p}.trial.poltheta_orig,angs(rot)-pi)...
                    | softEq(plat{p}.trial.poltheta_orig,angs(rot)+3*pi/2);
                currentXAxis = [plat{p}.transMat*rotMat*normalXAxis']';
                currentYAxis = [plat{p}.transMat*rotMat*normalYAxis']';
                
                % Plot Surface (sanity checking)
                figure(fig); if rot==1 clf; end
                plotTitle = [filename ' Plat # ' num2str(p)];
                set(gcf,'Name',plotTitle,'NumberTitle','off')
                
                % Plot Orignal (Rotated) Datapoints
                subplot(4,4,[1 2 5 6 9 10]); cla; hold on; grid on;
                title('Dataset Being Fit')
                if usenorm == 1
                    axlim = max(plat{p}.trial.polrho_norm);
                else
                    axlim = max(plat{p}.trial.polrho_orig);
                end
                axis([-axlim axlim -axlim axlim])
                plot3(tempRotPts(LccL,1),tempRotPts(LccL,2),plat{p}.trial.nspikes(LccL),'r*')
                plot3(tempRotPts(MccL,1),tempRotPts(MccL,2),plat{p}.trial.nspikes(MccL),'g*')
                plot3(tempRotPts(LpML,1),tempRotPts(LpML,2),plat{p}.trial.nspikes(LpML),'m*')
                plot3(tempRotPts(LmML,1),tempRotPts(LmML,2),plat{p}.trial.nspikes(LmML),'c*')
                plot3(TLline(:,1),TLline(:,2),[0;0],'r--')
                plot3(TMline(:,1),TMline(:,2),[0;0],'g--')
                plot3(TLpMline(:,1),TLpMline(:,2),[0;0],'m--')
                plot3(TLmMline(:,1),TLmMline(:,2),[0;0],'c--')
                %legend('Lcc', 'Mcc', 'L+M', 'L-M')
                surf(x,y,reshape(surface,size(x)))
                plot3(tempRotPts(:,1),tempRotPts(:,2),plat{p}.trial.nspikes,'ko')
                plot3(tempRotPts(L,1),tempRotPts(L,2),plat{p}.trial.nspikes(L),'k*')
                plot3(tempRotPts(L,1),tempRotPts(L,2),plat{p}.trial.nspikes(L),'b*')
                
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
                title('2-D Fit Using All Data')
                xlabel('Rotation (Degrees)')
                ylabel('Log Likelihood')
                
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
            
            if usenorm == 1
                tempRotPts = [rotMat*plat{p}.trial.stim_norm']';
            else
                tempRotPts = [rotMat*plat{p}.trial.stim_orig']';
            end
            
            surface = ComputeNakaRushtonJPW(params(bestIdx,:),[x(:) y(:)],'surface2');
            
            
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
            title('2-D Fit Using All Data')
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
            title('Fit Dataset')
            
            
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
            
            
            % Resturn analysis with plat data
            plat{p}.TwoDFitAllData.All.RotationInDeg = angs'/pi*180;
            plat{p}.TwoDFitAllData.All.LogLikelihood = platGOF;
            plat{p}.TwoDFitAllData.All.Params = params;
            plat{p}.TwoDFitAllData.Best.RotationInDeg = angs(bestIdx)/pi*180;
            plat{p}.TwoDFitAllData.Best.LogLikelihood = platGOF(bestIdx);
            plat{p}.TwoDFitAllData.Best.Params = params(bestIdx,:);
        end
        
    end
    
end


