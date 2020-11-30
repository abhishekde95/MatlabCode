
% This code is intended as an analysis for GLMP datafiles. It takes a plat
% structure and compares L-isolating and M-isolating stimuli of equal
% contrast.  This should replicate the cone weights of Shapley, etc.

% NOTE: Currently, this code works best with radial datasets.

% TO DO: Maybe include L/M ratio in output...
 
% 11/01/12      Created.    JPW

function [plat] = ConeWeightAnalysis(plat)

disp('Calculating Cone Weights by Comparing Like Contrasts in Cone Isolating Directions...')

fig = 350;
filename = plat{1}.datafile;


for p = 1:numel(plat)
    
    disp(['Evaluating Pmlatform # ',num2str(p),' of ',num2str(numel(plat)),'.'])
    
    
    % Identify stim along L and M isolating directions    
    LPosIsoL = plat{p}.par.poltheta_orig == pi/4;
    LNegIsoL = plat{p}.par.poltheta_orig == -3*pi/4;
    MPosIsoL = plat{p}.par.poltheta_orig == 3*pi/4;
    MNegIsoL = plat{p}.par.poltheta_orig == -pi/4;
    
    % Identify contrast of stimuli
    LPosCC = plat{p}.par.polrho_orig(LPosIsoL);
    LNegCC = plat{p}.par.polrho_orig(LNegIsoL);
    MPosCC = plat{p}.par.polrho_orig(MPosIsoL);
    MNegCC = plat{p}.par.polrho_orig(MNegIsoL);
    
    % Identify FR of stimuli
    LPosFR = plat{p}.par.meanfr(LPosIsoL);
    LNegFR = plat{p}.par.meanfr(LNegIsoL);
    MPosFR = plat{p}.par.meanfr(MPosIsoL);
    MNegFR = plat{p}.par.meanfr(MNegIsoL);
    

    % Exclude any contrast that does not have a partner
    [C,Lidx,Midx] = setxor(LPosCC,MPosCC);
    LPosCC(Lidx) = [];
    LPosFR(Lidx) = [];
    MPosCC(Midx) = [];
    MPosFR(Midx) = [];
    [C,Lidx,Midx] = setxor(LNegCC,MNegCC);
    LNegCC(Lidx) = [];
    LNegFR(Lidx) = [];
    MNegCC(Midx) = [];
    MNegFR(Midx) = [];
    
    
    % Fix row/column agreement in empty sets (there must be a better solution...)
    if isrow(LPosCC)
        LPosCC = LPosCC';
    end
    if isrow(LNegCC)
        LNegCC = LNegCC';
    end
    if isrow(MPosCC)
        MPosCC = MPosCC';
    end
    if isrow(MNegCC)
        MNegCC = MNegCC';
    end
    if isrow(LPosFR)
        LPosFR = LPosFR';
    end
    if isrow(LNegFR)
        LNegFR = LNegFR';
    end
    if isrow(MPosFR)
        MPosFR = MPosFR';
    end
    if isrow(MNegFR)
        MNegFR = MNegFR';
    end
    

    % Sort order of stim to compare like contrasts
    [LPosCC,sortOrd] = sort(LPosCC);
    LPosFR = LPosFR(sortOrd);
    [MPosCC,sortOrd] = sort(MPosCC);
    MPosFR = MPosFR(sortOrd);
    [LNegCC,sortOrd] = sort(LNegCC);
    LNegFR = LNegFR(sortOrd);
    [MNegCC,sortOrd] = sort(MNegCC);
    MNegFR = MNegFR(sortOrd);
    
    % Calculate Cone Weights Ratio (L/M)
    PosCW = LPosFR./MPosFR;
    NegCW = LNegFR./MNegFR;
    
    % Embed analysis in plat structure for return
    % Positive Cone Weights
    plat{p}.ConeWeightAnal.Pos.LIsoCC = LPosCC;
    plat{p}.ConeWeightAnal.Pos.MIsoCC = MPosCC;
    plat{p}.ConeWeightAnal.Pos.LIsoFR = LPosFR;
    plat{p}.ConeWeightAnal.Pos.MIsoFR = MPosFR;

    % Negative Cone Weights
    plat{p}.ConeWeightAnal.Neg.LIsoCC = LNegCC;
    plat{p}.ConeWeightAnal.Neg.MIsoCC = MNegCC;
    plat{p}.ConeWeightAnal.Neg.LIsoFR = LNegFR;
    plat{p}.ConeWeightAnal.Neg.MIsoFR = MNegFR;
    
    
    % Set up variables for plotting
    maxPos = max(max(LPosFR),max(MPosFR));
    maxNeg = max(max(LNegFR),max(MNegFR));
    
    % Display Results
    figure(fig); clf; hold on; grid on; fig=fig+1;
    plotTitle = [filename ' Plat # ' num2str(p)];
    set(gcf,'Name',plotTitle,'NumberTitle','off')
    
    subplot(1,2,1); cla; hold on; grid on;
    plot(LPosFR,MPosFR,'*')
    plot(0:.1:maxPos,0:.1:maxPos,'k--')
    xlabel('Response to L-Cone Isolating Stimuli (FR)')
    ylabel('Response to M-Cone Isolating Stimuli (FR)')
    title('Positive Cone-Isoating Contrasts')
    axis equal square
    
    subplot(1,2,2); cla; hold on; grid on;
    plot(LNegFR,MNegFR,'*')
    plot(0:.1:maxNeg,0:.1:maxNeg,'k--')
    xlabel('Response to L-Cone Isolating Stimuli (FR)')
    ylabel('Response to M-Cone Isolating Stimuli (FR)')
    title('Negative Cone-Isoating Contrasts')
    axis equal square


    
end

