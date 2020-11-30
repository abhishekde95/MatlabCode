% This code is intended as an analysis for GLMP datafiles. It compares
% within group (stimulus) variance to between group (stimulus) variance.
% Here, a group is a contrast along some ray emenating from the center, and
% it is compared with other contrasts along that ray.  The idea is that the
% the firing rate must clearly modulate in at least one direction, but not 
% necessarily all directions. This analysis find a metric by which we can 
% reject cells as too noisy or unmodulated, and keep cells that are fit for 
% analysis.

% TO DO: anova for whole group, or just along color directions?  
            %One-way vs two-way?

% 11/16/12      Created.    JPW

function plat = GLMP_AOV(plat)

disp('Performing a 1-way anova on color directions...')

filename = plat{1}.datafile;
fig = 700;

for p = 1:numel(plat)
    
    disp(['Evaluating Platform # ',num2str(p),' of ',num2str(numel(plat)),'.'])
    
    %anova takes a matrix, the groups are organized into columns
    allangs = unique(plat{p}.par.poltheta_orig);
    
    pvalues = nan(numel(allangs),1);
    
    for n = 1:numel(allangs)
        
        disp(['Checking new color direction. Current color direction = ',num2str(allangs(n)/pi*180),' degrees.'])
        
        % Index into stim along the same color direction
        radL = plat{p}.par.poltheta_orig==allangs(n);
        radIdx = find(radL);
        
        % Find maximum number of observations in a color direction
        nobservations = nan(1,numel(radIdx));
        for i = 1:numel(radIdx)
            nobservations(i) = numel(plat{p}.par.nspikes{radIdx(i)});
        end
        
        % Fill in missing nspikes
        maxobs = max(nobservations);
        observations = nan(numel(radIdx),maxobs);
        for i = 1:numel(radIdx)
            observations(i,1:nobservations(i)) = plat{p}.par.nspikes{radIdx(i)};
        end
        
        % 1-way anova
        pvalues(n) = anova1(observations',[],'off');

    end
    
    % Ssame analysis as above, but for the entire stimulus set
    nobservations = nan(1,numel(plat{p}.par.frs));
    for i = 1:numel(plat{p}.par.frs)
        nobservations(i) = numel(plat{p}.par.nspikes{i});
    end
    
    % Fill in missing nspikes
    maxobs = max(nobservations);
    observations = nan(numel(plat{p}.par.frs),maxobs);
    for i = 1:numel(plat{p}.par.frs)
        observations(i,1:nobservations(i)) = plat{p}.par.nspikes{i};
    end
    
    %1-way anova
    pval = anova1(observations',[],'off');
    

    % Return results in plat structure
    plat{p}.anova1.colordirs = pvalues;
    plat{p}.anova1.alldata = pval;
    
    if pval <= .01
        plat{p}.Analyze = 1;
    else
        plat{p}.Analyze = 0;
    end
    
    
end