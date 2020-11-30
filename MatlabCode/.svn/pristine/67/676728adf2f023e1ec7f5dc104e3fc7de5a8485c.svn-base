% This code is intended as an analysis for GLMP datafiles. It takes a plat 
% structure, looks at all datapoints collected in a single datafile, finds 
% the best 1D projection of the data (least squares), and returns the plat 
% structure with the best fitting spline and the angle (in degrees) of that 
% projection.


% 10/29/12      Created.    JPW

function [plat] = OneDSplineAnalysis(plat)

% Projecting onto unit vectors in the LM plane and fittig a 1-D spline

filename = plat{1}.datafile;

fig = 100;
nthetas = 360;
thetas = linspace(0,pi-pi/nthetas,nthetas);
smoothingparam = .995;
nplat = numel(plat);

for p = 1:nplat;

    % Clear data for each plat
    data = nan*ones(nthetas,1);

    % Plot Projection as Data is Rotated
    figure(fig); clf; hold on; fig=fig+1;
    plotTitle = [filename ' Plat # ' num2str(p)];
    set(gcf,'Name',plotTitle,'NumberTitle','off')
    subplot(1,2,1); hold on;
    title('Spline fit to Projections');
    ylabel('Firing Rate (sp/s)');
    xlabel('Normalized Contrast');
    
    
    for i = 1:length(thetas);
        
        unitvect = [cos(thetas(i)) sin(thetas(i))]';
        projs = plat{p}.par.stim_norm * unitvect;
        lastwarn('');
        pp  = csaps(projs,plat{p}.par.meanfr,smoothingparam);
        data(i) = sum((fnval(pp,projs) - plat{p}.par.meanfr).^2);
        if ~isempty(lastwarn)
            data(i) = nan;
        end
        
        % Plot spline fit to rotating projections
        subplot(1,2,1); grid on; hold on;
        axis([-1 1 0 max(plat{p}.par.meanfr)*1.1])
        fnplt(pp);
        plot(projs,plat{p}.par.meanfr,'k.'); 
        drawnow; cla
        
        % Plot "goodness of fits"
        keyboard
        subplot(1,2,2); axis([0 180 0 max(data)*1.1]) ; hold on; grid on;
        plot(linspace(0,180-180/nthetas,nthetas),data);
        xlabel('Theta (Deg)');
        ylabel('SSE from Spline Fit');
        title('Best Projection Vector');
        
    end
    
    bp = find(data == min(data));
    preftheta = thetas(bp);
    unitvect = [cos(thetas(bp)) sin(thetas(bp))]';
    prefproj = plat{p}.par.stim_norm(:,1:2) * unitvect;
    prefpp = csaps(prefproj,plat{p}.par.meanfr,smoothingparam);
    plat{p}.par.spline.oneDprefpp = prefpp;
    plat{p}.par.spline.OneDTheta = preftheta/pi*180;
    
    % Plot best projection
    subplot(1,2,1); hold on; grid on;
    fnplt(prefpp);
    plot(prefproj,plat{p}.par.meanfr,'k.')
    
    % Plot "goodness of fits"
    subplot(1,2,2); axis([0 180 0 max(data)*1.1]) ; hold on; grid on;
    plot(linspace(0,180-180/nthetas,nthetas),data);
    xlabel('Theta (Deg)');
    ylabel('SSE from Spline Fit');
    title('Best Projection Angle');
    plot(preftheta/pi*180,min(data),'*r')

end


