%% This code is intended to parameterize variability in cell responses...

function [plat] = GLMP_Response_Variability(plat)

% A new analysis script for GridLMPlane

% Created   4/11/12     JPW
% Making a variant for cell FR variability...   7/26/12     JPW
% Modified (Made script into function)  1/14/12     JPW


for p = 1:numel(plat)
    
    L = plat{p}.par.nsamps ~= 1;
    [b,bint,r,rint,stats] = regress(plat{p}.par.varnspikes(L),plat{p}.par.meannspikes(L));
    x = 0:.1:max(plat{p}.par.meannspikes(L))*1.1;
    y = b * x;
    
    plat{p}.ResponseVariability.LinearRegression.Mean2Var = b;
    
    figure(404); if p == 1; clf; end
    filename = plat{p}.datafile;
    plotTitle = [filename ' Plat # ' num2str(p)];
    set(gcf,'Name',plotTitle,'NumberTitle','off')
    subplot(numel(plat),1,p); hold on; grid on;
    plot(x,y,'--')
    legend(['Slope = ' num2str(b)])
    plot(plat{p}.par.meannspikes(L),plat{p}.par.varnspikes(L),'*k')
    title(['Plateau # ' num2str(p)])
    xlabel('Mean of Spike Count')
    ylabel('Variance in Spike Count')
    
end


