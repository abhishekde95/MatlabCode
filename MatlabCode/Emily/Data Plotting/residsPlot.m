%part of iterate and plot files
function hfig=residsPlot(data, modelParams, splitTitle, hax)
% residsPlot(data, modelParams, [splitTitle],[axes handle])
%
% Plots the residuals from a fit of the 2D Watson temporal contrast
% threshold surface.
%
% INPUT
%      data: nx4 matrix [L,M,TF,OOG] 
%      modelparams: 13-element vector of model parameters
%      splitTitle: Title put at the top of the figure window
%      axes handle: optional argument for an axes handle where the plot
%      should go.
% OUTPUT
%      h: handle to the figure

    if (nargin < 4)
        hax = [];
    end
    if (nargin < 3)
        splitTitle = [];
    end
    %disp('resids plot function');
    L = data(:, 1);
    M = data(:, 2);
    TF = data(:,3);
    Loog = data(:,4);
    % Figure 1
    % Predicted vs actual thresholds scatter plot
    predr = LMTF_thresh_from_model(L, M,  TF, modelParams);
    [th,r] = cart2pol(L,M);
    try
    % Looking at resids from fit
    %predlr = log10(predr(~Loog));

    % 3D plot of residuals
    if (isempty(hax))
        hfig = figure; axes; hold on;
    else
        axes(hax);
    end
    [tmpl,tmpm] = pol2cart(th,log10(r)-log10(predr));
    OS = r>predr; % Overestimating sensitivity
    if sum(~Loog & OS)>0
        h(1) = plot3(tmpl(~Loog & OS),tmpm(~Loog & OS),data(~Loog & OS,3), 'b.','MarkerSize',10);
        plot3(-tmpl(~Loog & OS),-tmpm(~Loog & OS),data(~Loog & OS,3), 'b.','MarkerSize',10);    
    end
    if sum(Loog & OS)>0 % Only plotting LOOG residuals if the model overestimates sensitivity.
        plot3(tmpl(Loog & OS),tmpm(Loog & OS),data(Loog & OS,3), 'bo','MarkerSize',5);
        plot3(-tmpl(Loog & OS),-tmpm(Loog & OS),data(Loog & OS,3), 'bo','MarkerSize',5);    
    end
    US = ~Loog & r<predr; % Underestimating sensitivity
    if sum(US)>0
        h(2) = plot3(tmpl(US),tmpm(US),data(US,3), 'r.','MarkerSize',10);
        plot3(-tmpl(US),-tmpm(US),data(US,3), 'r.','MarkerSize',10);
    end
    plot3([0, 0], [0, 0], [min(TF), max(TF)],'LineWidth',2);
    set(gca, 'zscale', 'log','Zlim',[min(TF) max(TF)],'View',[0 90]);
    set(gca,'xlim',max(abs([tmpl;tmpm]))*[-1 1],'ylim',max(abs([tmpl;tmpm]))*[-1 1]);
    if (isempty(hax))
        labels = {'Overest. sens.','Underest. sens.'};
        legend(h(isgraphics(h)),labels(isgraphics(h)),'Location','NorthEastOutside');
    else
        axis normal;  % Needed by LMTFBrowser
    end
    xlabel('log10(L)'); ylabel('log10(M)'); zlabel('TF (Hz)');
    if ~isempty(splitTitle)
        set(gcf, 'name', splitTitle{1}, 'numbertitle', 'off');
    end
    catch
        keyboard
    end
end