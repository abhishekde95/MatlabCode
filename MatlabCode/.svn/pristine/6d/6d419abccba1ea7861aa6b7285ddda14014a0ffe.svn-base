%Part of iterateandplotfiles_modularplusdb. Tried to plot the multiple
%ellipses over one another - doesn't currently work well as far as I know
%but I'm not sure what's wrong with it. When we decided to not turn this
%into a figure for the paper, I stopped working on it.
function binnedDataPlot(Data, multiplot, dots, ellipse, fig_title)
disp('binned raw data plots');
%plotting the data into 6 bins, L = x axis, M = Y axis.
%OR plotting data into single bin L = x axis, M = Y axis
% ORIGINAL bins = logspace(log10(min(TF)), log10(max(TF)),nbins+1);
nbins = 4; figure; hold on;
bins = logspace(log10(1),log10(25), nbins+1);
if iscell(Data)
    dataSize = length(Data);
else
    dataSize = size(Data,2)/7;
end
for s = 1:dataSize
    if iscell(Data)
        dataSection = Data{s};
    else
        dataSection = Data;
    end
    uniqueXY = unique(dataSection(:,[5,6]), 'rows');
    for u = 1:size(uniqueXY,1)
        data_subset =  dataSection(dataSection(:,5) == uniqueXY(u,1) & dataSection(:,6) == uniqueXY(u,2), :);
        TF = data_subset(:, 3); Loog = logical(data_subset(:,4));
        if iscell(Data) && length(Data) > 5
            if s > 4
                color = 'k'; %humans
            else
                color = 'b'; %monkeys
            end
        else
            colorChoice = {'b','r', 'g', 'k', 'm'};
            try
                color = colorChoice{s};
            catch
                keyboard;
            end
        end
        titleString = sprintf('%s %s:%s', fig_title, uniqueXY(u,1), uniqueXY(u,2));
        set(gcf,'name',titleString,'numbertitle','off');
        actualPlotter(bins, TF, data_subset, Loog, color, nbins, multiplot, dots, ellipse)
    end
end
end

function actualPlotter(bins, TF, data, Loog, color, nbins, multiplot, dots, ellipse)
plotcounter = 1; 
if ~multiplot
    bins = bins(1:2);
end
for i = 1:length(bins)-1
    TFbounds = [bins(i) bins(i+1)];
    Dbounds = TF >= TFbounds(1) & TF <= TFbounds(2);
    miniLoog = data(Dbounds&~Loog,4); L = data(Dbounds&~Loog,1); M = data(Dbounds&~Loog, 2); file = data(Dbounds&~Loog,7);
    if ~isempty(L) && ~isempty(M)
        if multiplot
            subplot(ceil(sqrt(nbins)),ceil(sqrt(nbins)),plotcounter); hold on; %multiplot
        else
            plot(ceil(sqrt(nbins)),ceil(sqrt(nbins)));
        end
        dot_col = sprintf('%s.', color);
%         if strcmpi(color, 'r') || strcmpi(color, 'm')
%             oog_col = sprintf('%s.', 'k');
%         else
%             oog_col = sprintf('%s.', 'r');
%         end
        if dots
            plot(L,M, dot_col); plot(-L,-M,dot_col);
            plot(data(Dbounds&Loog,1),data(Dbounds&Loog,2),'r.');
            plot(-data(Dbounds&Loog,1),-data(Dbounds&Loog,2),'r.');
        end
        fitErrFn = @(params)ellipsefiterr(params,[L M],miniLoog);
        opt = optimoptions(@fmincon, 'Display', 'off', 'Algorithm','active-set');
        pts = linspace(0, 2*pi, 100); [v,d] = eig(cov([[L;-L] [M;-M]]));
    else
        plotcounter = plotcounter + 1;
        continue;
    end
    if (numel(d) == 4) % need d to be a 2x2 matrix otherwise next line crashes
        alpha_b = fmincon(fitErrFn,[sqrt(d(2,2)), sqrt(d(1,1)),atan2(v(2,2), v(1,2))],[],[],[],[],[.01 .01 0],[10 10 2*pi], [], opt);
    else
        plotcounter = plotcounter + 1;
        continue;
    end
    theta = pts+alpha_b(3); rho = (alpha_b(1)*alpha_b(2))./sqrt((alpha_b(2)*cos(pts)).^2+((alpha_b(1)*sin(pts)).^2));
    if ellipse
        polar(theta, rho, color);
    end
    title(['TF: ',num2str(geomean(TFbounds))]);
    axis square;
    if (sum(Dbounds) > 0)
        lim = max(max(abs(data(Dbounds,[1 2]))));
    else
        lim = max(abs([x;y]));
    end
    if isempty(lim)
        lim = max(abs([x;y]));
    end
    set(gca,'Xlim',lim*[-2 2],'Ylim',lim*[-2 2]);
    drawnow;
    plotcounter = plotcounter + 1;
end
end