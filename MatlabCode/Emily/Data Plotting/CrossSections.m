%This was another attempt to make a figure showing multiple overlapping
%ellipses and subject data. We ended up not using this either, but it was
%previously a part of iterateandplotfiles. It was overwritten by
%binneddataplot, which also ended up not being used.
function CrossSections(Data, Besttheta, Fpar, SplitTitle)
disp('Cross Sections plots');
figure; nbins = 9;
for s = 1:length(Data)
    plotcounter = 1;
    data = Data{s}; besttheta = Besttheta{s}; fpar = Fpar{s};
    f1 = @(omega)fpar(1)*abs(((1i*2*pi*fpar(5).*omega+1).^-fpar(3))-fpar(2)*((1i*2*pi*fpar(6).*omega+1).^-fpar(4)));
    f2 = @(omega)fpar(1+6)*abs(((1i*2*pi*fpar(5+6).*omega+1).^-fpar(3+6))-fpar(2+6)*((1i*2*pi*fpar(6+6).*omega+1).^-fpar(4+6)));
    TF = data(:, 3); Loog = logical(data(:,4));
    bins = logspace(log10(min(TF)), log10(max(TF)),nbins+1);
    color = ['k', 'r', 'b']; 
    name = SplitTitle{s}{1}; whoseData = strsplit(name, 'L');
    lineName = [whoseData{1}]; dataName = [whoseData{1} '''s data'];
    for i = 1:length(bins)-1
        %if s > 4
        %    color = 'm';
        %else
        %    color = 'b';
        %end
        TFbounds = [bins(i) bins(i+1)];
        Dbounds = TF >= TFbounds(1) & TF <= TFbounds(2);
        subplot(ceil(sqrt(nbins)),ceil(sqrt(nbins)),plotcounter); hold on;
        %plot(data(Dbounds&~Loog,1),data(Dbounds&~Loog,2),'k.', 'DisplayName', dataName);
        %plot(-data(Dbounds&~Loog,1),-data(Dbounds&~Loog,2),'k.','DisplayName', dataName);
        %plot(data(Dbounds&Loog,1),data(Dbounds&Loog,2),'r.', 'DisplayName', dataName);
        %plot(-data(Dbounds&Loog,1),-data(Dbounds&Loog,2),'r.','DisplayName', dataName);
        y1 = f1(mean(TFbounds)).^-1; %lum thresholds
        y2 = f2(mean(TFbounds)).^-1;
        thtmp = linspace(0,2*pi,100)';
        rtmp = (y1.*y2)./sqrt((y2.*cos(thtmp)).^2+(y1.*sin(thtmp)).^2); % radius of ellipse - thank you, Wikipedia
        [x,y] = pol2cart(thtmp+besttheta,rtmp);
        plot(x,y,'color', color(s), 'DisplayName', lineName); %numerical color call? 
        hold on;
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
        plotcounter = plotcounter + 1;
    end
    hold on;
end
legend(gca, 'show', 'Location', 'bestoutside');