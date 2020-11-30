%For ring size, choose to use ring radius (degrees) or percent of center radius (%)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%radmethod = 1; % Ring radius (degrees)
%radmethod = 2; % Percent of center radius (%)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%Choose method of calculating radius that gave peak response
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%peakmethod = 1; % Radius that gave the peak response
%peakmethod = 2; % X-axis value for peak of spline fit
%peakmethod = 3; % X-axis value for peak of Gaussian fit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%Load unit or LFP data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%data = 'Units_FixInMove_M2';  % Fixation point in ring, moving, Monkey 2
%data = 'Units_FixInStay_M2';  % Fixation point in ring, stationary, Monkey 2
%data = 'Units_FixOutStay_M2'; % Fixation point outside of ring, stationary, Monkey 2

%data = 'LFP_FixInMove_M1';    % Fixation point in ring, moving, Monkey 1
%data = 'LFP_FixInMove_M2';    % Fixation point in ring, moving, Monkey 2
%data = 'LFP_FixInStay_M1';    % Fixation point in ring, stationary, Monkey 1
%data = 'LFP_FixInStay_M2';    % Fixation point in ring, stationary, Monkey 2
%data = 'LFP_FixOutStay_M2';   % Fixation point outside of ring, stationary, Monkey 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                   


function [radii, peakGnear, peakGfar, peakCnear, peakCfar, ...
          popradii, popGnear, popGfar, popCnear, popCfar] = SMurray_folder_analysis(radmethod, peakmethod, data)

load(num2str(data));
n = size(Gznear,1);

%If radmethod = 2, change ring size to percent of center radius (%)
if radmethod == 2 %percent of center radius (%)
    for i = 1:size(radii,1)
        centerrad = radii(i,4);
        radii(i,:) = radii(i,:) ./ centerrad .* 100;
    end
end

%ring size corresponding to peak response
peakGnear = NaN(size(Gznear,1),1);
peakGfar = NaN(size(Gznear,1),1);
peakCnear = NaN(size(Cznear,1),1);
peakCfar = NaN(size(Cznear,1),1);

for i = 1:n
    
    
    %Plot gray background data
    figure
    subplot(1, 2, 1)
    hold on
    errorbar(radii(i,:), Gznear(i,:), GznearSEM(i,:), 'm', 'linewidth', 1) %Gray background, Near
    errorbar(radii(i,:), Gzfar(i,:),  GzfarSEM(i,:),  'k', 'linewidth', 1) %Gray background, Far
    legend('Gray: Near', 'Gray: Far', 'Location', 'SouthEast')
    if radmethod == 1
        xlabel('Annulus radius (deg)')
    elseif radmethod == 2
        xlabel('Annulus size (% of middle ring radius)')
    end
    ylabel('Z-scored activity')
    [xval, fitval, peak, peakval] = fitPeak(peakmethod, radii(i,:), Gznear(i,:)); %calculate radius that gave peak response
    peakGnear(i,1) = peak;
    plot(peak,peakval,'m*','MarkerSize',10)
    plot(xval, fitval, '-m', 'linewidth', 2)
    [xval, fitval, peak, peakval] = fitPeak(peakmethod, radii(i,:), Gzfar(i,:)); %calculate radius that gave peak response
    peakGfar(i,1) = peak;
    plot(peak,peakval,'k*','MarkerSize',10)
    plot(xval, fitval, '-k', 'linewidth', 2)
    hold off
    
    
    %Plot corridor background data
    subplot(1, 2, 2)
    hold on
    errorbar(radii(i,:), Cznear(i,:), CznearSEM(i,:), 'r', 'linewidth', 1) %Corridor background, Near
    errorbar(radii(i,:), Czfar(i,:),  CzfarSEM(i,:),  'b', 'linewidth', 1) %Corridor background, Far
    legend('Corr: Near', 'Corr: Far', 'Location', 'SouthEast')
    if radmethod == 1
        xlabel('Annulus radius (deg)')
    elseif radmethod == 2
        xlabel('Annulus size (% of middle ring radius)')
    end
    ylabel('Z-scored activity')
    [xval, fitval, peak, peakval] = fitPeak(peakmethod, radii(i,:), Cznear(i,:)); %calculate radius that gave peak response
    peakCnear(i,1) = peak;
    plot(peak,peakval,'r*','MarkerSize',10)
    plot(xval, fitval, '-r', 'linewidth', 2)
    [xval, fitval, peak, peakval] = fitPeak(peakmethod, radii(i,:), Czfar(i,:)); %calculate radius that gave peak response
    peakCfar(i,1) = peak;
    plot(peak,peakval,'b*','MarkerSize',10)
    plot(xval, fitval, '-b', 'linewidth', 2)
    hold off
    
end
close(gcf-n+1:gcf)


















%plot peaks of size tuning curves
figure; hold on;
axis([min(min(radii)) max(max(radii)) min(min(radii)) max(max(radii))])
axis square
plot(peakGfar, peakGnear, 'ok')
plot(peakCfar, peakCnear, 'or')
xlabel('Far: peak of size tuning curve (deg)')
ylabel('Near: peak of size tuning curve (deg)')
legend('Gray', 'Corridor', 'Location', 'SouthEast')
title(['Size tuning curve peaks: ' num2str(data)]);
plot([min(min(radii)) max(max(radii))], [min(min(radii)) max(max(radii))], '-k')
hold off


%plot function fit of population
figure; subplot(1, 2, 1); hold on;
popradii = reshape(radii, 1, size(radii,1) * size(radii,2)); %all ring sizes, population
popGnear = reshape(Gznear, 1, size(Gznear,1) * size(Gznear,2));
[Gnearxval, Gnearfitval, ~, ~] = fitPeak(peakmethod, popradii, popGnear); %calculate radius that gave peak response
if peakmethod == 1
    plot(Gnearxval, Gnearfitval, '.m', 'linewidth', 2)
else
    plot(Gnearxval, Gnearfitval, '-m', 'linewidth', 2)
end
popGfar = reshape(Gzfar, 1, size(Gzfar,1) * size(Gzfar,2));
[Gfarxval, Gfarfitval, ~, ~] = fitPeak(peakmethod, popradii, popGfar); %calculate radius that gave peak response
if peakmethod == 1
    plot(Gfarxval, Gfarfitval, '.k', 'linewidth', 2)
else
    plot(Gfarxval, Gfarfitval, '-k', 'linewidth', 2)
end
legend('Near', 'Far', 'Location', 'SouthEast')
if radmethod == 1
    xlabel('Annulus radius (deg)')
elseif radmethod == 2
    xlabel('Annulus size (% of middle ring radius)')
end
ylabel('Z-scored activity')
title('Gray background')
hold off
subplot(1, 2, 2); hold on;
popCnear = reshape(Cznear, 1, size(Cznear,1) * size(Cznear,2));
[Cnearxval, Cnearfitval, ~, ~] = fitPeak(peakmethod, popradii, popCnear); %calculate radius that gave peak response
if peakmethod == 1
    plot(Cnearxval, Cnearfitval, '.r', 'linewidth', 2)
else
    plot(Cnearxval, Cnearfitval, '-r', 'linewidth', 2)
end
popCfar = reshape(Czfar, 1, size(Czfar,1) * size(Czfar,2));
[Cfarxval, Cfarfitval, ~, ~] = fitPeak(peakmethod, popradii, popCfar); %calculate radius that gave peak response
if peakmethod == 1
    plot(Cfarxval, Cfarfitval, '.b', 'linewidth', 2)
else
    plot(Cfarxval, Cfarfitval, '-b', 'linewidth', 2)
end
legend('Near', 'Far', 'Location', 'SouthEast')
if radmethod == 1
    xlabel('Annulus radius (deg)')
elseif radmethod == 2
    xlabel('Annulus size (% of middle ring radius)')
end
ylabel('Z-scored activity')
title('Corridor background')
set(gcf, 'Name', ['Population size tuning curve: ' num2str(data)])
hold off


end