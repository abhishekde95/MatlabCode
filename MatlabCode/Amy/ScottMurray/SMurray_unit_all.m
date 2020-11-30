%plots of single unit data
%combines all blocks for one day's data set

clear all
[fname, pathname] = uigetfile('*.nex', 'Select a NeuroExplorer file');
day = fname([1:8 11:14]);
allfiles = dir(pathname);
row = 1;
for i = 1:size(allfiles,1)
    if size(allfiles(i).name,2) == 14
        if allfiles(i).name([1:8 11:14]) == day
            allblocks(row,:) = allfiles(i).name;
            row = row + 1;
        end
    end
end

rowG = 1;
rowC = 1;
rowB = 1;
rowP = 1;
interstimG = [];
interstimC = [];
interstimB = [];
interstimP = [];

for i = 1:size(allblocks,1)
    stro = nex2stro(strcat(pathname, allblocks(i,:)));
    aw = [0.05; 0.05]; %analysis window
    [bkgnd, near, far, uniqueors, spikename, interstim, nearPre farPre] = SMurray_unit_spikes(stro, aw);
    
    if bkgnd == 2 %gray background
        nearG(rowG:rowG+size(near,1)-1, :, :) = near;
        farG(rowG:rowG+size(far,1)-1, :, :) = far;
        nearPreG(rowG:rowG+size(nearPre,1)-1, :, :) = nearPre;
        farPreG(rowG:rowG+size(farPre,1)-1, :, :) = farPre;
        rowG = rowG + size(near,1);
        interstimG = [interstimG; interstim];
    elseif bkgnd == 0 %corridor background
        nearC(rowC:rowC+size(near,1)-1, :, :) = near;
        farC(rowC:rowC+size(far,1)-1, :, :) = far;
        nearPreC(rowC:rowC+size(nearPre,1)-1, :, :) = nearPre;
        farPreC(rowC:rowC+size(farPre,1)-1, :, :) = farPre;
        rowC = rowC + size(near,1);
        interstimC = [interstimC; interstim];
    elseif bkgnd == 1 %brick background
        nearB(rowB:rowB+size(near,1)-1, :, :) = near;
        farB(rowB:rowB+size(far,1)-1, :, :) = far;
        nearPreB(rowB:rowB+size(nearPre,1)-1, :, :) = nearPre;
        farPreB(rowB:rowB+size(farPre,1)-1, :, :) = farPre;
        rowB = rowB + size(near,1);
        interstimB = [interstimB; interstim];
    elseif bkgnd == 5 %Ponzo background
        nearP(rowP:rowP+size(near,1)-1, :, :) = near;
        farP(rowP:rowP+size(far,1)-1, :, :) = far;
        nearPreP(rowP:rowP+size(nearPre,1)-1, :, :) = nearPre;
        farPreP(rowP:rowP+size(farPre,1)-1, :, :) = farPre;
        rowP = rowP + size(near,1);
        interstimP = [interstimP; interstim];
    end
    
end





    
for k = 1:size(near,3)
    
    %raw data plot: FIRING RATE
    figure;
    %gray background
    if exist('nearG','var')
        subplot(2, 2, 1)
        hold on;
        errorbar(uniqueors/10, mean(nearG(:,:,k)), std(nearG(:,:,k))./sqrt(size(nearG,1)), '.-m', 'linewidth', 2) %degrees
        errorbar(uniqueors/10, mean(farG(:,:,k)), std(farG(:,:,k))./sqrt(size(farG,1)), '.-k', 'linewidth', 2) %degrees
        plot(uniqueors/10, ones(length(uniqueors))*mean(interstimG(:,k)), '--k', 'linewidth', 1) %baseline
        plot(uniqueors/10, mean(nearPreG(:,:,k)), '.--m', 'linewidth', 0.5) %baseline
        plot(uniqueors/10, mean(farPreG(:,:,k)), '.--k', 'linewidth', 0.5) %baseline
        xlabel('Annulus radius (deg)');
        ylabel('Firing rate (sp/sec)');
        title(['GRAY, spike data: ' num2str(spikename(k,:))]);
        hold off
    end
    %corridor background
    if exist('nearC','var')
        subplot(2, 2, 3)
        hold on;
        errorbar(uniqueors/10, mean(nearC(:,:,k)), std(nearC(:,:,k))./sqrt(size(nearC,1)), '.-m', 'linewidth', 2) %degrees
        errorbar(uniqueors/10, mean(farC(:,:,k)), std(farC(:,:,k))./sqrt(size(farC,1)), '.-k', 'linewidth', 2) %degrees
        plot(uniqueors/10, ones(length(uniqueors))*mean(interstimC(:,k)), '--k', 'linewidth', 1) %baseline
        plot(uniqueors/10, mean(nearPreC(:,:,k)), '.--m', 'linewidth', 0.5) %baseline
        plot(uniqueors/10, mean(farPreC(:,:,k)), '.--k', 'linewidth', 0.5) %baseline
        xlabel('Annulus radius (deg)');
        ylabel('Firing rate (sp/sec)');
        title(['CORRIDOR, spike data: ' num2str(spikename(k,:))]);
        hold off
    end

    
    %z-scored data plot: FIRING RATE
    %gray background
    if exist('nearG','var')
        nearGz = NaN(size(nearG));
        farGz = NaN(size(farG));
        for i = 1:size(nearG,1)
            both(1,:) = nearG(i,:,k);
            both(2,:) = farG(i,:,k);
            [zscores] = zscorefunc(both, uniqueors);
            nearGz(i,:,k) = zscores(1,:);
            farGz(i,:,k)  = zscores(2,:);
        end
        subplot(2, 2, 2)
        hold on;
        errorbar(uniqueors/10, nanmean(nearGz(:,:,k)), nanstd(nearGz(:,:,k))./sqrt(sum(~isnan(nearGz(:,:,k)))), '.-m', 'linewidth', 2) %degrees
        errorbar(uniqueors/10, nanmean(farGz(:,:,k)), nanstd(farGz(:,:,k))./sqrt(sum(~isnan(farGz(:,:,k)))), '.-k', 'linewidth', 2) %degrees
        xlabel('Annulus radius (deg)');
        ylabel('Firing rate (z-scored)');
        title(['GRAY, Z-scored spike data: ' num2str(spikename(k,:))]);
        hold off
    end
    if exist('nearC','var')
        nearCz = NaN(size(nearC));
        farCz = NaN(size(farC));
        for i = 1:size(nearC,1)
            both(1,:) = nearC(i,:,k);
            both(2,:) = farC(i,:,k);
            [zscores] = zscorefunc(both, uniqueors);
            nearCz(i,:,k) = zscores(1,:);
            farCz(i,:,k)  = zscores(2,:);
        end
        subplot(2, 2, 4)
        hold on;
        errorbar(uniqueors/10, nanmean(nearCz(:,:,k)), nanstd(nearCz(:,:,k))./sqrt(sum(~isnan(nearCz(:,:,k)))), '.-m', 'linewidth', 2) %degrees
        errorbar(uniqueors/10, nanmean(farCz(:,:,k)), nanstd(farCz(:,:,k))./sqrt(sum(~isnan(farCz(:,:,k)))), '.-k', 'linewidth', 2) %degrees
        xlabel('Annulus radius (deg)');
        ylabel('Firing rate (z-scored)');
        title(['CORRIDOR, Z-scored spike data: ' num2str(spikename(k,:))]);
        hold off
    end

    
    
    %raw data plot: RESPONSE (firing rate - spontaneous)
    figure;
    %gray background
    if exist('nearG','var')
        nearRespG = nearG - nearPreG;
        farRespG = farG - farPreG;
        subplot(2, 2, 1)
        hold on;
        errorbar(uniqueors/10, mean(nearRespG(:,:,k)), std(nearRespG(:,:,k))./sqrt(size(nearRespG,1)), '.-m', 'linewidth', 2) %degrees
        errorbar(uniqueors/10, mean(farRespG(:,:,k)), std(farRespG(:,:,k))./sqrt(size(farRespG,1)), '.-k', 'linewidth', 2) %degrees
        plot(uniqueors/10, zeros(length(uniqueors)), '--k', 'linewidth', 1) %baseline
        xlabel('Annulus radius (deg)');
        ylabel('Response (sp/sec)');
        title(['GRAY, spike data: ' num2str(spikename(k,:))]);
        hold off
    end
    %corridor background
    if exist('nearC','var')
        nearRespC = nearC - nearPreC;
        farRespC = farC - farPreC;
        subplot(2, 2, 3)
        hold on;
        errorbar(uniqueors/10, mean(nearRespC(:,:,k)), std(nearRespC(:,:,k))./sqrt(size(nearRespC,1)), '.-m', 'linewidth', 2) %degrees
        errorbar(uniqueors/10, mean(farRespC(:,:,k)), std(farRespC(:,:,k))./sqrt(size(farRespC,1)), '.-k', 'linewidth', 2) %degrees
        plot(uniqueors/10, zeros(length(uniqueors)), '--k', 'linewidth', 1) %baseline
        xlabel('Annulus radius (deg)');
        ylabel('Response (sp/sec)');
        title(['CORRIDOR, spike data: ' num2str(spikename(k,:))]);
        hold off
    end

    
    %z-scored data plot: RESPONSE (firing rate - spontaneous)
    %gray background
    if exist('nearG','var')
        nearRespGz = NaN(size(nearRespG));
        farRespGz = NaN(size(farRespG));
        for i = 1:size(nearRespG,1)
            both(1,:) = nearRespG(i,:,k);
            both(2,:) = farRespG(i,:,k);
            [zscores] = zscorefunc(both, uniqueors);
            nearRespGz(i,:,k) = zscores(1,:);
            farRespGz(i,:,k)  = zscores(2,:);
        end
        subplot(2, 2, 2)
        hold on;
        errorbar(uniqueors/10, nanmean(nearRespGz(:,:,k)), nanstd(nearRespGz(:,:,k))./sqrt(sum(~isnan(nearRespGz(:,:,k)))), '.-m', 'linewidth', 2) %degrees
        errorbar(uniqueors/10, nanmean(farRespGz(:,:,k)), nanstd(farRespGz(:,:,k))./sqrt(sum(~isnan(farRespGz(:,:,k)))), '.-k', 'linewidth', 2) %degrees
        xlabel('Annulus radius (deg)');
        ylabel('Response (z-scored)');
        title(['GRAY, Z-scored spike data: ' num2str(spikename(k,:))]);
        hold off
    end
    if exist('nearC','var')
        nearRespCz = NaN(size(nearRespC));
        farRespCz = NaN(size(farRespC));
        for i = 1:size(nearRespC,1)
            both(1,:) = nearRespC(i,:,k);
            both(2,:) = farRespC(i,:,k);
            [zscores] = zscorefunc(both, uniqueors);
            nearRespCz(i,:,k) = zscores(1,:);
            farRespCz(i,:,k)  = zscores(2,:);
        end
        subplot(2, 2, 4)
        hold on;
        errorbar(uniqueors/10, nanmean(nearRespCz(:,:,k)), nanstd(nearRespCz(:,:,k))./sqrt(sum(~isnan(nearRespCz(:,:,k)))), '.-m', 'linewidth', 2) %degrees
        errorbar(uniqueors/10, nanmean(farRespCz(:,:,k)), nanstd(farRespCz(:,:,k))./sqrt(sum(~isnan(farRespCz(:,:,k)))), '.-k', 'linewidth', 2) %degrees
        xlabel('Annulus radius (deg)');
        ylabel('Response (z-scored)');
        title(['CORRIDOR, Z-scored spike data: ' num2str(spikename(k,:))]);
        hold off
    end


end
