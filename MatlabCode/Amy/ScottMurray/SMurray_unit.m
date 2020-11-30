%plots of single unit data

clear all
stro = nex2stro
aw = [0.05; 0.05]; %analysis window
[bkgnd, near, far, uniqueors, spikename, interstim, nearPre farPre] = SMurray_unit_spikes(stro, aw);

    
for k = 1:size(near,3)
    
    %raw data plot: FIRING RATE
    figure;
    subplot(2, 2, 1)
    hold on;
    errorbar(uniqueors/10, mean(near(:,:,k)), std(near(:,:,k))./sqrt(size(near,1)), '.-m', 'linewidth', 2) %degrees
    errorbar(uniqueors/10, mean(far(:,:,k)), std(far(:,:,k))./sqrt(size(far,1)), '.-k', 'linewidth', 2) %degrees
    plot(uniqueors/10, ones(length(uniqueors))*mean(interstim(:,k)), '--k', 'linewidth', 1) %baseline
    plot(uniqueors/10, mean(nearPre(:,:,k)), '.--m', 'linewidth', 0.5) %baseline
    plot(uniqueors/10, mean(farPre(:,:,k)), '.--k', 'linewidth', 0.5) %baseline
    xlabel('Annulus radius (deg)');
    ylabel('Firing rate (sp/sec)');
    title(['Raw spike data: ' num2str(spikename(k,:))]);
    hold off
    
    %z-scored data plot: FIRING RATE
    nearz = NaN(size(near));
    farz = NaN(size(far));
    for i = 1:size(near,1)
        both(1,:) = near(i,:,k);
        both(2,:) = far(i,:,k);
        [zscores] = zscorefunc(both, uniqueors);
        nearz(i,:,k) = zscores(1,:);
        farz(i,:,k)  = zscores(2,:);
    end
    subplot(2, 2, 2)
    hold on;
    errorbar(uniqueors/10, nanmean(nearz(:,:,k)), nanstd(nearz(:,:,k))./sqrt(sum(~isnan(nearz(:,:,k)))), '.-m', 'linewidth', 2) %degrees
    errorbar(uniqueors/10, nanmean(farz(:,:,k)), nanstd(farz(:,:,k))./sqrt(sum(~isnan(farz(:,:,k)))), '.-k', 'linewidth', 2) %degrees
    xlabel('Annulus radius (deg)');
    ylabel('Firing rate (z-scored)');
    title(['Z-scored spike data: ' num2str(spikename(k,:))]);
    hold off
    
    %raw data plot: RESPONSE (firing rate - spontaneous)
    nearResp = near - nearPre;
    farResp = far - farPre;
    subplot(2, 2, 3)
    hold on;
    errorbar(uniqueors/10, mean(nearResp(:,:,k)), std(nearResp(:,:,k))./sqrt(size(nearResp,1)), '.-m', 'linewidth', 2) %degrees
    errorbar(uniqueors/10, mean(farResp(:,:,k)), std(farResp(:,:,k))./sqrt(size(farResp,1)), '.-k', 'linewidth', 2) %degrees
    plot(uniqueors/10, zeros(length(uniqueors)), '--k', 'linewidth', 1) %baseline
    xlabel('Annulus radius (deg)');
    ylabel('Response (sp/sec)');
    title(['Raw spike data: ' num2str(spikename(k,:))]);
    hold off
    
    %z-scored data plot: RESPONSE (firing rate - spontaneous)
    nearRespz = NaN(size(nearResp));
    farRespz = NaN(size(farResp));
    for i = 1:size(nearResp,1)
        both(1,:) = nearResp(i,:,k);
        both(2,:) = farResp(i,:,k);
        [zscores] = zscorefunc(both, uniqueors);
        nearRespz(i,:,k) = zscores(1,:);
        farRespz(i,:,k)  = zscores(2,:);
    end
    subplot(2, 2, 4)
    hold on;
    errorbar(uniqueors/10, nanmean(nearRespz(:,:,k)), nanstd(nearRespz(:,:,k))./sqrt(sum(~isnan(nearRespz(:,:,k)))), '.-m', 'linewidth', 2) %degrees
    errorbar(uniqueors/10, nanmean(farRespz(:,:,k)), nanstd(farRespz(:,:,k))./sqrt(sum(~isnan(farRespz(:,:,k)))), '.-k', 'linewidth', 2) %degrees
    xlabel('Annulus radius (deg)');
    ylabel('Response (z-scored)');
    title(['Z-scored spike data: ' num2str(spikename(k,:))]);
    hold off
    
end
