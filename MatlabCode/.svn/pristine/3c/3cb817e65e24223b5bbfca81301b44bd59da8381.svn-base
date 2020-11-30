%plots of single unit data
%combines all blocks for one day's data set
%plots each day of data in a folder

clear all
[fname, pathname] = uigetfile('*.nex', 'Select a NeuroExplorer file');
folder = pathname;
allfiles = dir(folder);
for i = 1:length(allfiles)
    if allfiles(i).name(1) == '.'
        continue
    else
        break
    end
end
startfile = i;
allnames = NaN(length(allfiles)-startfile+1, 1);
for i = startfile:length(allfiles)
    allnames(i-startfile+1, 1:8) = allfiles(i).name(1:8);
end
daynames = unique(allnames, 'rows');
rowGall = 1;
rowCall = 1;
for i = 1:size(daynames,1)
    day = char(daynames(i,:));
    clear nearG farG nearPreG farPreG nearC farC nearPreC farPreC nearB farB nearPreB farPreB
    clear nearP farP nearPreP farPreP
    rowG = 1;
    rowC = 1;
    rowB = 1;
    rowP = 1;
    interstimG = [];
    interstimC = [];
    interstimB = [];
    interstimP = [];
    
    for j = startfile:size(allfiles,1)
        if allfiles(j).name(1:8) == day
            stro = nex2stro(strcat(folder, allfiles(j).name));
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
    end
    
    for k = 1:size(near,3)
        
        %raw data plot: FIRING RATE
        figure;
        set(gcf,'Name', day);
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
        set(gcf,'Name', day);
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

        
        
        
        
        
        
        
        
        
        %store radii of rings used for each unit
        radii(rowGall, :) = uniqueors/10;
        name(rowGall, :) = spikename(k,:);
        
        %store mean data and z-scored data for each unit (GRAY)
        Gnear(rowGall, :) = mean(nearG(:,:,k));
        Gfar(rowGall, :) = mean(farG(:,:,k));
        GnearResp(rowGall, :) = mean(nearRespG(:,:,k));
        GfarResp(rowGall, :) = mean(farRespG(:,:,k));
        GnearSEM(rowGall, :) = std(nearG(:,:,k))./sqrt(size(nearG,1));
        GfarSEM(rowGall, :)  = std(farG(:,:,k))./sqrt(size(farG,1));
        GnearRespSEM(rowGall, :) = std(nearRespG(:,:,k))./sqrt(size(nearRespG,1));
        GfarRespSEM(rowGall, :)  = std(farRespG(:,:,k))./sqrt(size(farRespG,1));
        
        Gznear(rowGall, :) = nanmean(nearGz(:,:,k));
        Gzfar(rowGall, :) = nanmean(farGz(:,:,k));
        GznearResp(rowGall, :) = nanmean(nearRespGz(:,:,k));
        GzfarResp(rowGall, :) = nanmean(farRespGz(:,:,k));
        GznearSEM(rowGall, :) = nanstd(nearGz(:,:,k))./sqrt(sum(~isnan(nearGz(:,:,k))));
        GzfarSEM(rowGall, :)  = nanstd(farGz(:,:,k))./sqrt(sum(~isnan(farGz(:,:,k))));
        GznearRespSEM(rowGall, :) = nanstd(nearRespGz(:,:,k))./sqrt(sum(~isnan(nearRespGz(:,:,k))));
        GzfarRespSEM(rowGall, :)  = nanstd(farRespGz(:,:,k))./sqrt(sum(~isnan(farRespGz(:,:,k))));
        
        rowGall = rowGall + 1;
        
        %store mean data and z-scored data for each unit (CORRIDOR)
        Cnear(rowCall, :) = mean(nearC(:,:,k));
        Cfar(rowCall, :) = mean(farC(:,:,k));
        CnearResp(rowCall, :) = mean(nearRespC(:,:,k));
        CfarResp(rowCall, :) = mean(farRespC(:,:,k));
        CnearSEM(rowCall, :) = std(nearC(:,:,k))./sqrt(size(nearC,1));
        CfarSEM(rowCall, :)  = std(farC(:,:,k))./sqrt(size(farC,1));
        CnearRespSEM(rowCall, :) = std(nearRespC(:,:,k))./sqrt(size(nearRespC,1));
        CfarRespSEM(rowCall, :)  = std(farRespC(:,:,k))./sqrt(size(farRespC,1));
        
        Cznear(rowCall, :) = nanmean(nearCz(:,:,k));
        Czfar(rowCall, :) = nanmean(farCz(:,:,k));
        CznearResp(rowCall, :) = nanmean(nearRespCz(:,:,k));
        CzfarResp(rowCall, :) = nanmean(farRespCz(:,:,k));
        CznearSEM(rowCall, :) = nanstd(nearCz(:,:,k))./sqrt(sum(~isnan(nearCz(:,:,k))));
        CzfarSEM(rowCall, :)  = nanstd(farCz(:,:,k))./sqrt(sum(~isnan(farCz(:,:,k))));
        CznearRespSEM(rowCall, :) = nanstd(nearRespCz(:,:,k))./sqrt(sum(~isnan(nearRespCz(:,:,k))));
        CzfarRespSEM(rowCall, :)  = nanstd(farRespCz(:,:,k))./sqrt(sum(~isnan(farRespCz(:,:,k))));
        
        rowCall = rowCall + 1;
   
    end
    
end









%plot of average of z-scores of FIRING RATES of all units
ordradii = 1:1:length(uniqueors); %ordinals for increasing ring size
figure
subplot(1, 2, 1)
hold on
errorbar(ordradii, mean(Gznear), std(Gznear)./sqrt(size(Gznear,1)), '.-m', 'linewidth', 2) 
errorbar(ordradii, mean(Gzfar), std(Gzfar)./sqrt(size(Gzfar,1)), '.-k', 'linewidth', 2) 
xlabel('Annulus radius ordinal');
ylabel('Firing rate (z-scored)');
title('Population: GRAY z-scored FIRING RATE data');
hold off
subplot(1, 2, 2)
hold on
errorbar(ordradii, mean(Cznear), std(Cznear)./sqrt(size(Cznear,1)), '.-m', 'linewidth', 2) 
errorbar(ordradii, mean(Czfar), std(Czfar)./sqrt(size(Czfar,1)), '.-k', 'linewidth', 2) 
xlabel('Annulus radius ordinal');
ylabel('Firing rate (z-scored)');
title('Population: CORRIDOR z-scored FIRING RATE data');
hold off

figure
hold on
errorbar(ordradii, mean(Gznear), std(Gznear)./sqrt(size(Gznear,1)), '.-m', 'linewidth', 2) 
errorbar(ordradii, mean(Gzfar), std(Gzfar)./sqrt(size(Gzfar,1)), '.-k', 'linewidth', 2) 
errorbar(ordradii, mean(Cznear), std(Cznear)./sqrt(size(Cznear,1)), '.-r', 'linewidth', 2) 
errorbar(ordradii, mean(Czfar), std(Czfar)./sqrt(size(Czfar,1)), '.-b', 'linewidth', 2) 
xlabel('Annulus radius ordinal');
ylabel('Firing rate (z-scored)');
title('Population: z-scored FIRING RATE data');
legend('Left', 'Right', 'Near', 'Far')
hold off

nocontext = Gznear;
sizeN = size(Gznear,1) + 1;
sizeF = size(Gzfar,1) + sizeN - 1;
nocontext(sizeN:sizeF, 1:7) = Gzfar;
figure
hold on
errorbar(ordradii, mean(nocontext), std(nocontext)./sqrt(size(nocontext,1)), '.-k', 'linewidth', 2) 
errorbar(ordradii, mean(Cznear), std(Cznear)./sqrt(size(Cznear,1)), '.-r', 'linewidth', 2) 
errorbar(ordradii, mean(Czfar), std(Czfar)./sqrt(size(Czfar,1)), '.-b', 'linewidth', 2) 
xlabel('Annulus radius ordinal');
ylabel('Firing rate (z-scored)');
title('Population: z-scored FIRING RATE data');
legend('No context', 'Near', 'Far')
hold off
















%plot of average of z-scores of RESPONSES of all units
figure
subplot(1, 2, 1)
hold on
errorbar(ordradii, mean(GznearResp), std(GznearResp)./sqrt(size(GznearResp,1)), '.-m', 'linewidth', 2)
errorbar(ordradii, mean(GzfarResp), std(GzfarResp)./sqrt(size(GzfarResp,1)), '.-k', 'linewidth', 2) 
xlabel('Annulus radius ordinal');
ylabel('Response (z-scored)');
title('Population: GRAY z-scored RESPONSE data');
hold off
subplot(1, 2, 2)
hold on
errorbar(ordradii, mean(CznearResp), std(CznearResp)./sqrt(size(CznearResp,1)), '.-m', 'linewidth', 2) 
errorbar(ordradii, mean(CzfarResp), std(CzfarResp)./sqrt(size(CzfarResp,1)), '.-k', 'linewidth', 2) 
xlabel('Annulus radius ordinal');
ylabel('Response (z-scored)');
title('Population: CORRIDOR z-scored RESPONSE data');
hold off

figure
hold on
errorbar(ordradii, mean(GznearResp), std(GznearResp)./sqrt(size(GznearResp,1)), '.-m', 'linewidth', 2) 
errorbar(ordradii, mean(GzfarResp), std(GzfarResp)./sqrt(size(GzfarResp,1)), '.-k', 'linewidth', 2) 
errorbar(ordradii, mean(CznearResp), std(CznearResp)./sqrt(size(CznearResp,1)), '.-r', 'linewidth', 2)
errorbar(ordradii, mean(CzfarResp), std(CzfarResp)./sqrt(size(CzfarResp,1)), '.-b', 'linewidth', 2)
xlabel('Annulus radius ordinal');
ylabel('Response (z-scored)');
title('Population: z-scored RESPONSE data');
legend('Left', 'Right', 'Near', 'Far')
hold off

nocontextResp = GznearResp;
sizeN = size(GznearResp,1) + 1;
sizeF = size(GzfarResp,1) + sizeN - 1;
nocontextResp(sizeN:sizeF, 1:7) = GzfarResp;
figure
hold on
errorbar(ordradii, mean(nocontextResp), std(nocontextResp)./sqrt(size(nocontextResp,1)), '.-k', 'linewidth', 2) 
errorbar(ordradii, mean(CznearResp), std(CznearResp)./sqrt(size(CznearResp,1)), '.-r', 'linewidth', 2)
errorbar(ordradii, mean(CzfarResp), std(CzfarResp)./sqrt(size(CzfarResp,1)), '.-b', 'linewidth', 2) 
xlabel('Annulus radius ordinal');
ylabel('Response (z-scored)');
title('Population: z-scored RESPONSE data');
legend('No context', 'Near', 'Far')
hold off





%{
save Units_FixInStay_Ap radii name ...
                        Gnear Gfar GnearResp GfarResp ...
                        GnearSEM GfarSEM GnearRespSEM GfarRespSEM ...
                        Gznear Gzfar GznearResp GzfarResp ...
                        GznearSEM GzfarSEM GznearRespSEM GzfarRespSEM ...
                        Cnear Cfar CnearResp CfarResp ...
                        CnearSEM CfarSEM CnearRespSEM CfarRespSEM ...
                        Cznear Czfar CznearResp CzfarResp ...
                        CznearSEM CzfarSEM CznearRespSEM CzfarRespSEM
      
%}     