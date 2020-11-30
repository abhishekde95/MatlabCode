%%
%time course analysis
%compares neural response across different time analysis windows
%combines all blocks for one day's data set
%analyzes all days of data in a folder

%%
%{
%calculate spike data
%average firing rate for each small analysis window

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
    rowG = 1;
    rowC = 1;
    rowB = 1;
    
    %combining blocks of trials
    for j = startfile:size(allfiles,1)
        if allfiles(j).name(1:8) == day
            stro = nex2stro(strcat(folder, allfiles(j).name));
            [bkgnd, near, far, uniqueors, spikename] = SMurray_unit_spikes_all(stro); %spike rates per small analysis windows
            if bkgnd == 2 %gray background
                nearG(rowG:rowG+size(near,1)-1, :, :, :) = near;
                farG(rowG:rowG+size(far,1)-1, :, :, :) = far;
                rowG = rowG + size(near,1);
            elseif bkgnd == 0 %corridor background
                nearC(rowC:rowC+size(near,1)-1, :, :, :) = near;
                farC(rowC:rowC+size(far,1)-1, :, :, :) = far;
                rowC = rowC + size(near,1);
            elseif bkgnd == 1 %brick background
                nearB(rowB:rowB+size(near,1)-1, :, :, :) = near;
                farB(rowB:rowB+size(far,1)-1, :, :, :) = far;
                rowB = rowB + size(near,1);
            end
        end
    end
    
    for k = 1:size(near,4)
        %store radii of rings used for each unit
        radii(rowGall, :) = uniqueors/10;
        name(rowGall, :) = spikename(k,:);
        
        %store mean data and z-scored data for each unit (GRAY)
        Gnear(rowGall, :, :) = mean(nearG(:,:,:,k));
        Gfar(rowGall, :, :) = mean(farG(:,:,:,k));
        GnearSEM(rowGall, :, :) = std(nearG(:,:,:,k))./sqrt(size(nearG,1));
        GfarSEM(rowGall, :, :)  = std(farG(:,:,:,k))./sqrt(size(farG,1));
        rowGall = rowGall + 1;
        
        %store mean data and z-scored data for each unit (CORRIDOR)
        Cnear(rowCall, :, :) = mean(nearC(:,:,:,k));
        Cfar(rowCall, :, :) = mean(farC(:,:,:,k));
        CnearSEM(rowCall, :, :) = std(nearC(:,:,:,k))./sqrt(size(nearC,1));
        CfarSEM(rowCall, :, :)  = std(farC(:,:,:,k))./sqrt(size(farC,1));
        rowCall = rowCall + 1;
    end 
end

%}
%%
%{
%save data 

save Units_time_FixInMove_Ap radii name ...
                             Gnear Gfar GnearSEM GfarSEM ...
                             Cnear Cfar CnearSEM CfarSEM
 %}                        
%%
%load data

clear all; clc
load Units_time_FixInMove_Ap

%% 
%Normalize data to max in either the near or far condition

%Gray background
for i = 1:size(Gnear, 1) %each unit
    M = max([max(Gnear(i,:,:),[],3); max(Gfar(i,:,:),[],3)], [], 1);
    for j = 1:size(Gnear(i,:,:), 3)
        GnearNorm(i,:,j) = Gnear(i,:,j)./M;
        GfarNorm(i,:,j) = Gfar(i,:,j)./M;
    end
end

%Corridor background
for i = 1:size(Cnear, 1) %each unit
    M = max([max(Cnear(i,:,:),[],3); max(Cfar(i,:,:),[],3)], [], 1);
    for j = 1:size(Cnear(i,:,:), 3)
        CnearNorm(i,:,j) = Cnear(i,:,j)./M;
        CfarNorm(i,:,j) = Cfar(i,:,j)./M;
    end
end

%%
%Population average of normalized unit data

%Gray background
GnearPop = reshape(mean(GnearNorm,1),6,7); %rows: analysis time windows, columns: radii
GfarPop = reshape(mean(GfarNorm,1),6,7); 

%Corridor background
CnearPop = reshape(mean(CnearNorm,1),6,7); %rows: analysis time windows, columns: radii
CfarPop = reshape(mean(CfarNorm,1),6,7); 

%%
%plot each time analysis window

figure
set(gcf,'Name','50 ms analysis windows (starting from stim-onset)');
for i = 1:size(GnearPop,1)
    %Gray background
    subplot(2,size(GnearPop,1),i)
    hold on
    plot(radii(i,:), GnearPop(i,:), '.-b', 'linewidth', 2) 
    plot(radii(i,:), GfarPop(i,:), '.-r', 'linewidth', 2) 
    legend('Gray Near', 'Gray Far')
    hold off
    
    %Corridor background
    subplot(2,size(CnearPop,1),i+size(GnearPop,1))
    hold on
    plot(radii(i,:), CnearPop(i,:), '.-b', 'linewidth', 2) 
    plot(radii(i,:), CfarPop(i,:), '.-r', 'linewidth', 2) 
    legend('Corr Near', 'Corr Far')
    hold off
end