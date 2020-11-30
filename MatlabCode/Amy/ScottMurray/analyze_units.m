%% Scott Murray's data analysis for single unit fits and population data

clear;close all
load Units_FixInMove_M2

%% Gray Background First

%% Normalize data to max in either the near or far condition
% Fit the data, y suffix is the fitted data, x suffix is for fit params

for i = 1:size(Gnear, 1)
    M = max([max(Gnear(i,:)) max(Gfar(i,:))]);
    Gnear(i,:) = Gnear(i,:)./M; 
    Gfar(i,:) = Gfar(i,:)./M;
    
    [Gray_Neary(i,:) Gray_Nearx(i,:)] = gauss_fit(radii(i,:),Gnear(i,:));
    [Gray_Fary(i,:) Gray_Farx(i,:)] = gauss_fit(radii(i,:),Gfar(i,:));
end

%% Scroll through and look at fits
for i = 1:size(Gnear,1)
    plot(radii(i,:),Gray_Neary(i,:),'m')
    hold on
    plot(radii(i,:),Gray_Fary(i,:),'k')
    legend ('Near','Far')
    set(gcf,'position',[920 502 560 420])
    plot(radii(i,:),Gnear(i,:),'m*')  
    plot(radii(i,:),Gfar(i,:),'k*')
    pause
    clf
end
%% Scatterplot of shift parameters
figure(1)
plot(Gray_Farx(:,2),Gray_Nearx(:,2),'*')
hold on
plot([1.5 3],[1.5 3])
xlabel('Far')
ylabel('Near')
axis([1.5 3 1.5 3])
title('Gray Background')
set(gcf,'position',[920 502 560 420])

%% Corridor Background

%% Normalize data to max in either the near or far condition
% Fit the data, y suffix is the fitted data, x suffix is for fit params

for i = 1:size(Cnear, 1)
    M = max([max(Cnear(i,:)) max(Cfar(i,:))]);
    Cnear(i,:) = Cnear(i,:)./M; 
    Cfar(i,:) = Cfar(i,:)./M;
    
    [Corr_Neary(i,:) Corr_Nearx(i,:)] = gauss_fit(radii(i,:),Cnear(i,:));
    [Corr_Fary(i,:) Corr_Farx(i,:)] = gauss_fit(radii(i,:),Cfar(i,:));
end

%% Scroll through and look at fits
figure(2)
for i = 1:size(Cnear,1)
    plot(radii(i,:),Corr_Neary(i,:),'r')
    hold on
    plot(radii(i,:),Corr_Fary(i,:),'b')
    legend ('Near','Far')
    set(gcf,'position',[920 502 560 420])
    plot(radii(i,:),Cnear(i,:),'r*')  
    plot(radii(i,:),Cfar(i,:),'b*')
    pause
    clf
end
%% Scatterplot of shift parameters
figure(2)
plot(Corr_Farx(:,2),Corr_Nearx(:,2),'*')
hold on
plot([1.5 3],[1.5 3])
xlabel('Far')
ylabel('Near')
axis([1.5 3 1.5 3])
title('Corridor Background')
set(gcf,'position',[920 502 560 420])

%% shift all of the radii by estimated far center
for i = 1:size(radii,1)
    R(i,:) = radii(i,:) - Gray_Farx(i,2);
end
%% calculate mean in bins for Gray
steps = -1:.2:1.2;
for i = 1 : length(steps)-1
    f=find(R< steps(i+1) & R >= steps(i));
    if length(f) > 0
        nearBin(i) = mean(Gnear(f));
        farBin(i) = mean(Gfar(f));
        midBin(i) = mean([steps(i+1) steps(i)]);
        meanBin(i) = mean(R(f));
    end
end
        
figure(3)
plot(midBin,nearBin,'m-o')
hold on
plot(midBin,farBin,'k-o')
axis([-1 1 0 1])
legend('Near','Far')
title('Gray Background')
set(gcf,'position',[920 502 560 420])

%% calculate mean in bins for Corridor
for i = 1 : length(steps)-1
    f=find(R< steps(i+1) & R >= steps(i));
    if length(f) > 0
        nearBin(i) = mean(Cnear(f));
        farBin(i) = mean(Cfar(f));
        midBin(i) = mean([steps(i+1) steps(i)]);
        meanBin(i) = mean(R(f));
    end
end

figure(4)
plot(midBin,nearBin,'r-o')
hold on
plot(midBin,farBin,'b-o')
axis([-1 1 0 1])
legend('Near','Far')
title('Corridor Background')
set(gcf,'position',[920 502 560 420])

