%analyze all data using function 'SMurray_folder_analysis'


%For ring size, choose to use ring radius (degrees) or percent of center radius (%)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%radmethod = 1; % Ring radius (degrees)
radmethod = 2; % Percent of center radius (%)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Choose method of calculating radius that gave peak response
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%peakmethod = 1; % Radius that gave the peak response
%peakmethod = 2; % X-axis value for peak of spline fit
peakmethod = 3; % X-axis value for peak of Gaussian fit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Choose to analyze unit data or LFP data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%whichdata = 1; % unit data
whichdata = 2; % LFP data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





if whichdata == 1 %unit data
    data = 'Units_FixInMove_M2'; %Fixation point in ring, moving, Monkey 2
    [FM2radii, FM2Gnear, FM2Gfar, FM2Cnear, FM2Cfar ...
        FM2popradii, FM2popGnear, FM2popGfar, FM2popCnear, FM2popCfar] = SMurray_folder_analysis(radmethod, peakmethod, data);
    data = 'Units_FixInStay_M2'; %Fixation point in ring, stationary, Monkey 2
    [FS2radii, FS2Gnear, FS2Gfar, FS2Cnear, FS2Cfar ...
        FS2popradii, FS2popGnear, FS2popGfar, FS2popCnear, FS2popCfar] = SMurray_folder_analysis(radmethod, peakmethod, data);
    data = 'Units_FixOutStay_M2'; %Fixation point outside of ring, stationary, Monkey 2
    [FO2radii, FO2Gnear, FO2Gfar, FO2Cnear, FO2Cfar ...
        FO2popradii, FO2popGnear, FO2popGfar, FO2popCnear, FO2popCfar] = SMurray_folder_analysis(radmethod, peakmethod, data);
    
    
    %scatterplot of peaks
    figure; hold on;
    minradii = [min(min(FM2radii)); min(min(FS2radii)); min(min(FO2radii))];
    maxradii = [max(max(FM2radii)); max(max(FS2radii)); max(max(FO2radii))];
    axis([min(minradii) max(maxradii) min(minradii) max(maxradii)]);
    axis square
    xlabel('Far: peak of size tuning curve (deg)')
    ylabel('Near: peak of size tuning curve (deg)')
    title('Single unit data: Monkey 2')
    plot(FM2Gfar, FM2Gnear, '+k','MarkerSize',8);
    plot(FM2Cfar, FM2Cnear, '+r','MarkerSize',8);
    plot(FS2Gfar, FS2Gnear, 'ok','MarkerSize',8);
    plot(FS2Cfar, FS2Cnear, 'or','MarkerSize',8);
    plot(FO2Gfar, FO2Gnear, 'xk','MarkerSize',8);
    plot(FO2Cfar, FO2Cnear, 'xr','MarkerSize',8);
    legend('Fix Move M2: Gray', 'Fix Move M2: Corridor', ...
           'Fix Stay M2: Gray', 'Fix Stay M2: Corridor', ...
           'Fix Out M2: Gray', 'Fix Out M2: Corridor', ...
           'Location', 'SouthEast')
    plot([min(minradii) max(maxradii)], [min(minradii) max(maxradii)], '-k')
    hold off
    
    
    %histogram of differences in peaks
    FM2Gdiff = FM2Gnear - FM2Gfar;
    FM2Cdiff = FM2Cnear - FM2Cfar;
    FS2Gdiff = FS2Gnear - FS2Gfar;
    FS2Cdiff = FS2Cnear - FS2Cfar;
    FO2Gdiff = FO2Gnear - FO2Gfar;
    FO2Cdiff = FO2Cnear - FO2Cfar;
    Gdiff = [FM2Gdiff; FS2Gdiff; FO2Gdiff];
    Cdiff = [FM2Cdiff; FS2Cdiff; FO2Cdiff];
    figure; subplot(1,2,1); hold on;
    axis([min([Gdiff; Cdiff]) max([Gdiff; Cdiff]) 0 Inf]); 
    histbins = (min([Gdiff; Cdiff]): 2 :max([Gdiff; Cdiff]));
    hist(Gdiff, histbins);
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor','k','EdgeColor','w')
    xlabel('Difference in spline peaks, Near - Far (deg)')
    ylabel('Units')
    title('Single unit: Gray background')
    hold off
    subplot(1,2,2); hold on;
    axis([min([Gdiff; Cdiff]) max([Gdiff; Cdiff]) 0 Inf]);
    hist(Cdiff, histbins);
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor','r','EdgeColor','w')
    xlabel('Difference in spline peaks, Near - Far (deg)')
    ylabel('Units')
    title('Single unit: Corridor background')
    hold off
 
    

elseif whichdata == 2 %LFP data
    data = 'LFP_FixInMove_M1'; %Fixation point in ring, moving, Monkey 1
    [FM1radii, FM1Gnear, FM1Gfar, FM1Cnear, FM1Cfar ...
        FM1popradii, FM1popGnear, FM1popGfar, FM1popCnear, FM1popCfar] = SMurray_folder_analysis(radmethod, peakmethod, data);
    data = 'LFP_FixInMove_M2'; %Fixation point in ring, moving, Monkey 2
    [FM2radii, FM2Gnear, FM2Gfar, FM2Cnear, FM2Cfar ...
        FM2popradii, FM2popGnear, FM2popGfar, FM2popCnear, FM2popCfar] = SMurray_folder_analysis(radmethod, peakmethod, data);
    data = 'LFP_FixInStay_M1'; %Fixation point in ring, stationary, Monkey 1
    [FS1radii, FS1Gnear, FS1Gfar, FS1Cnear, FS1Cfar ...
        FS1popradii, FS1popGnear, FS1popGfar, FS1popCnear, FS1popCfar] = SMurray_folder_analysis(radmethod, peakmethod, data);
    data = 'LFP_FixInStay_M2'; %Fixation point in ring, stationary, Monkey 2
    [FS2radii, FS2Gnear, FS2Gfar, FS2Cnear, FS2Cfar ...
        FS2popradii, FS2popGnear, FS2popGfar, FS2popCnear, FS2popCfar] = SMurray_folder_analysis(radmethod, peakmethod, data);
    data = 'LFP_FixOutStay_M2'; %Fixation point outside of ring, stationary, Monkey 2
    [FO2radii, FO2Gnear, FO2Gfar, FO2Cnear, FO2Cfar ...
        FO2popradii, FO2popGnear, FO2popGfar, FO2popCnear, FO2popCfar] = SMurray_folder_analysis(radmethod, peakmethod, data);
    
    
    %scatterplot of peaks
    figure; subplot(1,2,1); hold on;
    minradii1 = [min(min(FM1radii)); min(min(FS1radii))];
    maxradii1 = [max(max(FM1radii)); max(max(FS1radii))];
    axis([min(minradii1)*.8 max(maxradii1)*1.2 min(minradii1)*.8 max(maxradii1)*1.2]);
    axis square
    xlabel('Far: peak of size tuning curve (deg)')
    ylabel('Near: peak of size tuning curve (deg)')
    title('LFP data: Monkey 1')
    plot(FM1Gfar, FM1Gnear, '^k','MarkerSize',8);
    plot(FM1Cfar, FM1Cnear, '^r','MarkerSize',8);
    plot(FS1Gfar, FS1Gnear, '>k','MarkerSize',8);
    plot(FS1Cfar, FS1Cnear, '>r','MarkerSize',8);
    legend('Fix Move: Gray', 'Fix Move: Corridor', ...
           'Fix Stay: Gray', 'Fix Stay: Corridor', ...
           'Location', 'SouthEast')
    plot([min(minradii1)*.8 max(maxradii1)*1.2], [min(minradii1)*.8 max(maxradii1)*1.2], '-k')
    hold off
    
    subplot(1,2,2); hold on;
    minradii2 = [min(min(FM2radii)); min(min(FS2radii)); min(min(FO2radii))];
    maxradii2 = [max(max(FM2radii)); max(max(FS2radii)); max(max(FO2radii))];
    axis([min(minradii2)*.8 max(maxradii2)*1.2 min(minradii2)*.8 max(maxradii2)*1.2]);
    axis square
    xlabel('Far: peak of size tuning curve (deg)')
    ylabel('Near: peak of size tuning curve (deg)')
    title('LFP data: Monkey 2')
    plot(FM2Gfar, FM2Gnear, 'vk','MarkerSize',8);
    plot(FM2Cfar, FM2Cnear, 'vr','MarkerSize',8);
    plot(FS2Gfar, FS2Gnear, '<k','MarkerSize',8);
    plot(FS2Cfar, FS2Cnear, '<r','MarkerSize',8);
    plot(FO2Gfar, FO2Gnear, 'ok','MarkerSize',8);
    plot(FO2Cfar, FO2Cnear, 'or','MarkerSize',8);
    legend('Fix Move: Gray', 'Fix Move: Corridor', ...
           'Fix Stay: Gray', 'Fix Stay: Corridor', ...
           'Fix Out: Gray', 'Fix Out: Corridor', ...
           'Location', 'SouthEast')
    plot([min(minradii2)*.8 max(maxradii2)*1.2], [min(minradii2)*.8 max(maxradii2)*1.2], '-k')
    hold off
     
    
    %histogram of differences in peaks
    FM1Gdiff = FM1Gnear - FM1Gfar;
    FM1Cdiff = FM1Cnear - FM1Cfar;
    FS1Gdiff = FS1Gnear - FS1Gfar;
    FS1Cdiff = FS1Cnear - FS1Cfar;
    FM2Gdiff = FM2Gnear - FM2Gfar;
    FM2Cdiff = FM2Cnear - FM2Cfar;
    FS2Gdiff = FS2Gnear - FS2Gfar;
    FS2Cdiff = FS2Cnear - FS2Cfar;
    FO2Gdiff = FO2Gnear - FO2Gfar;
    FO2Cdiff = FO2Cnear - FO2Cfar;
    Gdiff = [FM1Gdiff; FS1Gdiff; FM2Gdiff; FS2Gdiff; FO2Gdiff];
    Cdiff = [FM1Cdiff; FS1Cdiff; FM2Cdiff; FS2Cdiff; FO2Cdiff];
    figure; subplot(1,2,1); hold on;
    axis([min([Gdiff; Cdiff]) max([Gdiff; Cdiff]) 0 Inf]); 
    histbins = (min([Gdiff; Cdiff]): 2 :max([Gdiff; Cdiff]));
    hist(Gdiff, histbins);
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor','k','EdgeColor','w')
    xlabel('Difference in spline peaks, Near - Far (deg)')
    ylabel('days')
    title('LFP: Gray background')
    hold off
    subplot(1,2,2); hold on;
    axis([min([Gdiff; Cdiff]) max([Gdiff; Cdiff]) 0 Inf]); 
    hist(Cdiff, histbins);
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor','r','EdgeColor','w')
    xlabel('Difference in spline peaks, Near - Far (deg)')
    ylabel('days')
    title('LFP: Corridor background')
    hold off

end
    


