 %% Long since deprecated and replaced by iterateAndPlotFiles_modularPlusDB
%Section 1: Initializing
textFileName = 'AbhishekLMTF.txt'; %edge case?
figTitle = strsplit(textFileName, '.');
startpath = '/Volumes/NO BACKUP/NexFiles/Abhishek'; %where to look for the files themselves
flist = flatten(fnamesFromTxt2(fullfile(nexfilepath,'nexfilelists','Greg', 'LMTF', textFileName)));
rfxyData = [];
data = [];
isoData = [];
tfMap = {};
%%
%Section 2: going through file list loop
for file = 1:length(flist)
    stro = notnex2stro(findfile(flist{file},startpath)); % <-- process the information in the nex file and put it into the "stro" structure
    % Below, just figuring out what information is in what column
    Lstim_idx = strcmp(stro.sum.trialFields(1,:), 'stim_idx');
    Llcc = strcmp(stro.sum.trialFields(1,:), 'lcc');
    Lmcc = strcmp(stro.sum.trialFields(1,:), 'mcc');
    Ltf = strcmp(stro.sum.trialFields(1,:), 'tf');
    
    %Getting the threshold points
    [stim_idxs,init_stim_trial_idxs] = unique(stro.trial(:,Lstim_idx),'last');
    questmodes = stro.trial(init_stim_trial_idxs, Llcc|Lmcc);
    tfs = stro.trial(init_stim_trial_idxs,Ltf);
    
    %Out of gamut checking
    funds = reshape(stro.sum.exptParams.fundamentals,length(stro.sum.exptParams.fundamentals)/3,3);
    if (size(stro.sum.exptParams.mon_spd,1) == 303)
        spds = SplineSpd([380:4:780]',reshape(stro.sum.exptParams.mon_spd,length(stro.sum.exptParams.mon_spd)/3,3),[380:5:780]');
    else
        spds = SplineSpd([380:2:780]',reshape(stro.sum.exptParams.mon_spd,length(stro.sum.exptParams.mon_spd)/3,3),[380:5:780]');
    end
    M = funds'*spds;
    bkgndrgb = stro.sum.exptParams.bkgndrgb;
    tmp = horzcat(stro.sum.exptParams.stim_x, stro.sum.exptParams.stim_y);
    rfxyData = [rfxyData; tmp];
    [in_gamut,scalar] = gamutCheck([questmodes zeros(size(questmodes,1),1)]', bkgndrgb, M, 'both');
    questmodes(~in_gamut) = questmodes(~in_gamut).*scalar(~in_gamut);
    data = [data; repmat(file,length(tfs),1) questmodes tfs ~in_gamut']; %all data
end
%%
%Section 3: creating a cell array of data
%sorted by Receptive Field location
uniqueRFXY = unique(rfxyData,'rows'); %unique Receptive Field
idxMap = NaN(size(rfxyData,1), size(uniqueRFXY, 1)); %map of indexed Receptive Field
for member = 1:size(uniqueRFXY,1)
    idxMap(:,member) = uniqueRFXY(member,1) == rfxyData(:,1) & uniqueRFXY(member,2) == rfxyData(:,2);
end
completeIndexedData = {size(uniqueRFXY,1)}; %complete indexed data
for index = 1:size(uniqueRFXY,1)
    currentColumn = idxMap(:,index);
    isDataMember = ismember(data(:,1), find(currentColumn));
    completeIndexedData{index} = data(isDataMember,:); 
end

%%
% Section 4: loop through each member of array to plot the data in individual files
% At this point "tempData" has n rows and 4 columns. The n rows correspond to n
% different stimulus directions (L,M,TF combinations) that were tested. The
% four columns are (1) L-cone contrast, (2) M-cone contrast, (3) Temporal
% frequency, (4) out of gamut (1 or 0)
% -----------------------------------------
% Now we're transitioning from collecting the data from the various files
% to fitting the data with the surface model(s).
for cellMember = 1:size(uniqueRFXY,1) 
    tempData = completeIndexedData{cellMember};
    %collecting temporal frequency information
    tf = tempData(:,4);
    upperTF = 19;
    lowerTF = 1;
    tfVector = tf <= upperTF & tf >= lowerTF; %collecting tfs within range
    tfMap{cellMember} = tempData(tfVector,[2,3,5]);
    
    % This line is just re-representing L- and M-cone contrasts in polar coordinates.
    LB = [0 0 1 1 .001 .001];
    UB = [100 1 20 20 .03 .03];
    x = tempData(:,2);
    y = tempData(:,3);
    miniLoog = tempData(:, end);
    
    initparams = [40 .1 9 3 .005 .002]; % An initial guess for the parameter values
    options = optimset('MaxFunEvals',50000,'MaxIter',50000,'TolFun',10^-8);
    thetas = linspace(0,pi/2,24);
    toterr = [];
    fpars = [];
    fpar = [initparams, initparams];
    %Fitting the data many times (each time rotating the data by a small amount)
    for i = 1:length(thetas)
        rotmat = [cos(thetas(i)) -sin(thetas(i)); sin(thetas(i)) cos(thetas(i))];
        % Rotating data clockwise = rotating fit axis counterclockwise.
        xytmp = [x,y]*rotmat;
        options = optimset('Display', 'off'); %turns off the fmincon output
        [fpar,fv] = fmincon(@(params) tf_fiterr2(params,[xytmp(:,1) xytmp(:,2) tf],miniLoog),fpar,[],[],[],[],[LB LB],[UB UB],[],options);
        toterr(i) = fv;
        fpars(i,:) = fpar;
    end
    bestrotidx = find(toterr == min(toterr));
    reshape(fpars(bestrotidx,:),[6 2])  % In case the user wants to see the fitted parameters
    fpar = fpars(bestrotidx,:);
    % fits to one of the several rotations of the data that were tried.
    f1 = @(omega)fpar(1)*abs(((1i*2*pi*fpar(5).*omega+1).^-fpar(3))-fpar(2)*((1i*2*pi*fpar(6).*omega+1).^-fpar(4)));
    f2 = @(omega)fpar(1+6)*abs(((1i*2*pi*fpar(5+6).*omega+1).^-fpar(3+6))-fpar(2+6)*((1i*2*pi*fpar(6+6).*omega+1).^-fpar(4+6)));
    if (f1(1) < f2(1)) % "If f1 is luminance" then exchange f1 and f2, forcing f1 to be chromatic
        fpar= fpar([[7:12]';[1:6]']);
        f1 = @(omega)fpar(1)*abs(((1i*2*pi*fpar(5).*omega+1).^-fpar(3))-fpar(2)*((1i*2*pi*fpar(6).*omega+1).^-fpar(4)));
        f2 = @(omega)fpar(1+6)*abs(((1i*2*pi*fpar(5+6).*omega+1).^-fpar(3+6))-fpar(2+6)*((1i*2*pi*fpar(6+6).*omega+1).^-fpar(4+6)));
    end
    %extract isodetection data
    isoRFXY = uniqueRFXY(cellMember,:);
    isoData(cellMember,:) = [isoRFXY fpar thetas(bestrotidx)];
    
    % plot data and name that plot as the txt file's name
    %'DEBUGTXT: Plotting the raw data'
    figure; axes; hold on;
    plot3(x(logical(~miniLoog)),y(logical(~miniLoog)),tf(logical(~miniLoog)),'ko','MarkerFaceColor','black','MarkerSize',5);
    plot3(-x(logical(~miniLoog),1),-y(logical(~miniLoog)),tf(logical(~miniLoog)),'ko','MarkerFaceColor','black','MarkerSize',5);
    plot3(x(logical(miniLoog),1),y(logical(miniLoog)),tf(logical(miniLoog)),'ro','MarkerFaceColor','red','MarkerSize',5);
    plot3(-x(logical(miniLoog),1),-y(logical(miniLoog)),tf(logical(miniLoog)),'ro','MarkerFaceColor','red','MarkerSize',5);
    
    % Adjusting the axes so we can see everything
    lim = max(max(abs([x y])));
    set(gca,'Xlim',1.1*[-lim lim]);
    set(gca,'Zscale','log');
    set(gca,'Ylim',1.1*[-lim lim]);
    axis square
    xlabel('L-cone contrast');
    ylabel('M-cone contrast');
    zlabel('TF (Hz)');
    set(gca,'View',[135 20]);
    axis vis3d
    
    %'DEBUGTXT: Plotting the fit'
    [xx,yy,zz] = meshgrid(linspace(-max(abs(x)),max(abs(x)),20),...
        linspace(-max(abs(y)),max(abs(y)),20),...
        linspace(min(tf),max(tf),20));
    a = abs(f1(zz)).^-1; % chromatic threshold
    b = abs(f2(zz)).^-1; % luminance threshold
    thtmp = atan2(yy,xx)+thetas(bestrotidx); % clockwise rotation from [L,M] to [a,b]
    rtmp = (a.*b)./sqrt((b.*cos(thtmp)).^2+(a.*sin(thtmp)).^2); % radius of ellipse - thank you, Wikipedia
    
    V = sqrt(xx.^2+yy.^2)-rtmp;
    FV = isosurface(xx,yy,zz,V,0);
    h = patch(FV);
    %set(h,'FaceColor','green','EdgeColor','none');
    set(h,'FaceColor','green','FaceAlpha',.5','EdgeColor','none','EdgeAlpha',0);
    set(gcf,'Renderer','painters');
    % Good views
    set(gca,'View',[135 12]);
    set(gca,'View',[225 12]);
    plotX = num2str(uniqueRFXY(cellMember,1));
    plotY = num2str(uniqueRFXY(cellMember,2));
    plotName = ['X: ' plotX ' and Y: ' plotY];
    set(gcf, 'name', plotName, 'numbertitle', 'off');
    title([figTitle{1}, ' ', plotX, ',', plotY]);
    drawnow;
end

%%
%Section 5: plot isodetection data
r = 2.5;
figure;
hold on;
for i = 1:size(isoData, 1)
    x0 = isoData(i,1);
    y0 = isoData(i,2);
    theta = isoData(i,end);
    x1 = x0 - r*(cos(theta));
    x2 = x0 + r*(cos(theta));
    y1 = y0 - r*(sin(theta));
    y2 = y0 + r*(sin(theta));
    plot([x1 x2], [y1 y2], 'k');
    drawnow;
end
set(gcf, 'name', 'isoData', 'numbertitle', 'off');
title([figTitle{1}, ' isoData']);

%%
%Section 6: plot data as modulated by temporal frequency
figure;
hold on;
tempOLmtx = [];
tempSLmtx = [];
tempOMmtx = [];
tempSMmtx = [];
for i = 1:size(uniqueRFXY, 1);
    shift = 20;
    %plot xy as a triangle
    tempRFX = uniqueRFXY(i,1);
    tempRFY = uniqueRFXY(i,2);
    plot(tempRFX, tempRFY, '*b');
    %plot each tf data point as xy dots around the triangles
    for j = 1:size(tfMap{i},1)
        %make two sets of numbers - original and symmetric - and add RFXY
        origL = round(tfMap{i}(j,1),4) + tempRFX; %L cone contrast is being treated as the X coord
        symmL = 2*(tempRFX) - origL;
        origM = round(tfMap{i}(j,2),4) + tempRFY; %M cone contrast as Y coord
        symmM = 2*(tempRFY) - origM;
        %add them to a temporary matrix so they can be plotted later
        tempOLmtx = [tempOLmtx; origL];
        tempSLmtx = [tempSLmtx; symmL];
        tempOMmtx = [tempOMmtx; origM];
        tempSMmtx = [tempSMmtx; symmM];
        %then plot them in the oog v. not oog colors using the following:
        if tfMap{i}(j,3)==0
           plot(origL, origM, '.k');
           plot(symmL, symmM, '.k');
        else
           plot(origL, origM, '.r');
           plot(symmL, symmM, '.r');
        end
    end
    %then fit a curve to the set of data points, but for real - FIXME
    x = [tempOLmtx; tempSLmtx]; %x-coords of ellipse
    y = [tempOMmtx; tempSMmtx]; %y-coords
    shiftedX = x - tempRFX; %shifting x coords so they are centered around 0 for ellipsefiterr
    shiftedY = y - tempRFY; %shifting y coords
    Loog = tfMap{i}(:,3);
    fullLoog = [Loog; Loog]; %to account for doubling x and y by adding OLmtx on top of SLmtx (for example)
    initparams = [std(shiftedX) std(shiftedY) acos((shiftedX./norm(shiftedX))'*shiftedY./norm(shiftedY))]; %first guess for fminsearch parameters - FIXME
    options = optimset('MaxFunEvals',50000,'MaxIter',50000,'TolFun',10^-5,'Display','none'); 
    minFnFit = fminsearch(@(params) ellipsefiterr(params,[shiftedX shiftedY],fullLoog),initparams, options); 
    angles = linspace(0,2*pi,100); 
    R = (minFnFit(2)^2-minFnFit(1)^2)*cos(2*angles-2*minFnFit(3))+minFnFit(1)^2+minFnFit(2)^2; 
    Q = sqrt(2)*minFnFit(1)*minFnFit(2)*sqrt(R); 
    r = Q./R; 
    [tmpx,tmpy] = pol2cart(angles,r); 
    plot(tempRFX+tmpx,tempRFY+tmpy,'m-'); %adds rfx/y to vectors of tmpx/y to shift back to original position around the rf
    drawnow;
end
 plotName = 'TF Data surrounding respective X"s and Y"s';
 set(gcf, 'name', plotName, 'numbertitle', 'off');
