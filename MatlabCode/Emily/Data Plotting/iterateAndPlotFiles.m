%5/8/18 deprecated. Long since replaced by iterateAndPlotFiles_modularPlusDB
%all of Apollo's file names, in this case - working on getting rid of this
textFileList = {'ApolloTmpListLMTF.txt'};
%initializing
startpath = 'C:/NO BACKUP/NexFiles/Greg/Apollo'; %where to look for the files themselves
data = [];
isoData = [];
RFXYData = [];

%for each txt file, iterate through individual nex files
for iterator = 1:length(textFileList)
    % collect data files from individual text file
    fileName = textFileList{iterator};
    splitTitle = strsplit(fileName, '.');
    flist = flatten(fnamesFromTxt2(fullfile(nexfilepath,'nexfilelists', 'Greg','LMTF', fileName)));
    rfxyDataTemp = [];
    % put that data into the correct format using the functions
    % Iterate over the list of .nex files.
    for i = 1:length(flist)
        stro = notnex2stro(findfile(flist{i},startpath)); % <-- process the information in the nex file and put it into the "stro" structure
        % Below, just figuring out what information is in what column
        Lstim_idx = strcmp(stro.sum.trialFields(1,:), 'stim_idx');
        Llcc = strcmp(stro.sum.trialFields(1,:), 'lcc');
        Lmcc = strcmp(stro.sum.trialFields(1,:), 'mcc');
        Ltf = strcmp(stro.sum.trialFields(1,:), 'tf');
        
        % Getting the threshold points
        [stim_idxs,init_stim_trial_idxs] = unique(stro.trial(:,Lstim_idx),'last');
        questmodes = stro.trial(init_stim_trial_idxs, Llcc|Lmcc);
        tfs = stro.trial(init_stim_trial_idxs,Ltf);
        
        % Out of gamut checking
        funds = reshape(stro.sum.exptParams.fundamentals,length(stro.sum.exptParams.fundamentals)/3,3);
        if (size(stro.sum.exptParams.mon_spd,1) == 303)
            spds = SplineSpd([380:4:780]',reshape(stro.sum.exptParams.mon_spd,length(stro.sum.exptParams.mon_spd)/3,3),[380:5:780]');
        else
            spds = SplineSpd([380:2:780]',reshape(stro.sum.exptParams.mon_spd,length(stro.sum.exptParams.mon_spd)/3,3),[380:5:780]');
        end
        M = funds'*spds;
        bkgndrgb = stro.sum.exptParams.bkgndrgb;
        tmp = horzcat(stro.sum.exptParams.stim_x, stro.sum.exptParams.stim_y);
        rfxyDataTemp = [rfxyDataTemp; tmp];
        [in_gamut,scalar] = gamutCheck([questmodes zeros(size(questmodes,1),1)]', bkgndrgb, M, 'both');
        questmodes(~in_gamut) = questmodes(~in_gamut).*scalar(~in_gamut);
        data = [data; questmodes tfs ~in_gamut'];
    end
    RFXYData = transpose(unique(rfxyDataTemp));
    Loog = logical(data(:,end));
    
    % At this point "data" has n rows and 4 columns. The n rows correspond to n
    % different stimulus directions (L,M,TF combinations) that were tested. The
    % four columns are (1) L-cone contrast, (2) M-cone contrast, (3) Temporal
    % frequency, (4) out of gamut (1 or 0)
    % -----------------------------------------
    % Now we're transitioning from collecting the data from the various files
    % to fitting the data with the surface model.
    
    % This line is just re-representing L- and M-cone contrasts in polar coordinates.
    tf = data(:,3);
    LB = [0   0  1  1 .001 .001];
    UB = [100 1 10 10  .03  .03];
    x = data(:,1);
    y = data(:,2);
    
    initparams = [40 .1 5 3 .005 .002]; % An initial guess for the parameter values
    options = optimset('MaxFunEvals',50000,'MaxIter',50000,'TolFun',10^-8);
    thetas = linspace(0,pi/2,12);
    toterr = [];
    fpars = [];
 %   fpar = [initparams, initparams];
    PLOTINTERMEDIATES = 1;
    % Fitting the data many times (each time rotating the data by a small amount)
    for i = 1:length(thetas)
        rotmat = [cos(thetas(i)) -sin(thetas(i)); sin(thetas(i)) cos(thetas(i))];
        % Rotating data clockwise = rotating fit axis counterclockwise.
        xytmp = [x,y]*rotmat;
        
        % --- GDLH trying new, 1-D initial guesses ---
        [th,r] = cart2pol(xytmp(:,1),xytmp(:,2));
        L1 = mod(abs(th),pi) < .17; % a ~20 deg wedge centered on 0
        [fpar1,~] = fmincon(@(params) tf_fiterr(params,data(L1,3),1./r(L1)),initparams,[],[],[],[],LB,UB,[],options);
        L2 =  mod(abs(th-pi/2),pi) < .17; % a ~20 deg wedge centered on pi/2
        [fpar2,~] = fmincon(@(params) tf_fiterr(params,data(L2,3),1./r(L2)),initparams,[],[],[],[],LB,UB,[],options);
 
        if (PLOTINTERMEDIATES)
            omega = logspace(log10(1),log10(25),100);
            figure; axes; hold on;
            plot(data(L1,3),r(L1),'ko','MarkerFaceColor','black')
            plot(data(L2,3),r(L2),'ro','MarkerFaceColor','red')
            f1 = (1i*2*pi*10.^fpar1(5).*omega+1).^-fpar1(3);
            f2 = (1i*2*pi*10.^(fpar1(5)+fpar1(6)).*omega+1).^-(fpar1(3)+fpar1(4));
            f = fpar1(1)*abs(f1-fpar1(2)*f2);
            linestyle = '-';
            plot(omega,abs(1./f),linestyle,'linewidth',2,'color','black');
            f1 = (1i*2*pi*10.^fpar2(5).*omega+1).^-fpar2(3);
            f2 = (1i*2*pi*10.^(fpar2(5)+fpar2(6)).*omega+1).^-(fpar2(3)+fpar2(4));
            f = fpar2(1)*abs(f1-fpar2(2)*f2);
            plot(omega,1./f,linestyle,'linewidth',2,'color','red');
            drawnow;
        end
        % ---------------------------------------        
        % Now the 2-D fit
        options = optimset('Display', 'off'); %turns off the fmincon output
        [fpar,fv] = fmincon(@(params) tf_fiterr2(params,[xytmp(:,1) xytmp(:,2) tf],Loog),[fpar1 fpar2],[],[],[],[],[LB LB],[UB UB],[],options);
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
  %  if (thetas(bestrotidx) < pi/4) % GDLH ????
         fpar= fpar([[7:12]';[1:6]']);
        f1 = @(omega)fpar(1)*abs(((1i*2*pi*fpar(5).*omega+1).^-fpar(3))-fpar(2)*((1i*2*pi*fpar(6).*omega+1).^-fpar(4)));
        f2 = @(omega)fpar(1+6)*abs(((1i*2*pi*fpar(5+6).*omega+1).^-fpar(3+6))-fpar(2+6)*((1i*2*pi*fpar(6+6).*omega+1).^-fpar(4+6)));
    end
    %extract isodetection data and populate cell array
    isoData(iterator,:) = horzcat(RFXYData, fpar, bestrotidx);
    
    % plot data and name that plot as the txt file's name
    % Plotting the raw data
    figure; axes; hold on;
    plot3(x(~Loog),y(~Loog),tf(~Loog),'ko','MarkerFaceColor','black','MarkerSize',5);
    plot3(-x(~Loog,1),-y(~Loog),tf(~Loog),'ko','MarkerFaceColor','black','MarkerSize',5);
    plot3(x(Loog,1),y(Loog),tf(Loog),'ro','MarkerFaceColor','red','MarkerSize',5);
    plot3(-x(Loog,1),-y(Loog),tf(Loog),'ro','MarkerFaceColor','red','MarkerSize',5);
    
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
    
    % Plotting the fit
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
    set(gcf, 'name', splitTitle{1}, 'numbertitle', 'off');
    drawnow;
end    %move on to the next txt file in the list and start again

%Looking at resids from fit
predr = LMTF_thresh_from_model(data(:,1), data(:, 2),  data(:, 3), [fpar'; thetas(bestrotidx)]);
figure
plot(data(:,3), log10(predr) - log10(r), '.');
hold on;
plot([min(data(:,3)), max(data(:, 3))], [0,0]);
set(gca, 'xscale', 'log');
%3d plot of resids
[predL, predM] = pol2cart(atan2(data(:, 2), data(:, 1)), predr);
figure
resids = [predL(~Loog) - data(~Loog, 1), predM(~Loog) - data(~Loog, 2)];
plot3(predL - data(:, 1), predM - data(:, 2), data(:, 3), '.');
set(gca, 'zscale', 'log');
hold on;
plot3([0, 0], [0, 0], [min(data(:,3)), max(data(:, 3))]);
%plot3(-predL + data(:, 1), -predM + data(:, 2), data(:, 3), '.');




