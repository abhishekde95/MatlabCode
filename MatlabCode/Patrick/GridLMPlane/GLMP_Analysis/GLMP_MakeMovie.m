% This code is intended to make movies of GLMP datasets for JPW's General Exam.

function plat = GLMP_MakeMovie(plat,usenorm)

if nargin < 2
    usenorm = 0;
end

for p = 1%1:numel(plat)
    
    % User Defined Variables
    STARTAZANGLE = 45;
    ENDAZANGLE = 720;
    STARTELEVANGLE = 90;
    ENDELEVANGLE = 25;
    nViewShots = 500;
    
    
    % Plot figure
    if usenorm == 1
        x = plat{p}.par.Lcc_norm;
        y = plat{p}.par.Mcc_norm;
    else
        x = plat{p}.par.Lcc_orig;
        y = plat{p}.par.Mcc_orig;
    end
    z = plat{p}.par.meannspikes;
    
    % Linear Interpolation
    F = TriScatteredInterp(x,y,z);
    [qx,qy] = meshgrid(min(x):.01:max(x),min(y):.01:max(y));
    qz = F(qx,qy);
    
    % 2-D Naka-Rushton
    if ~isfield(plat{p},'TwoDFitAllData')
        plat = TwoDFit_AllData(plat);
    end
    
    if usenorm == 1
        xx = plat{p}.trial.Lcc_norm;
        yy = plat{p}.trial.Mcc_norm;
    else
        xx = plat{p}.trial.Lcc_orig;
        yy = plat{p}.trial.Mcc_orig;
    end

    zz = plat{p}.trial.nspikes;
    NRparams = plat{p}.TwoDFitAllData.Best.Params;
    rotMat = plat{p}.TwoDFitAllData.Best.RotMat;
    tempRotPts = [rotMat*[xx yy]']';
    xx = tempRotPts(:,1);
    yy = tempRotPts(:,2);
    if usenorm == 1
        axlim =max(plat{p}.trial.polrho_norm);
    else
        axlim = max(plat{p}.trial.polrho_orig);
    end
    [XX,YY] = meshgrid(linspace(-axlim,axlim,50));
    surface = ComputeNakaRushtonJPW(NRparams,[XX(:) YY(:)],'surface1');
    
    % Plot
    figure(1000+p); clf;
    axes('XTick',[],'YTick',[],'ZTick',[],'View',[0 STARTELEVANGLE])
    hold on;
    filename = plat{p}.datafile;
    plotTitle = [filename 'P' num2str(p)];
    set(gcf,'Name',plotTitle,'NumberTitle','off')    
    surfc(XX,YY,reshape(surface,size(XX)))
    plot3(xx,yy,zz,'*')    
    %xlim([min(x) max(x)]);
    %ylim([min(y) max(y)]);
    %axis equal square vis3d
    axis vis3d
    %subplot(1,3,3)
    %colorbar

    
    
    % Non-User defined variables
    viewangles = zeros(2,nViewShots);
    viewangles(1,:) = [repmat(STARTAZANGLE,1,nViewShots/5),linspace(STARTAZANGLE,ENDAZANGLE,4*nViewShots/5)];
    viewangles(2,:) = [linspace(STARTELEVANGLE, ENDELEVANGLE, nViewShots/4), repmat(ENDELEVANGLE,1,3*nViewShots/4)];
  %  keyboard
    for i = 1:length(viewangles)
        set(gca,'View',[viewangles(1,i) viewangles(2,i)])
        M(i) = getframe(gcf);
    end
    
    repeat = 1;     %default = 1
    pSearch = 1;    %default = 0
    bSearch = 1;    %default = 1
    reference = 1;  %default = 0
    pixRange = 10;  %default = 10
    iFrame = 8;     %default = 8
    pFrame = 10;    %default = 10
    bFrame = 25;    %default = 25
    options = [repeat, pSearch, bSearch, reference, pixRange, iFrame, pFrame, bFrame];
    mpgwrite(M, white, [plotTitle, '.mpg'], options);
    
%     figure(100); clf; hold on;
%     plot3(xx,yy,zz,'*')
%     surfc(XX,YY,reshape(surface,size(XX)))
%     colorbar('location','EastOutside')

   
end

