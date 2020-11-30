% Writing some new script to see how the cone-signal combination and
% spatial integration can be parsed out in this new simulation. This might
% be based on WNsimulate.m
% Author - Abhishek De, 4/18
close all; clearvars;
plot_counter = 1;
load mon_spd.mat % monitor spectral power distribution
load fundamentals.mat % cone fundamentals
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]); %1st column - L, 2nd- M, 3rd- S
mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
mon_spd = spline([380:4:780], mon_spd', [380:5:780]); % fitting a cubic spline
M = fundamentals'*mon_spd'; % matrix that converts RGB phosphor intensites to L,M,S cone fundamentals

Subunit1wtsLMS = [0.5; -0.5; 0.5]; % LMS cone wts
Subunit1wtsLMS = Subunit1wtsLMS/norm(Subunit1wtsLMS);
Subunit2wtsLMS = [-0.5; 0.5; -0.5]; % LMS cone wts
Subunit2wtsLMS = Subunit2wtsLMS/norm(Subunit2wtsLMS);
totframes = 300000;
RefreshRate = 75;  % Stim refresh rate (Hz)
dtbin = .1; % binsize for Poisson spike generation
x = linspace(-1,1,51);
[X,Y,Z] = meshgrid(x,x,x);
THETA = 0:5:360; THETA(end) = []; % in degrees
THETA = THETA*pi/180;
dirx = cos(THETA); diry = sin(THETA);
nreversals = 10;
TFR = 20; % in spikes/sec
changeincontrast = 0.1;
OOG  = 10;
stepsize_scalefactor = 0.75;
drive_multiplicativefactor = 5000;
Subunit1wts = inv(inv(M'))*Subunit1wtsLMS;
Subunit2wts = inv(inv(M'))*Subunit2wtsLMS;
Subunit1wts = Subunit1wts/norm(Subunit1wts);
Subunit2wts = Subunit2wts/norm(Subunit2wts);
input = randn(2*numel(Subunit1wts),totframes);
figuretitles = {'linear CI, linear SI';'quadratic CI, linear SI';'linear CI, elliptical SI';'quadratic CI, elliptical SI';'linear CI; hyperbolic SI';'quadratic CI, hyperbolic SI'};

for mode = 1:6
    % The drive depends on the cone weights and the LMS cone activations from the RGB lights
    drive  = getcumdrive3(Subunit1wtsLMS,Subunit2wtsLMS,M*input(1:3,:),M*input(4:6,:),mode);
    r = (max(0,drive_multiplicativefactor*drive)).^2; % spiking rectified quadratic nonlinearity
    rbig = repmat(r/RefreshRate*dtbin,1./dtbin,1); % make Poisson spike trail
    sp = sum(rand(size(rbig))<rbig)';
    nspikes = sum(sp);

    STA = mean(input(:,logical(sp)),2);
    k1(1,1,:) = Subunit1wts; k2(1,1,:) = Subunit2wts;
    k3(1,1,:) = STA(1:3); k4(1,1,:) = STA(4:6);
    truefilter = cat(2,repmat(k1,[100 50 1]),repmat(k2,[100 50 1]));
    retrievedfilter = cat(2,repmat(k3,[100 50 1]),repmat(k4,[100 50 1]));
    truefilter = 0.5 + 0.5*truefilter/(max(abs(truefilter(:)))+0.001);
    retrievedfilter = 0.5 + 0.5*retrievedfilter/(max(abs(retrievedfilter(:)))+0.001);
    STS = sum(input(:,logical(sp)),2);
    tmp = STS(:)*STS(:)';
    STCross = cov(input(:,logical(sp))');
    STCs = (nspikes.*STCross-tmp)/(nspikes*(nspikes-1));
    P = eye(size(STCs)) - STA(:)*inv(STA(:)'*STA(:))*STA(:)'; % WHAT DOES THIS LINE MEAN
    STCs = P*STCs*P';
    [tmp,d] = eig(STCs);
    [~,idx] = sort(diag(d)); % storing all the eigenvalues
    PC1 = tmp(:,idx(end)); PC2 = tmp(:,idx(end-1));
    k1(1,1,:) = PC1(1:3); k2(1,1,:) = PC1(4:6);
    k3(1,1,:) = PC2(1:3); k4(1,1,:) = PC2(4:6);
    PC1 = cat(2,repmat(k1,[100 50 1]),repmat(k2,[100 50 1]));
    PC2 = cat(2,repmat(k3,[100 50 1]),repmat(k4,[100 50 1]));
    PC1 = 0.5 + 0.5*PC1/(max(abs(PC1(:)))+0.001);
    PC2 = 0.5 + 0.5*PC2/(max(abs(PC2(:)))+0.001);
    
    % Need to do some kind of permutation test
    PC1eig = [];
    for jj = 1:1000
        sp1 = circshift(sp,randi(numel(sp)));
        STS1 = sum(input(:,logical(sp1)),2);
        tmp1 = STS1(:)*STS1(:)';
        STA1 = mean(input(:,logical(sp1)),2);
        STCross1 = cov(input(:,logical(sp1))');
        STCs1 = (nspikes.*STCross1-tmp1)/(nspikes*(nspikes-1));
        P1 = eye(size(STCs1)) - STA1(:)*inv(STA1(:)'*STA1(:))*STA1(:)'; % WHAT DOES THIS LINE MEAN
        STCs1 = P1*STCs1*P1';
        [tmp1,d1] = eig(STCs1);
        [~,idx1] = sort(diag(d1)); % storing all the eigenvalues
        PC1eig = [PC1eig; d1(idx1(end),idx1(end))];
    end
    
    % plotting the filters
    figure(plot_counter),set(gcf,'Name',figuretitles{mode});
    subplot(251);image(truefilter); set(gca,'XTick',[],'YTick',[]); title('True filter');axis square;
    subplot(252);image(retrievedfilter); set(gca,'XTick',[],'YTick',[]); title('STA'); axis square;
    subplot(253);image(PC1); set(gca,'XTick',[],'YTick',[]); title('PC1'); axis square;
    subplot(254);image(PC2); set(gca,'XTick',[],'YTick',[]); title('PC2'); axis square;
    subplot(255); plot(sort(diag(d)),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1]); hold on;
    plot([1 6],[prctile(PC1eig,99) prctile(PC1eig,99)],'k','Linewidth',2); ylabel('eigenvalues'); axis square; hold off;
    
    muS1spike = STA(1:3);
    muS2spike = STA(4:6);
    covS1spike = cov(input(1:3,logical(sp))'); covS1 = eye(3);
    covS2spike = cov(input(4:6,logical(sp))'); covS2 = eye(3);
    p1 = mvnpdf([X(:) Y(:) Z(:)],[0 0 0],covS1); p1 = reshape(p1,size(X));
    p2 = mvnpdf([X(:) Y(:) Z(:)],muS1spike',covS1spike); p2 = reshape(p2,size(X));
    p3 = mvnpdf([X(:) Y(:) Z(:)],[0 0 0],covS2); p3 = reshape(p3,size(X));
    p4 = mvnpdf([X(:) Y(:) Z(:)],muS2spike',covS2spike); p4 = reshape(p4,size(X));
    ratio1 = p2./p1;
    ratio2 = p4./p3;
    val1 = prctile(ratio1(:),[10 50 90]);
    val2 = prctile(ratio2(:),[10 50 90]);
    figure(plot_counter)
    for jj = 1:3
        fv1 = isosurface(X,Y,Z,ratio1,val1(jj));
        fv2 = isosurface(X,Y,Z,ratio2,val2(jj));
        subplot(256); plot3(fv1.vertices(:,1),fv1.vertices(:,2),fv1.vertices(:,3),'o','MarkerSize',2,'LineWidth',0.5); hold on;
        subplot(257); plot3(fv2.vertices(:,1),fv2.vertices(:,2),fv2.vertices(:,3),'o','MarkerSize',2,'LineWidth',0.5); hold on;
    end
    STA(1:3) = STA(1:3)/norm(STA(1:3));
    STA(4:6) = STA(4:6)/norm(STA(4:6));
    subplot(256); plot3([-STA(1) STA(1)],[-STA(2) STA(2)],[-STA(3) STA(3)],'-k','Linewidth',4); xlabel('R'); ylabel('G'); zlabel('B'); hold off;
    subplot(257); plot3([-STA(4) STA(4)],[-STA(5) STA(5)],[-STA(6) STA(6)],'-k','Linewidth',4); xlabel('R'); ylabel('G'); zlabel('B'); hold off;
    
    Subunit1light = Subunit1wts; % in terms of RGB
    Subunit2light = Subunit2wts; % in terms of RGB
    Subunit1light = Subunit1light/norm(Subunit1light);
    Subunit2light = Subunit2light/norm(Subunit2light);
    staircaseterminationpts = zeros(size(THETA));
    oogidxs = zeros(size(THETA));
    for ii = 1:numel(THETA)
        % need to implement a staircase procedure here
        lightS1 = repmat(dirx(ii),[3 1]).*Subunit1light;
        lightS2 = repmat(diry(ii),[3 1]).*Subunit2light;
        contrast = [1];
        resp = [];
        reversals = 0;
        dirspecificOOG = min([abs(OOG/dirx(ii)) abs(OOG/diry(ii))]);
        f = 1;
        while reversals <= nreversals
            factor = repmat(contrast(end),[3 1]);
            drive = getcumdrive3(Subunit1wtsLMS,Subunit2wtsLMS,M*(factor.*lightS1),M*(factor.*lightS2),mode);
            FR = (max(0,drive_multiplicativefactor*drive)).^2; % spiking non-linearity;
            resp = [resp; FR];
            if resp(end)<=TFR
                contrast = [contrast; contrast(end)*(1+f*changeincontrast)];
            else
                contrast = [contrast; contrast(end)*(1-f*changeincontrast)];
            end
            
            % checking for reversals
            if numel(resp)>=2
                if sign((resp(end)-TFR))*sign((resp(end-1)-TFR))== -1
                    reversals = reversals + 1;
                    f = f*stepsize_scalefactor;
                end
            end
            if contrast(end) >= dirspecificOOG
                oogidxs(ii) = 1;
                break;
            end
        end
        staircaseterminationpts(ii) = contrast(end);
    end
    oogidxs = logical(oogidxs);
    
    % Plotting the results
    [x,y] = pol2cart(THETA,staircaseterminationpts);
    figure(plot_counter); subplot(258), plot(x(~oogidxs),y(~oogidxs),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1]); hold on;
    if sum(oogidxs)>0
        plot(upsample(x(oogidxs),2),upsample(y(oogidxs),2),'k');
    end
    xlabel('S1 contrast'); ylabel('S2 contrast'); grid on; axis equal; axis square;
    set(gca,'Xlim',[-10 10],'Ylim',[-10 10]); hold off;

    % Simulating something similar to Patrick's experiment
    [Liso,Miso] = meshgrid(linspace(-50,50,101));
    Siso = zeros(size(Liso));
    LMstim = [Liso(:) Miso(:) Siso(:)]';
    driveS1 = getcumdrive3(Subunit1wtsLMS,Subunit2wtsLMS,LMstim,zeros(size(LMstim)),mode);
    FR1 = (max(0,drive_multiplicativefactor*driveS1)).^2;
    FR1 = reshape(FR1,size(Liso));
    driveS2 = getcumdrive3(Subunit1wtsLMS,Subunit2wtsLMS,zeros(size(LMstim)),LMstim,mode);
    FR2 = (max(0,drive_multiplicativefactor*driveS2)).^2;
    FR2 = reshape(FR2,size(Liso));
    figure(plot_counter); subplot(259); contourf(Liso,Miso,FR1); xlabel('M'); ylabel('L'); axis square;
    figure(plot_counter); subplot(2,5,10); contourf(Liso,Miso,FR2); xlabel('M'); ylabel('L'); axis square;
    plot_counter = plot_counter + 1;

end

