% Am writing this script to verify if I am thinking conceptually correctly
% I want to see if the isoresponse contour is a line even when the cone-signal combination is non-linear.
% Author - Abhishek De, 4/18
% The answer to this exercise is that the isoresponse contour in my
% experiment is dependent on both signal integration across space and
% cone-signal integration within each part of the space
close all; clearvars;
plot_counter = 1;
% Creating the M matrix from monitor SPDs and cone fundamentals
load fundamentals.mat
load mon_spd.mat
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]); %1st column - L, 2nd - M, 3rd - S
mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
mon_spd = spline([380:4:780], mon_spd', [380:5:780]); % fitting a cubic spline
M = fundamentals'*mon_spd'; % matrix that converts RGB phosphor intensites to L,M,S cone fundamentals
M = inv(M');
numcells = cell(3,1);
numcells{1}.S1wts = [0.5; 0.5; 0.5]; numcells{1}.S2wts = [-0.5; -0.5; -0.5]; % Luminance cell
numcells{2}.S1wts = [0.5; -0.5; 0]; numcells{2}.S2wts = [-0.5; 0.5; 0]; % L-M cell DO cell
numcells{3}.S1wts = [0.5; -0.5; -0.5]; numcells{3}.S2wts = [-0.5; 0.5; 0.5]; % L-M-S cell, orange-cyan DO cell
numcells{4}.S1wts = [-0.5; 0.5; -0.5]; numcells{4}.S2wts = [0.5; -0.5; 0.5]; % L-M+S cell, lime-magenta DO cell
numcells{5}.S1wts = [0.5; 0.5; -0.5]; numcells{5}.S2wts = [-0.5; -0.5; 0.5]; % -L-M+S cell, BY DO cell
% Here I will be creating the directions to be tested in the Subunit1-Subunit2 space
THETA = 0:10:360; THETA(end) = []; % in degrees
THETA = THETA*pi/180;
dirx = cos(THETA); diry = sin(THETA);
nreversals = 10;
TFR = 60; % in spikes/sec
changeincontrast = 0.1;
OOG  = 10;
mode = [1;2;3;4;5;6];

% Cone weights for each subunits
for kk = 1:numel(numcells)
    Subunit1wtsLMS = numcells{kk}.S1wts; % in LMS
    Subunit2wtsLMS = numcells{kk}.S2wts; % in LMS
    Subunit1wts = inv(M)*Subunit1wtsLMS; % in RGB
    Subunit2wts = inv(M)*Subunit2wtsLMS; % in RGB
    Subunit1wts = Subunit1wts/norm(Subunit1wts);
    Subunit2wts = Subunit2wts/norm(Subunit2wts);
    Subunit1light = Subunit1wts; % in terms of RGB
    Subunit2light = Subunit2wts; % in terms of RGB
    Subunit1light = Subunit1light/norm(Subunit1light);
    Subunit2light = Subunit2light/norm(Subunit2light);
    figure(plot_counter); set(gcf,'Name','Isoresponse contour: spatial integration');
    subplottitle = {'lin CI lin SI';'quad CI lin SI';'lin CI ellip SI';'quad CI ellip SI';'lin CI hype SI';'quad CI hype SI'};
    
    for jj = 1:numel(mode)
        stepsize_scalefactor = 0.75;
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
                FR = calcresponse(Subunit1wts,Subunit2wts,factor.*lightS1,factor.*lightS2,mode(jj));
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
        subplot(2,3,jj); plot(x(~oogidxs),y(~oogidxs),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1]); hold on;
        if sum(oogidxs)>0
            plot(upsample(x(oogidxs),2),upsample(y(oogidxs),2),'k');
        end
        xlabel('S1 contrast'); ylabel('S2 contrast'); grid on; axis equal; axis square;
        set(gca,'Xlim',[-7 7],'Ylim',[-7 7]); title(subplottitle{jj}); hold off;
    end
    plot_counter = plot_counter + 1;
end

