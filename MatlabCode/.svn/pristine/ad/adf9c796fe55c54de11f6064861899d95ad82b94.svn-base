% New script - daylight_props_7.m
% Author - Abhishek De, 1/18
% Part 1- Till now I have only built analysis asking how would the signals from the
% subunits of a DO combine signals but now I am asking how do the lights
% from a munsell edge itself change as a with change in illuminants
close all; clearvars;
plot_counter = 1;
wave = 400:10:720; % Taking this bizzare range to match the wavelength ranges of the natural images 
dayBasis = ieReadSpectra('cieDaylightBasis',wave); % Daylight spectra basis functions from isetbio
num_spectras = 100;
x = linspace(0.25,0.40,num_spectras);
y = 2.870*x - 3.000*(x.*x) - 0.275;
coeff1 = (-1.3515-1.7703*x+5.9114*y)./(0.0241+0.2562*x-0.7341*y);
coeff2 = (0.0300-31.4424*x+30.0717*y)./(0.0241+0.2562*x-0.7341*y);
coeffs = cat(2,ones(num_spectras,1),coeff1',coeff2'); % Limiting the coefficients between 0 and 1
illuminants = coeffs * dayBasis';

load fundamentals.mat;
load mon_spd;
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]); %1st column - L, 2nd- M, 3rd- S
mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);  % MONITOR SPECTRAL DISTRIBUTION IN R,G,B
mon_spd = spline([380:4:780], mon_spd', [380:5:780]); % fitting a cubic spline
lo = find(380:5:780==400);
hi = find(380:5:780==720);
fundamentals = fundamentals(lo:2:hi,:); % Starting the fundamentals from 390 nm
mon_spd = mon_spd(:,lo:2:hi);
M = fundamentals'*mon_spd';

load munsell380_800_1.mat % rows - wavelength, columns - different reflectances
lo1= find(380:1:780==400);
hi1 = find(380:1:780==720);
munsell = munsell(lo1:10:hi1,:)'; % all the reflectances are column vectors
totmunsells = size(munsell,1);
munsellsurfs = randi(totmunsells,[1000 2]);
RGBamp1 = cell(1,size(munsellsurfs,1)); 
RGBamp2 = cell(1,size(munsellsurfs,1));
RGBtriplet1 = cell(1,size(munsellsurfs,1));
RGBtriplet2 = cell(1,size(munsellsurfs,1));
LMSamp1 = cell(1,size(munsellsurfs,1));
LMSamp2 = cell(1,size(munsellsurfs,1));
LMStriplet1 = cell(1,size(munsellsurfs,1));
LMStriplet2 = cell(1,size(munsellsurfs,1));
ratioeigRGB = zeros(size(munsellsurfs,1),1);
ratioeigLMS = zeros(size(munsellsurfs,1),1);
Areargb = [];
Arealms = [];
plotting_option = 0; % 0- don't plot, 1- do plot
R = ceil(sqrt(size(munsellsurfs,1)));
for jj = 1:size(munsellsurfs,1)
    RGBamp1tmp = []; RGBamp2tmp = []; RGBtriplet1tmp = []; RGBtriplet2tmp = [];
    LMSamp1tmp = []; LMSamp2tmp = []; LMStriplet1tmp = []; LMStriplet2tmp = [];
    for ii = 1:size(illuminants,1)
        surf1 = munsell(munsellsurfs(jj,1),:); % reflectance spectra of surface 1
        surf2 = munsell(munsellsurfs(jj,2),:); % reflectance spectra of surface 2
        lightfromsurf1 = illuminants(ii,:).*surf1;
        lightfromsurf2 = illuminants(ii,:).*surf2;
        % RGB calculations
        RGBamp1tmp = [RGBamp1tmp; norm(mon_spd*lightfromsurf1')];
        RGBamp2tmp = [RGBamp2tmp; norm(mon_spd*lightfromsurf2')];
        RGBtriplet1tmp = [RGBtriplet1tmp mon_spd*lightfromsurf1'];
        RGBtriplet2tmp = [RGBtriplet2tmp mon_spd*lightfromsurf2'];
        % LMS calculations
        LMSamp1tmp = [LMSamp1tmp; norm(fundamentals'*lightfromsurf1')];
        LMSamp2tmp = [LMSamp2tmp; norm(fundamentals'*lightfromsurf2')];
        LMStriplet1tmp = [LMStriplet1tmp fundamentals'*lightfromsurf1'];
        LMStriplet2tmp = [LMStriplet2tmp fundamentals'*lightfromsurf2'];
    end
    RGBamp1{jj} = RGBamp1tmp; 
    RGBamp2{jj} = RGBamp2tmp;
    RGBtriplet1{jj} = RGBtriplet1tmp;
    RGBtriplet2{jj} = RGBtriplet2tmp;
    LMSamp1{jj} = LMSamp1tmp; 
    LMSamp2{jj} = LMSamp2tmp;
    LMStriplet1{jj} = LMStriplet1tmp;
    LMStriplet2{jj} = LMStriplet2tmp;
    
    try % in case the points are collinear
        [k,areatmprgb] = convhull(RGBamp1tmp,RGBamp2tmp); % getting the convex hull in the RGB amplitude space
    catch
        areatmprgb = 0;
    end
    try 
        [k,areatmplms] = convhull(LMSamp1tmp,LMSamp2tmp); % getting the convex hull in the LMS excitation space 
    catch 
                areatmplms = 0;
    end
    Areargb = [Areargb; areatmprgb];
    Arealms = [Arealms; areatmplms];
    if plotting_option
        figure(plot_counter);
        subplot(R,R,jj); plot(RGBamp1tmp,RGBamp2tmp,'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1]); hold on;
        plot(RGBamp1tmp(k),RGBamp2tmp(k),'r-');
        set(gca,'XTick',[],'YTick',[]); hold off; drawnow;
    end
    
    tmpvar = [RGBamp1tmp RGBamp2tmp]';
    crossprodRGB = tmpvar*tmpvar';
    [~,eigvalRGB] = eig(crossprodRGB);
    ratioeigRGB(jj) = max(diag(eigvalRGB))/min(diag(eigvalRGB));
    
    tmpvar = [LMSamp1tmp LMSamp2tmp]';
    crossprodLMS = tmpvar*tmpvar';
    [~,eigvalLMS] = eig(crossprodLMS);
    ratioeigLMS(jj) = max(diag(eigvalLMS))/min(diag(eigvalLMS));
    
    clear  RGBamp1tmp RGBamp2tmp RGBtriplet1tmp RGBtriplet2tmp
    clear  LMSamp1tmp LMSamp2tmp LMStriplet1tmp LMStriplet2tmp
end
[c1,p1] = corr(log10(Areargb),log10(ratioeigRGB));
[c2,p2] = corr(log10(Arealms),log10(ratioeigLMS));
plot_counter = plot_counter + 1;
figure(plot_counter),subplot(321);hist(log10(Areargb),30); xlabel('Area within convex hull'); ylabel('Frequency'); title('RGB');
subplot(322);hist(log10(Arealms),30); xlabel('Area within convex hull'); ylabel('Frequency'); title('LMS');
subplot(323),hist(log10(ratioeigRGB),30); xlabel('Ratio'), ylabel('Frequency'); title('Ratio of max to min eig val RGB');
subplot(324),hist(log10(ratioeigLMS),30); xlabel('Ratio'), ylabel('Frequency'); title('Ratio of max to min eig val LMS');
subplot(325),plot(log10(Areargb),log10(ratioeigRGB),'o','MarkerSize',2,'LineWidth',0.5,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1]); hold on; xlabel('Area'); ylabel('Ratio'); title('RGB');
subplot(326),plot(log10(Arealms),log10(ratioeigLMS),'o','MarkerSize',2,'LineWidth',0.5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0]); hold on; xlabel('Area'); ylabel('Ratio'); title('LMS');
subplot(325), hold on; lsline; hold off;
subplot(326), hold on; lsline; hold off;
plot_counter = plot_counter + 1;

% Sorting and plotting 5 most and least linear contour
[B,I] = sort(Areargb);
numcontour = 16;
mostlinear = I(1:numcontour);
leastlinear = I(end-numcontour-1:end);
numsubplots = ceil(sqrt(numcontour));
for ii = 1:numcontour
    % linear trajectories 
    figure(plot_counter), subplot(numsubplots,numsubplots,ii),plot(RGBamp1{mostlinear(ii)}(:),RGBamp2{mostlinear(ii)}(:),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1]);
    set(gca,'XTick',[],'YTick',[]); drawnow;
    figure(plot_counter+1), subplot(numsubplots,numsubplots,ii),plot(LMSamp1{mostlinear(ii)}(:),LMSamp2{mostlinear(ii)}(:),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[0 1 0]);
    set(gca,'XTick',[],'YTick',[]);drawnow;
    figure(plot_counter+2), subplot(numsubplots,numsubplots,ii), plot(munsell(munsellsurfs(mostlinear(ii),1),:),'b','Linewidth',2); hold on; plot(munsell(munsellsurfs(mostlinear(ii),2),:),'r','Linewidth',2);
    set(gca,'XTick',[],'YTick',[]);hold off; drawnow;
    
    % nonlinear trajectories
    figure(plot_counter+3), subplot(numsubplots,numsubplots,ii),plot(RGBamp1{leastlinear(ii)}(:),RGBamp2{leastlinear(ii)}(:),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1]);
    set(gca,'XTick',[],'YTick',[]);drawnow;
    figure(plot_counter+4), subplot(numsubplots,numsubplots,ii),plot(LMSamp1{leastlinear(ii)}(:),LMSamp2{leastlinear(ii)}(:),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[0 1 0]);
    set(gca,'XTick',[],'YTick',[]);drawnow;
    figure(plot_counter+5), subplot(numsubplots,numsubplots,ii), plot(munsell(munsellsurfs(leastlinear(ii),1),:),'b','Linewidth',2); hold on; plot(munsell(munsellsurfs(leastlinear(ii),2),:),'r','Linewidth',2);
    set(gca,'XTick',[],'YTick',[]);hold off; drawnow;
end
plot_counter = plot_counter + 6;

% plotting the RGB and LMS triplets
for ii = 1:5
    % RGB triplets
    figure(plot_counter),subplot(2,5,ii),plot3(RGBtriplet1{mostlinear(ii)}(1,:),RGBtriplet1{mostlinear(ii)}(2,:),RGBtriplet1{mostlinear(ii)}(3,:),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1]); hold on;
    plot3(RGBtriplet2{mostlinear(ii)}(1,:),RGBtriplet2{mostlinear(ii)}(2,:),RGBtriplet2{mostlinear(ii)}(3,:),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0]); title('Lin');hold off;
    figure(plot_counter),subplot(2,5,5+ii),plot3(RGBtriplet1{leastlinear(ii)}(1,:),RGBtriplet1{leastlinear(ii)}(2,:),RGBtriplet1{leastlinear(ii)}(3,:),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1]); hold on;
    plot3(RGBtriplet2{leastlinear(ii)}(1,:),RGBtriplet2{leastlinear(ii)}(2,:),RGBtriplet2{leastlinear(ii)}(3,:),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0]); title('Non-lin');hold off;
    % LMS triplets
    figure(plot_counter+1),subplot(2,5,ii),plot3(LMStriplet1{mostlinear(ii)}(1,:),LMStriplet1{mostlinear(ii)}(2,:),LMStriplet1{mostlinear(ii)}(3,:),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[0 1 0]); hold on;
    plot3(LMStriplet2{mostlinear(ii)}(1,:),LMStriplet2{mostlinear(ii)}(2,:),LMStriplet2{mostlinear(ii)}(3,:),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0]); title('Lin');hold off;
    figure(plot_counter+1),subplot(2,5,5+ii),plot3(LMStriplet1{leastlinear(ii)}(1,:),LMStriplet1{leastlinear(ii)}(2,:),LMStriplet1{leastlinear(ii)}(3,:),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[0 1 0]); hold on;
    plot3(LMStriplet2{leastlinear(ii)}(1,:),LMStriplet2{leastlinear(ii)}(2,:),LMStriplet2{leastlinear(ii)}(3,:),'o','MarkerSize',4,'LineWidth',0.5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0]); title('Non-lin');hold off;
end
plot_counter = plot_counter + 2;

%% Part 2 - This might be a follow up analysis to look at the cone excitation trajectories in the multi-dimensional space 

Volume = [];
LMStriplet = [LMStriplet1 LMStriplet2];
Crossprod = 0;
for ii = 1:numel(LMStriplet)
    try % in case the points are collinear
        [~,Volumetmp] = convhull(LMStriplet{ii}(1,:),LMStriplet{ii}(2,:),LMStriplet{ii}(3,:)); % getting the convex hull in the LMS cone excitation space
    catch
        Volumetmp = 0;
    end
    Volume = [Volume; Volumetmp];
    Crossprod = Crossprod + LMStriplet{ii}*LMStriplet{ii}';
end
figure(plot_counter); hist(Volume,30); 
xlabel('Volume'), ylabel('Frequency'),title('Volume enclosed by the hull'); hold off;
plot_counter = plot_counter + 1;
[vec,eigval] = eig(Crossprod);
% The result of this analysis is that most of the LMS excitation trajectories appear to be lines in this 2D space 

%% Part 3 - The next exercise to see what are the cone weights of the cells which will allow linear and non-linear computation
% This part is in continuation of the analysis done in part 1
numconditions = 2;
mean_cone_exc = fundamentals'*dayBasis(:,1);
cmap = [1 0 0; 0 1 0; 0 0 1];

figuretitle = ['linear'; 'nonlin'];
figure(plot_counter);
for ii = 1:numconditions
    disp(ii);
    cone_weights_V1_subunit = [];
    OC_cell = 0; % Orange-Cyan
    LM_cell = 0; % Lime-Magenta
    Lum_cell = 0; % Luminance cell
    BY_cell = 0; % Blue- Yellow Double opponent cell
    for jj = 1:100%size(munsellsurfs,1)
        surf1 = munsell(munsellsurfs(jj,1),:); % reflectance spectra of surface 1
        surf2 = munsell(munsellsurfs(jj,2),:); % reflectance spectra of surface 2
        [model,fval,success] = V1cellfit2(illuminants,surf1,surf2,fundamentals,mean_cone_exc,ii);
%         disp(jj);
        model = model./norm(model);
        cone_weights_V1_subunit = [cone_weights_V1_subunit; model'];
        L = cone_weights_V1_subunit(end,1);
        M = cone_weights_V1_subunit(end,2);
        S = cone_weights_V1_subunit(end,3);
        if sign(L) == -1*sign(M)
            if sign(S) == sign(L)
                LM_cell = LM_cell + 1;
            elseif sign(S) == sign(M)
                OC_cell = OC_cell + 1;
            end
        elseif sign(L) == sign(M)
            if sign(L) == sign(S)
                Lum_cell = Lum_cell + 1;
            elseif sign(L) == -1*sign(S)
                BY_cell = BY_cell + 1;
            end
        end
    end
    subplot(2,2,2*(ii-1)+1),plot3(cone_weights_V1_subunit(:,1),cone_weights_V1_subunit(:,2),cone_weights_V1_subunit(:,3),'o','MarkerSize',5,'LineWidth',0.1,'MarkerFaceColor',cmap(ii,:));
    set(gca,'Xlim',[-1 1],'YLim',[-1 1],'Zlim',[-1 1]); xlabel('L'), ylabel('M'), zlabel('S'); title('Cone Inputs to V1 subunit');
    line([-1 1],[0 0],[0 0],'LineWidth',4); line([0 0],[-1 1],[0 0],'LineWidth',4); line([0 0],[0 0],[-1 1],'LineWidth',4); grid on; hold off;
    subplot(2,2,2*(ii-1)+2),bar([LM_cell, OC_cell, Lum_cell, BY_cell]);set(gca,'XTick',[1 2 3 4],'XTickLabel',{'LM','OC','Lum','BY'}); title(figuretitle(ii,:));
end

plot_counter = plot_counter + 1;

