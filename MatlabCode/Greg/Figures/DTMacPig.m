% Figures to the DTMacPig paper
%
% Section 1: Psychophysics. DTspot data. Detection data for blue and green
% phosphors as a function of retinal eccentricity.
%
% Section 1.1: As section 1, but calculating means and SDs for plotting
%
% Section 1.2: Predicted blue to green threshold ratios, for use with section
% 1 (or 1.1), above.
%
% Section 1.3: Comparison of blue to green threshold ratios at 15 and 25
% Hz. (from Apollo).
%
% Section 2: Comparison of three sets of cone fundamentals (2 deg, 10 deg
% and synth)
%
% Section 3: Example L+M planes with three degrees of macular pigment
% densities.
%
% Section 3.1: Example L+M planes with three lens densities
%
% Section 4: S-cone weight as a function of macular pigment density.
%
% Section 4.1: S-cone weights as a function of lens density.
%
% Section 5: S-cone weight as a function of eccentricity (converted to
% macular pigment density).
%
% Section 6: Example L+M ellipsoid with three degrees of macular pigment
% densities.
%
% Section 7: Direction of long axis of ellipsoids as a function of 
% lens density.
%
% Section 8: Distribution of S-cone weights among L+M cells with several
% different sets of cone fundamentals.
%
% Section 9: DTscot population analyses. Code largely lifted from
% DTScotAnalysis.m
%
% Section 9.05: DTscot population analysis. Normalized green threshold vs.
% normalized blue threshold.
% 
% Section 9.1: As above, but quantifying intensity in scotopic trolands.
%
% Section 9.2: barplot of scotopic B/G + inferred lens density
%
% Section 9.3: Scatterplot of scotopic B/G (inferred lens density) and
% photopic B/G (inferred lens density)
%
% Section 10: Comparing monkey luminance (from the literature) to human
% vlambda.
%
% Section 11: Comparing 15 Hz, 5� DTNT data across observers
%
% Section 11.1: Comparing 15 Hz, 5� DTNT data across observers, including
% directions of blue and green phosphor directions.
%
% Section 11.2: Bootstrap estimates of standard error of beta and
% calculating non-luminance contributions
%
% Section 12: 15 Hz, 5� DTNT surface slice through blue/green plane +
% macpig data
%
% Section 12.1: 15 Hz, 5� DTNT surface slice through blue/green plane +
% macpig data. Plotted in a space where blue is along one axis and green is
% along the orthogonal axis.
%
% Secton 12.2: 15 Hz 5� DTNT surface showing plane spaning the blue and
% green phosphor directions.
%
% Section 13: 25 Hz DTNT planes with different fundamentals
%
% Section 13.1: 25 Hz DTNT S-cone weights as a bar plot
%
% Section 14: Monitor spectra at different intensities (for Reviewer 2)
% -------------------
% Analyses that don't have accompanying figures
%
% Section A: Are the macaque B/G threshold ratios consistent with a scaled
% version of the human ratios? (Can we rule out the possibility that humans
% and monkeys have a similar spatial distribution of macular pigment?)
%
% Section B: Which is more sensitive to changes in L/M cone ratio: 
% green vs blue gun or red vs green gun?
%
% Section C: Computing the luminance of the ViewPixx monitor (and Dell 4)
%
% Section D: Does changing the peak of the S-cone absorbtion spectrum
% change our results at all?
%
% Section E: What lens density does Wald, 1945 see at 400 nm?
%
% Section F: Trying to reconcile Freya and Apollo's photopic and scotopic
% B/G ratios.
%
% Section G: Looking at the monitor spectra from the DTMacPig datasets. Was
% the blue phosphor unusually bright for some observers (e.g. Freya)?
%
% Section H: Added for Reviewer #2. Adjusting photopic B/G according to
% estimated lens density (from scotopic B/G). Does this decrease the
% variance?
%
% Section I: A big model looking at a variety of factors that could
% contribute to the photopic (and scotopic?) B/G with the approximately
% measured values. Factors are:
%     1) L:M ratio
%     2) S-cone contribution
%     3) Rod contribution
%     4) Lens optical density
%     5) photopigment optical density
%     6) Lambda max of S-cone (430 nm)
% Not even considering mac pig since we're only dealing with 7 deg ecc.
%
% Section J: Getting estimates of L and M from B/G ratios. Also how much
% do the non-dominant mechanisms contribute to the detection of B and G
% from the Cole fit.
%
% Section K: Looking at the peaks of photopigment optical densities from
% the CVRL website (Stockman). These are based on nomogram fits.
%
% Section J: Comparing Sedna's 25 Hz isdetection surface to isoresponse
% surfaces from L+M neurons. minimal preprocssing. Tey either align or they
% don't.
%%
% Section 1: Detection data for blue and green phosphors as a function of 
% retinal eccentricity. OUTDATED. Plots every point.

if (ispc)
    filelistpath = 'C:\NO BACKUP\NexFiles\nexfilelists\Greg\DTEM';
else
    filelistpath = '/Volumes/NO BACKUP/NexFiles/nexfilelists/Greg/DTEM';
end
HUMANS = 0;

if (HUMANS)
    listnames = {'GregMacPig.txt','ZackMacPig.txt','LeahMacPig.txt'};
else
    listnames = {'ApolloMacPig.txt'};
end

figure('Units','inches','PaperPosition',[0 0 8.5 11],'Position',[0 0 8.5 11],'PaperOrientation','Portrait','DefaultAxesUnits','inches');
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultTextFontName','Helvetica','DefaultTextFontSize',12);
colordef white;
load T_xyz1931;
Vlambda = T_xyz1931(2,:);
AXESWIDTH = 1.5;
MARGIN = .25;

h = [];
for i = 1:length(listnames)
    filenames = fnamesFromTxt2([filelistpath,filesep,char(listnames{i})]);
    data = [];
    for fileidx = 1:size(filenames,1)
        stro = nex2stro(findfile(char(filenames{fileidx,:})));
        [thresh, colorDirs, sfs] = DTquestUnpack(stro, 'mode'); close(gcf);
        normcolordirs = colorDirs./repmat(sqrt(sum(colorDirs.^2,2)),1,3);
        guns = reshape(stro.sum.exptParams.mon_spect,81,3);
        bkgndrgb = [stro.sum.exptParams.bkgnd_r stro.sum.exptParams.bkgnd_g stro.sum.exptParams.bkgnd_b];
        M = reshape(stro.sum.exptParams.m_mtx,3,3);
        bkgndlms = M*bkgndrgb';
        bkgndlum = bkgndrgb*guns'*Vlambda';
        threshrgb = inv(M)*(repmat(bkgndlms',2,1).*(1+repmat(thresh/100,[1 3]).*normcolordirs([1 2],:)))';
        thresholds = (threshrgb'*guns'*Vlambda')./bkgndlum;
        set(gcf,'Name',stro.sum.fileName);
        ecc = stro.sum.exptParams.rf_x;
        data = [data; ecc/10 thresholds'];
    end

    % Panel 1
    axes('position',[1+(i-1)*(AXESWIDTH+MARGIN) 4 AXESWIDTH AXESWIDTH]); hold on;
    plot(data(:,1),(data(:,2)-1)*100,'go','MarkerFaceColor','green','MarkerSize',3);
    plot(data(:,1),(data(:,3)-1)*100,'bo','MarkerFaceColor','blue','MarkerSize',3);
    if (i == 1);
        ylabel({'Detection threshold','(% contrast)'},'FontSize',12);
    end
    set(gca,'XLim',[0 8],'YLim',[1 12],'Yscale','log');
    set(gca,'YTick',[1 3 10],'YTickLabel',[1 3 10],'FontSize',10);
  
   % pp_g  = csaps(data(:,1),log10((data(:,2)-1)*100),.5);
   % pp_b  = csaps(data(:,1),log10((data(:,3)-1)*100),.5);
   % x = linspace(0,8,100);
   % plot(x,10.^(fnval(x,pp_g)),'g-','Linewidth',2)
   % plot(x,10.^(fnval(x,pp_b)),'b-','Linewidth',2)

    % Panel 2
    h(i) = axes('position',[1+(i-1)*(AXESWIDTH+MARGIN) 2 AXESWIDTH AXESWIDTH]); hold on;
    plot(data(:,1),(data(:,3)-1)./(data(:,2)-1),'ko','MarkerFaceColor','black','MarkerSize',3);
    [r,p] = corrcoef ([data(:,1),data(:,3)./data(:,2)])
    if (i == 1)
        ylabel({'Blue threshold','Green threshold'},'FontSize',12);
    end
    set(gca,'XLim',[0 8],'YLim',[0.4 1.25],'Yscale','log');
    set(gca,'YTick',[.5 .7 1])

    plot(x,10.^fnval(x,pp_b)./10.^fnval(x,pp_g),'k-','LineWidth',2)
    xlabel('Eccentricity (deg)','FontSize',12);
    if (HUMANS == 0)
        text(5,.3,['Monkey ',listnames{i}(1)],'fontsize',7);
    else
        text(5,.3,['Human ',listnames{i}(1)],'fontsize',7);
    end
end

%%
% Section 1.1: Detection data for blue and green phosphors as a function of 
% retinal eccentricity. Plotting means instead of raw data points.
% 2/11/14 GDLH Changing from 1931 V(lambda) to V*(lambda)

if (ispc)
    filelistpath = 'C:\NO BACKUP\NexFiles\nexfilelists\Greg\DTEM';
else
    filelistpath = '/Volumes/NO BACKUP/NexFiles/nexfilelists/Greg/DTEM';
end

figure('Units','inches','PaperPosition',[0 0 8.5 11],'Position',[0 0 8.5 11],'PaperOrientation','Portrait','DefaultAxesUnits','inches');
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultTextFontName','Helvetica','DefaultTextFontSize',12);
colordef white;
%load T_xyz1931;
%Vlambda = T_xyz1931(2,:);
load Vlambda_star
Vlambda = Vlambda_star(:,2)';

AXESWIDTH = 1.5;
XMARGIN = .25;
YMARGIN = .75;
plotcounter = 0;

corrdata = [];
bgdata = [];
for HUMANS = 0:1
    if (HUMANS)
        listnames = {'GregMacPig.txt','ZackMacPig.txt','LeahMacPig.txt'};
    else
        listnames = {'ApolloMacPig.txt','FreyaMacPig.txt','SednaMacPig.txt','KaliMacPig.txt','NutMacPig.txt',};
    end
    
    h = [];
    for i = 1:length(listnames)
        filenames = fnamesFromTxt2([filelistpath,filesep,char(listnames{i})]);
        
        data = [];
        for fileidx = 1:size(filenames,1)
            stro = nex2stro(findfile(char(filenames{fileidx,:})));
            [thresh, colorDirs, sfs] = DTquestUnpack(stro, 'mode'); close(gcf);
            normcolordirs = colorDirs./repmat(sqrt(sum(colorDirs.^2,2)),1,3);
            guns = reshape(stro.sum.exptParams.mon_spect,81,3);
            bkgndrgb = [stro.sum.exptParams.bkgnd_r stro.sum.exptParams.bkgnd_g stro.sum.exptParams.bkgnd_b];
            M = reshape(stro.sum.exptParams.m_mtx,3,3);
            bkgndlms = M*bkgndrgb';
            bkgndlum = bkgndrgb*guns'*Vlambda';
            threshrgb = inv(M)*(repmat(bkgndlms',2,1).*(1+repmat(thresh/100,[1 3]).*normcolordirs([1 2],:)))';
            deltagun = threshrgb-repmat(bkgndrgb',1,2);
            % Quantifying thresholds in 2 deg luminance contrast units
            thresholds = [nan; nan];
            for j = 1:size(deltagun,2)
                thresholds(j) =  (guns*deltagun(:,j))'*Vlambda';
                thresholds(j) =   thresholds(j)./bkgndlum;
            end
            ecc = stro.sum.exptParams.rf_x;
            data = [data; ecc/10 thresholds'];
        end
        
        % Collecting data for stats
        [r,p] = corr([data(:,1),data(:,3)./data(:,2)],'type','Spearman');
        corrdata = [corrdata; r(1,2) p(1,2)];
        bg = data(data(:,1) == 7,[2 3]);
        bgdata = [bgdata; mean(log10(bg(:,2)./bg(:,1)))]
        
        % Panel 1
        axes('position',[1+(3-mod(plotcounter,4))*(AXESWIDTH+XMARGIN) 3+floor(plotcounter/4)*(4+YMARGIN) AXESWIDTH AXESWIDTH]); hold on;
        %axes('position',[1+(i-1)*(AXESWIDTH+MARGIN) 3+4*HUMANS AXESWIDTH AXESWIDTH]); hold on;
        uniqueeccs = unique(data(:,1));
        for j = 1:length(uniqueeccs)
            L = logical(data(:,1) == uniqueeccs(j));
            mn = mean(data(L,2)*100);
            sem = std(data(L,2)*100)./sqrt(sum(L));
            plot(uniqueeccs(j),mn,'go','MarkerFaceColor','green','MarkerSize',3);
            plot([uniqueeccs(j) uniqueeccs(j)],mn+[-sem sem],'g-','LineWidth',1);
            mn = mean(data(L,3)*100);
            sem = std(data(L,3)*100)./sqrt(sum(L));
            plot(uniqueeccs(j),mn,'bo','MarkerFaceColor','blue','MarkerSize',3);
            plot([uniqueeccs(j) uniqueeccs(j)],mn+[-sem sem],'b-','LineWidth',1);
        end
        
        set(gca,'XLim',[0 8],'YLim',[1 12],'Yscale','log');
        if (mod(plotcounter,4)==3)
            ylabel({'Detection threshold','(% contrast)'},'FontSize',12);
            set(gca,'YTick',[1 3 10],'YTickLabel',[1 3 10],'FontSize',10);
        else
            set(gca,'YTickLabel',[]);
        end
        pp_g  = csaps(data(:,1),log10(data(:,2)*100),.5);
        pp_b  = csaps(data(:,1),log10(data(:,3)*100),.5);
        x = linspace(0,8,100);
        plot(x,10.^(fnval(x,pp_g)),'g-','Linewidth',2);
        plot(x,10.^(fnval(x,pp_b)),'b-','Linewidth',2);
        if (~HUMANS)
            set(gca,'Color',[.9 .9 .9]);
        end
        
        % ---- Panel 2 -----
        h(i) = axes('position',[1+(3-mod(plotcounter,4))*(AXESWIDTH+XMARGIN) 1+floor(plotcounter/4)*(4+YMARGIN) AXESWIDTH AXESWIDTH]); hold on;
        for j = 1:length(uniqueeccs)
            L = logical(data(:,1) == uniqueeccs(j));
            ratios = data(L,3)./data(L,2);
            mn = mean(log10(ratios));
            sem = std(log10(ratios))./sqrt(sum(L));
            plot(uniqueeccs(j),10.^mn,'ko','MarkerFaceColor','black','MarkerSize',3);
            plot([uniqueeccs(j) uniqueeccs(j)],10.^mn+[-sem sem],'k-','LineWidth',1);
        end
        if (mod(plotcounter,4)==3)
            ylabel({'Blue threshold','Green threshold'},'FontSize',12);
            set(gca,'YTick',[.5 .7 1])
        else
            set(gca,'YTickLabel',[]);
        end
        set(gca,'XLim',[0 8],'XTick',[0:2:8],'YLim',[0.5 1.5],'Yscale','log','FontSize',10);
        plot(x,10.^fnval(x,pp_b)./10.^fnval(x,pp_g),'k-','LineWidth',2)
        xlabel('Eccentricity (deg)','FontSize',12);
        if (10.^fnval(5,pp_b)./10.^fnval(5,pp_g) > 0.65)
            ypos = .55;
        else
            ypos = .7;
        end
        if (~HUMANS)
            text(5,ypos,['Monkey ',listnames{i}(1)],'fontsize',7);
            set(gca,'Color',[.9 .9 .9]);
        else
            text(5,ypos,['Human ',listnames{i}(1)],'fontsize',7);
        end
        
        % ---- Plotting predictions -----
        load('T_cones_myss2')
        spds = reshape(stro.sum.exptParams.mon_spect,81,3);
        M = T_cones_ss2'*spds;
        glum = Vlambda*spds(:,2);
        blum = Vlambda*spds(:,3);
        
        % 2 degree predictions
        LVLAMBDA = 1.98; % L contribution to vlambda - close enough to 1.98 which is what Sharpe correction says
        gpred = [LVLAMBDA 1]*M([1 2],2);
        bpred = [LVLAMBDA 1]*M([1 2],3);
        deg2pred = (gpred./bpred)./(glum./blum);
        plot([0 8],[deg2pred deg2pred],'r--');

        % 10 degree predictions
        load('T_cones_myss10')
        LVLAMBDA = 1; % L contribution to vlambda
        M = T_cones_ss10'*spds;
        gpred = [LVLAMBDA 1]*M([1 2],2);
        bpred = [LVLAMBDA 1]*M([1 2],3);
        deg10pred = (gpred./bpred)./(glum./blum);
        plot([0 8],[deg10pred deg10pred],'c--');
        
        % Synthetic predictions
        load('T_cones_synthgh2')
        M = T_cones_synthgh2'*spds;
        LVLAMBDA = 1; % L contribution to vlambda for "synthetic" predictions
        gpred = [LVLAMBDA 1]*M([1 2],2);
        bpred = [LVLAMBDA 1]*M([1 2],3);
        synthpred = (gpred./bpred)./(glum./blum);
        plot([0 8],[synthpred synthpred],'k--');
        plotcounter = plotcounter + 1;
        
        
    end
end

% debugging = a new set of synth preds.
load ('T_cones_myss10');
load ('den_lens_ss');
lenstransmittance = 1./(10.^(den_lens_ss));
fundnopig = T_cones_ss10./repmat(lenstransmittance,1,3);
fundnopig = fundnopig./repmat(max(fundnopig),size(fundnopig,1),1);
% Adding it back
lenstransmittance = 1./10.^(den_lens_ss./den_lens_ss(5)*1.5);
tmpfunds = fundnopig.*repmat(lenstransmittance,1,3);
tmpfunds = tmpfunds./repmat(max(tmpfunds),81,1);
M = tmpfunds'*spds;
LVLAMBDA = 1; % L contribution to vlambda for "synthetic" predictions
gpred = [LVLAMBDA 1]*M([1 2],2);
bpred = [LVLAMBDA 1]*M([1 2],3);
synthpred = (gpred./bpred)./(glum./blum) % .5565 with no lens
% 1.5 lens at 400 nm predicts B/G of 0.7 which is way too high!

%%
% % Section 1.2: Predictions of blum/glum
% % Need to run this right after Section 1, above
% Trying to figure outhow little (or how much negative)
% macular pigment density we need to get to the low B/G
% that we get to when we reduce the lens density to 1 at 400 nm.

% Standard synthetic prediction
load('T_cones_synthgh2')
M = T_cones_synthgh2'*spds;
LVLAMBDA = 1; % L contribution to vlambda for "synthetic" predictions
gpred = [LVLAMBDA 1]*M([1 2],2);
bpred = [LVLAMBDA 1]*M([1 2],3);
synthpred = (gpred./bpred)./(glum./blum);

% 10 degree predictions with zero macular pigment
% Stockman and Sharpe assume macular pigment has a peak density of 0.095
% for the 10� fundamentals
LVLAMBDA = 1;
load('T_cones_myss10')
load('den_mac_ss');
macpig = den_mac_ss;
macpigtransmittance = 1./(10.^(macpig./max(macpig)*0.095));
fundsnomacpig = T_cones_ss10./repmat(macpigtransmittance,1,3);
fundsnomacpig = fundsnomacpig./repmat(max(fundsnomacpig),81,1);
M = fundsnomacpig'*spds;
gpred = [LVLAMBDA 1]*M([1 2],2);
bpred = [LVLAMBDA 1]*M([1 2],3);
nomacpigpred = (gpred./bpred)./(glum./blum)

% negative macular pigment - finding answer by trial and error
macpigtransmittance = 1./(10.^(macpig./max(macpig)*-0.082));
funds = fundsnomacpig.*repmat(macpigtransmittance,1,3);
funds = funds./repmat(max(fundsnomacpig),81,1);
M = funds'*spds;
gpred = [LVLAMBDA 1]*M([1 2],2);
bpred = [LVLAMBDA 1]*M([1 2],3);
pred = (gpred./bpred)./(glum./blum)

% Need to reduce macular pigment by 0.095 to get (nominally) zero macular
% pigment in the 10� fundmentals. You have to reduce it by another 0.082
% to recapitulate the monkeys' B/G

%%
% Section 1.3
% Single panel - Apollo
% Comparison of blue to green threshold ratios at 15 and 25 Hz.

if (ispc)
    filelistpath = 'C:\NO BACKUP\NexFiles\nexfilelists\Greg\DTEM';
else
    filelistpath = '/Volumes/NO BACKUP/NexFiles/nexfilelists/Greg/DTEM';
end
listnames = {'ApolloMacPig.txt','ApolloMacPig25Hz.txt'};

h = []; data = [];
for i = 1:length(listnames)
    filenames = fnamesFromTxt2([filelistpath,filesep,char(listnames{i})]);
    for fileidx = 1:size(filenames,1)
        stro = nex2stro(findfile(char(filenames{fileidx,:})));
        [thresh, colorDirs, sfs] = DTquestUnpack(stro, 'mode'); close(gcf);
        normcolordirs = colorDirs./repmat(sqrt(sum(colorDirs.^2,2)),1,3);
        guns = reshape(stro.sum.exptParams.mon_spect,81,3);
        bkgndrgb = [stro.sum.exptParams.bkgnd_r stro.sum.exptParams.bkgnd_g stro.sum.exptParams.bkgnd_b];
        M = reshape(stro.sum.exptParams.m_mtx,3,3);
        bkgndlms = M*bkgndrgb';
        bkgndlum = bkgndrgb*guns'*Vlambda';
        threshrgb = inv(M)*(repmat(bkgndlms',2,1).*(1+repmat(thresh/100,[1 3]).*normcolordirs([1 2],:)))';
        deltagun = threshrgb-repmat(bkgndrgb',1,2);
        % Quantifying thresholds in 2 deg luminance contrast units
        thresholds = [nan; nan];
        for j = 1:size(deltagun,2)
            thresholds(j) =  (guns*deltagun(:,j))'*Vlambda';
            thresholds(j) =   thresholds(j)./bkgndlum;
        end
        ecc = stro.sum.exptParams.rf_x;
        TF = stro.trial(1,strcmp(stro.sum.trialFields(1,:),'gabor_speed'));
        data = [data; TF ecc/10 thresholds'];
    end
end

figure('Units','inches','PaperPosition',[0 0 8.5 11],'Position',[0 0 8.5 11],'PaperOrientation','Portrait','DefaultAxesUnits','inches');
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultTextFontName','Helvetica','DefaultTextFontSize',12);
colordef white;
load T_xyz1931;
Vlambda = T_xyz1931(2,:);
AXESWIDTH = 1.5;
MARGIN = .25;

% Panel 1
axes('position',[(AXESWIDTH+MARGIN) 4 AXESWIDTH AXESWIDTH]); hold on;
uniqueeccs = unique(data(:,2));
uniqueTFs = unique(data(:,1));
for i = length(uniqueTFs):-1:1
    LTF = logical(data(:,1) == uniqueTFs(i));
    ph = [];
    for j = 1:length(uniqueeccs)
        L = LTF & logical(data(:,2) == uniqueeccs(j));
        mn = mean(data(L,3)*100);
        sem = std(data(L,3)*100)./sqrt(sum(L));
        ph(1)=plot(uniqueeccs(j),mn,'go','MarkerFaceColor','green','MarkerSize',3);
        ph(2)=plot([uniqueeccs(j) uniqueeccs(j)],mn+[-sem sem],'g-','LineWidth',2);
        if (i == 2)
            set(ph,'Color',[.75 1 .75],'MarkerFaceColor',[.75 1 .75]);
        end
        mn = mean(data(L,4)*100);
        sem = std(data(L,4)*100)./sqrt(sum(L));
        ph(1)=plot(uniqueeccs(j),mn,'bo','MarkerFaceColor','blue','MarkerSize',3);
        ph(2)=plot([uniqueeccs(j) uniqueeccs(j)],mn+[-sem sem],'b-','LineWidth',2);
        if (i == 2)
            set(ph,'Color',[.75 .75 1],'MarkerFaceColor',[.75 .75 1]);
        end
    end
end
set(gca,'XLim',[0 8],'XTick',[0:2:8],'YLim',[1 12],'YTick',[1 3 10],'YTickLabel',[1 3 10],'Yscale','log','Color',[.9 .9 .9]);
set(gca,'YtickLabel',[]);

for i = length(uniqueTFs):-1:1
    LTF = logical(data(:,1) == uniqueTFs(i));
    pp_g(i) = csaps(data(LTF,2),log10(data(LTF,3)*100),.5);
    pp_b(i) = csaps(data(LTF,2),log10(data(LTF,4)*100),.5);
    x = linspace(0,8,100);
    ph=plot(x,10.^(fnval(x,pp_g(i))),'g-','Linewidth',2);
    if (i == 2)
        set(ph,'Color',[.75 1 .75],'MarkerFaceColor',[.75 1 .75]);
    end
    ph=plot(x,10.^(fnval(x,pp_b(i))),'b-','Linewidth',2);
    if (i == 2)
        set(ph,'Color',[.75 .75 1],'MarkerFaceColor',[.75 .75 1]);
    end
end

% Panel 2
h = axes('position',[AXESWIDTH+MARGIN 2 AXESWIDTH AXESWIDTH]); hold on;
for i = length(uniqueTFs):-1:1
    LTF = logical(data(:,1) == uniqueTFs(i));
    for j = 1:length(uniqueeccs)
        L = LTF & logical(data(:,2) == uniqueeccs(j));
        ratios = data(L,4)./data(L,3);
        mn = geomean(ratios);
        ph(1) = plot(uniqueeccs(j),mn,'ko','MarkerFaceColor','black','MarkerSize',3);
        ph(2) = plot(x,10.^fnval(x,pp_b(i))./10.^fnval(x,pp_g(i)),'k-','LineWidth',2)
        if (i == 2)
            set(ph,'Color',[.5 .5 .5],'MarkerFaceColor',[.5 .5 .5]);
        end
    end
end
set(gca,'XLim',[0 8],'XTick',[0:2:8],'YLim',[0.4 1.25],'Yscale','log','YTick',[.5 .7 1],'Color',[.9 .9 .9]);
set(gca,'YtickLabel',[]);


% ---- Plotting predictions -----
load('T_cones_myss2')
spds = reshape(stro.sum.exptParams.mon_spect,81,3);
M = T_cones_ss2'*spds;
glum = Vlambda*guns(:,2);
blum = Vlambda*guns(:,3);

% 2 degree predictions
LVLAMBDA = 1.5; % L contribution to vlambda (for 2 and 10 deg predictions)
gpred = [LVLAMBDA 1]*M([1 2],2);
bpred = [LVLAMBDA 1]*M([1 2],3);
deg2pred = (gpred./bpred)./(glum./blum);

% 10 degree predictions
load('T_cones_myss10')
M = T_cones_ss10'*spds;
gpred = [LVLAMBDA 1]*M([1 2],2);
bpred = [LVLAMBDA 1]*M([1 2],3);
deg10pred = (gpred./bpred)./(glum./blum);
plot([0 8],[deg2pred deg2pred],'k:');
plot([0 8],[deg10pred deg10pred],'k:');

% Synthetic predictions
load('T_cones_synthgh2')
M = T_cones_synthgh2'*spds;
LVLAMBDA = 1; % L contribution to vlambda for "synthetic" predictions
gpred = [LVLAMBDA 1]*M([1 2],2);
bpred = [LVLAMBDA 1]*M([1 2],3);
synthpred = (gpred./bpred)./(glum./blum);
plot([0 8],[synthpred synthpred],'k:');

% ----------------------------------


%% 
% Section 2
% A comparison of the three sets of cone fundamentals used in section 1

figure('Units','inches','PaperPosition',[0 0 8.5 11],'Position',[0 0 8.5 11],'PaperOrientation','Portrait','DefaultAxesUnits','inches');
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultTextFontName','Helvetica','DefaultTextFontSize',12);
colordef white;
axes('position',[2 2 4 4]); hold on;
load('T_cones_myss2');
load('T_cones_myss10');
load('T_cones_synthgh2');

wls = [380:5:780];
h = [];
h(:,3)=plot(wls,T_cones_synthgh2,'b-','LineWidth',3);
h(:,1)=plot(wls,T_cones_ss2,'k-','LineWidth',2);
h(:,2)=plot(wls,T_cones_ss10,'r-','LineWidth',2);
legend(h(1,:),{'SS 2 deg','SS 10 deg','Lens OD=1 @ 400 nm'},'location','southwest');
set(gca,'Yscale','log');
set(gca,'YLim',[10.^-6 1.5],'XLim',[380 780]);

%%
% Section 3
% Example planar L+M Neurothresh neuron

filename = 'K071309003.nex';  % planar luminance cell looks good
filename = 'K082109009.nex';  % planar luminance cell looks good
XLIM = .15; YLIM = .7;
AXESWIDTH = 2;
AXESMARGIN = .2;

load ('T_cones_synthgh2');
load ('den_mac_ss');
macpig = den_mac_ss;
macpigpeaks = [0 .11 .32];

% Getting data ready
stro = nex2stro(findfile(filename));
lmsidxs = [find(strcmp(stro.sum.trialFields(1,:),'lcont'))...
    find(strcmp(stro.sum.trialFields(1,:),'mcont'))...
    find(strcmp(stro.sum.trialFields(1,:),'scont'))];
% -------------------------------
% Converting cone contrasts in nex file to 10 deg fundmentals.
% Must go through excitations first - can't just transform contrasts.
% -------------------------------
fundamentals = stro.sum.exptParams.fundamentals;
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
mon_spd = stro.sum.exptParams.mon_spd;
mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
Moriginal = fundamentals'*mon_spd;
bkgndrgb = stro.sum.exptParams.bkgndrgb;
out = NTpreprocess(stro,0,Inf);
scaled = out(:,[2 3 4]).*repmat(out(:,5),1,3);
Loog = logical(out(:,7));
wls = [380:5:780];

figure('Units','inches','PaperPosition',[0 0 8.5 11],'Position',[0 0 8.5 11],'PaperOrientation','Portrait','DefaultAxesUnits','inches');
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultTextFontName','Helvetica','DefaultTextFontSize',12);
colordef white;

for i = 1:length(macpigpeaks)     
    % making new fundamentals
    macpigtransmittance = 1./(10.^(macpig./max(macpig).*macpigpeaks(i)));
    tmpfunds = T_cones_synthgh2.*repmat(macpigtransmittance,1,3);
    tmpfunds = tmpfunds./repmat(max(tmpfunds),81,1);
    Mnew = tmpfunds'*mon_spd;
    tmpscaled = ConvertConeContrastBasis(Moriginal, Mnew, bkgndrgb, scaled);
    Loog = logical(out(:,end));
    [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(tmpscaled, Loog);
    coneweights = (xformmat*planeparams)';
    if (sign(coneweights(1)) == -1)
        coneweights = -coneweights;
    end
    normconeweights = coneweights./sum(abs(coneweights)); 
    % Destructively changing tmpscaled here
    lmcoefs = mkbasis(coneweights([1 2]));
    xyscaled = tmpscaled(:,[1 2])*lmcoefs';
    xyscaled(:,2) = tmpscaled(:,3);
    
	axes('position',[1+(i-1)*(AXESWIDTH+AXESMARGIN) 3 AXESWIDTH AXESWIDTH])
    title(['Density @ 460 nm: ',num2str(macpigpeaks(i))]);
    hold on;
    axis square;

    % Plotting the data points
    h = plot(xyscaled(~Loog,1),xyscaled(~Loog,2),'ko');
    set(h,'Markersize',3,'Markerfacecolor','black');
    h = plot(-xyscaled(~Loog,1),-xyscaled(~Loog,2),'ko');
    set(h,'Markersize',3,'Markerfacecolor','black');
    xlabel([num2str(lmcoefs(1)./sum(lmcoefs),2),' L + ',num2str(lmcoefs(2)/sum(lmcoefs),2),' M']);
    set(get(gca,'XLabel'),'Color',[0 0 0]);
    if (i == 1)
        ylabel('S'); set(get(gca,'YLabel'),'Color',[0 0 0]);
    else
        set(gca,'YTick',[]);
    end
    
    % Plotting the OOG directions
    normfact = 1;
    h = plot([zeros(sum(Loog),1) xyscaled(Loog,1)./normfact]',[zeros(sum(Loog),1) xyscaled(Loog,2)./normfact]','-');
    set(h,'Color',[.8 .8 .8]);
    h = plot([zeros(sum(Loog),1) -xyscaled(Loog,1)./normfact]',[zeros(sum(Loog),1) -xyscaled(Loog,2)./normfact]','-');
    set(h,'Color',[.8 .8 .8]);
    
    % Adjusting the axes
    set(gca,'XLim',XLIM*[-1 1],'Ylim',YLIM*[-1 1]);
    
    % Plotting plane fits (which are actually lines in this projection)
    x = mean(abs(xyscaled))+[-.15 .15];
    y = (-norm(coneweights([1 2]))/coneweights(3))*x+(1/coneweights(3));
    h = plot(x,y,'g-','LineWidth',3);
    h = plot(-x,-y,'g-','LineWidth',3);
    text(.075,-.6,sprintf('S = %0.2f',normconeweights(:,3)),'FontSize',8)
end


%%
% Section 3.1
% An example luminance-tuned neuron from NeuroThresh rendered with various
% lens densities.

filename = 'S021010011.nex';  % planar luminance cell looks good
filename = 'K051311002.nex';
filename = 'K082109009.nex';  % planar luminance cell looks good
XLIM = .15; YLIM = .7;
AXESWIDTH = 2;
AXESMARGIN = .2;

load ('T_cones_myss10');
load ('den_lens_ss');
lenstransmittance = 1./(10.^(den_lens_ss));
fundnopig = T_cones_ss10./repmat(lenstransmittance,1,3);
fundnopig = fundnopig./repmat(max(fundnopig),size(fundnopig,1),1);
lensat400nm = [2 1 0];

% Getting data ready
stro = nex2stro(findfile(filename));
lmsidxs = [find(strcmp(stro.sum.trialFields(1,:),'lcont'))...
    find(strcmp(stro.sum.trialFields(1,:),'mcont'))...
    find(strcmp(stro.sum.trialFields(1,:),'scont'))];
% -------------------------------
% Converting cone contrasts in nex file to 10 deg fundmentals.
% Must go through excitations first - can't just transform contrasts.
% -------------------------------
fundamentals = stro.sum.exptParams.fundamentals;
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
mon_spd = stro.sum.exptParams.mon_spd;
mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
Moriginal = fundamentals'*mon_spd;
bkgndrgb = stro.sum.exptParams.bkgndrgb;
out = NTpreprocess(stro,0,Inf);
scaled = out(:,[2 3 4]).*repmat(out(:,5),1,3);
Loog = logical(out(:,7));
wls = [380:5:780];

figure('Units','inches','PaperPosition',[0 0 8.5 11],'Position',[0 0 8.5 11],'PaperOrientation','Portrait','DefaultAxesUnits','inches');
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultTextFontName','Helvetica','DefaultTextFontSize',12);
colordef white;

for i = 1:length(lensat400nm)     
    % making new fundamentals
    lenstransmittance = 1./(10.^(den_lens_ss/den_lens_ss(5).*lensat400nm(i)));
    tmpfunds = fundnopig.*repmat(lenstransmittance,1,3);
    tmpfunds = tmpfunds./repmat(max(tmpfunds),81,1);
    Mnew = tmpfunds'*mon_spd;
    tmpscaled = ConvertConeContrastBasis(Moriginal, Mnew, bkgndrgb, scaled);
    Loog = logical(out(:,end));
    [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(tmpscaled, Loog);
    coneweights = (xformmat*planeparams)';
    if (sign(coneweights(1)) == -1)
        coneweights = -coneweights;
    end
    normconeweights = coneweights./sum(abs(coneweights)); 
    % Destructively changing tmpscaled here
    lmcoefs = mkbasis(coneweights([1 2]));
    xyscaled = tmpscaled(:,[1 2])*lmcoefs';
    xyscaled(:,2) = tmpscaled(:,3);
    
	axes('position',[1+(i-1)*(AXESWIDTH+AXESMARGIN) 3 AXESWIDTH AXESWIDTH])
    title(['Density @ 400 nm: ',num2str(lensat400nm(i))]);
    hold on;
    axis square;

    % Plotting the data points
    h = plot(xyscaled(~Loog,1),xyscaled(~Loog,2),'ko');
    set(h,'Markersize',3,'Markerfacecolor','black');
    h = plot(-xyscaled(~Loog,1),-xyscaled(~Loog,2),'ko');
    set(h,'Markersize',3,'Markerfacecolor','black');
    xlabel([num2str(lmcoefs(1)./sum(lmcoefs),2),' L + ',num2str(lmcoefs(2)/sum(lmcoefs),2),' M']);
    set(get(gca,'XLabel'),'Color',[0 0 0]);
    if (i == 1)
        ylabel('S'); set(get(gca,'YLabel'),'Color',[0 0 0]);
    else
        set(gca,'YTick',[]);
    end
    
    % Plotting the OOG directions
    normfact = 1;
    h = plot([zeros(sum(Loog),1) xyscaled(Loog,1)./normfact]',[zeros(sum(Loog),1) xyscaled(Loog,2)./normfact]','-');
    set(h,'Color',[.8 .8 .8]);
    h = plot([zeros(sum(Loog),1) -xyscaled(Loog,1)./normfact]',[zeros(sum(Loog),1) -xyscaled(Loog,2)./normfact]','-');
    set(h,'Color',[.8 .8 .8]);
    
    % Adjusting the axes
    set(gca,'XLim',XLIM*[-1 1],'Ylim',YLIM*[-1 1]);
    
    % Plotting plane fits (which are actually lines in this projection)
    x = mean(abs(xyscaled))+[-.15 .15];
    y = (-norm(coneweights([1 2]))/coneweights(3))*x+(1/coneweights(3));
    h = plot(x,y,'g-','LineWidth',3);
    h = plot(-x,-y,'g-','LineWidth',3);
    text(.075,-.6,sprintf('S = %0.4f',normconeweights(:,3)),'FontSize',12)
end


%%
% Section 4)
% S-cone weight as a function of macular pigment density (using the
% synthetic cone fundamentals).

load ('T_cones_synthgh2');
fundnopig = T_cones_synthgh2;
wls = [380:5:780];

load ('den_mac_ss');
macpig = (den_mac_ss./max(den_mac_ss));
peak10deg = 0.095; % from SS paper
peak2deg = 0.35; % from SS paper

% ------- Debugging ---------
%load ('T_cones_smj10');
%mult10deg = 0.14; % from SMJ paper (Don't change, used below)
%macpigtransmittance = 1./(10.^(macpig.*mult10deg));
%fund10 = T_cones_smj10'./repmat(max(T_cones_smj10'),81,1);
%fundnomacpig = fund10./repmat(macpigtransmittance,1,3);
%T_cones_synthgh1 = fundnomacpig./repmat(max(fundnomacpig),81,1);
% ------- Debugging ---------

% First gathering all the data
data = [];
[fnames, spikeIdx] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\NT\NTlum2.txt');
for cellcounter = 1:size(fnames,1)
    filename = findfile(char(fnames{cellcounter}));
    paradigmID = getparadigmID(filename);
    if paradigmID == 103
        NT = nex2stro(filename);
    end
    out = NTpreprocess(NT,0,Inf);
    fundamentals = NT.sum.exptParams.fundamentals;
    fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
    mon_spd = NT.sum.exptParams.mon_spd;
    mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
    mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
    M = fundamentals'*mon_spd;
    Loog = logical(out(:,7));
    ecc = sqrt((NT.sum.exptParams.rf_x/10).^2 + (NT.sum.exptParams.rf_y/10).^2);
    data(cellcounter).M = M;
    data(cellcounter).data = out;
    data(cellcounter).bkgndrgb = NT.sum.exptParams.bkgndrgb;
    data(cellcounter).ecc = ecc;
    data(cellcounter).monspd = mon_spd;
end

% Now, computing cone weights as a function of macular pigment density
macpigpeaks = linspace(0,peak2deg,10);
for i = 1:length(data)
    i
    cc_original = data(i).data(:,[2:4]).*repmat(data(i).data(:,5),1,3);
    Loog = data(i).data(:,end);
    allconeweights = [];
    for j = 1:length(macpigpeaks)
        j
        macpigtransmittance = 1./(10.^(macpig.*macpigpeaks(j)));
        tmpfund = fundnopig.*repmat(macpigtransmittance,1,3);
        tmpfund = tmpfund./repmat(max(tmpfund),81,1);
        M = tmpfund'*mon_spd;
        scaled = ConvertConeContrastBasis(data(i).M, M, data(i).bkgndrgb, cc_original);
        [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(scaled, Loog);
        coneweights = (xformmat*planeparams)';
        if (sign(coneweights(2)) == -1)
            coneweights = -coneweights;
        end
        allconeweights = [allconeweights; coneweights];
    end
    data(i).coneweights = allconeweights;
end

% Computing means and standard deviations
allconeweights = cat(3,data(:).coneweights);
mn = []; stdev = [];
for j = 1:length(macpigpeaks)
    tmpconeweights = squeeze(allconeweights(j,:,:));
    normconeweights = tmpconeweights./repmat(sum(abs(tmpconeweights)),3,1);
    mn(j) = mean(normconeweights(3,:));
    stdev(j) = std(normconeweights(3,:));
end

% Doing the plotting
% At what macular pigment density do we have zero S-cone weight?

figure('Units','inches','PaperPosition',[0 0 8.5 11],'Position',[0 0 8.5 11],'PaperOrientation','Portrait','DefaultAxesUnits','inches');
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultTextFontName','Helvetica','DefaultTextFontSize',12);
colordef white;
axes('position',[2 2 4 4]); hold on;

tcrit = tinv(.025,length(mn)-1);
CI_upper = tcrit*stdev./sqrt(size(normconeweights,2));
h = errorbar(macpigpeaks,mn, stdev./sqrt(size(normconeweights,2)));
set(h,'LineWidth',2,'Color','black');
h = plot(macpigpeaks,mn,'ko','LineWidth',1,'MarkerFaceColor','black');
xlabel('Macular pigment density at 460 nm','FontSize',12);
ylabel('Mean S-cone weight','FontSize',12);
b = regress(mn',[macpigpeaks' ones(length(mn),1)]);

plot(-b(2)/b(1),0,'y*');
plot([0 -b(2)/b(1)],[0 0],'k-');
plot(-b(2)/b(1)*[1 1],[0 -.06],'k-');
set(gca,'XLim',[-.05 0.4],'YLim',[-.06 0.06])

% Finding where the CI crosses zero
blower = regress([mn-CI_upper]',[macpigpeaks' ones(length(mn),1)]);
bupper= regress([mn+CI_upper]',[macpigpeaks' ones(length(mn),1)]);
%plot(-blower(2)/blower(1)*[1 1],[-.06 0],'y-');
%plot(-bupper(2)/bupper(1)*[1 1],[-.06 0],'y-');

plot(-blower(2)/blower(1)*[1 1],(b(1)*blower(2)/blower(1))*[1 -1]+b(2),'y-');
plot(-bupper(2)/bupper(1)*[1 1],(b(1)*bupper(2)/bupper(1))*[1 -1]+b(2),'y-');

plot(peak10deg,-0.055,'rv','MarkerFaceColor','red','MarkerSize',6);
plot(peak2deg,-0.055,'gv','MarkerFaceColor','green','MarkerSize',6);

%%
% Section 4.1)
% S-cone weight as a function of lens density 
% (starting with the Stockman and Sharpe 10 degree fundamentals)

load ('T_cones_myss10.mat');
wls = [380:5:780];
load den_lens_ss;
lenstransmittance = 1./(10.^(den_lens_ss));
fundnopig = T_cones_ss10./repmat(lenstransmittance,1,3);
fundnopig = fundnopig./repmat(max(fundnopig),size(fundnopig,1),1);

load ('T_cones_smj10.mat');
wls = [380:5:780];
load den_lens_ss;
lenstransmittance = 1./(10.^(den_lens_ss));
fundnopig = T_cones_smj10./repmat(lenstransmittance',3,1);
fundnopig = fundnopig./repmat(max(fundnopig,[],2),1,size(fundnopig,2));


% First gathering all the data
data = [];
if(ispc)
    filelistpath = 'N:\NexFiles\nexfilelists\Greg\NT';
else
   filelistpath = '/Volumes/NO BACKUP/NexFiles/nexfilelists/Greg/NT';
end
[fnames, spikeIdx] = fnamesFromTxt2([filelistpath, filesep, 'NTlum2.txt']);
for cellcounter = 1:size(fnames,1)
    filename = findfile(char(fnames{cellcounter}));
    paradigmID = getparadigmID(filename);
    if paradigmID == 103
        NT = nex2stro(filename);
    end
    out = NTpreprocess(NT,0,Inf);
    fundamentals = NT.sum.exptParams.fundamentals;
    fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
    mon_spd = NT.sum.exptParams.mon_spd;
    mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
    mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
    M = fundamentals'*mon_spd;
    Loog = logical(out(:,7));
    ecc = sqrt((NT.sum.exptParams.rf_x/10).^2 + (NT.sum.exptParams.rf_y/10).^2);
    data(cellcounter).M = M;
    data(cellcounter).data = out;
    data(cellcounter).bkgndrgb = NT.sum.exptParams.bkgndrgb;
    data(cellcounter).ecc = ecc;
    data(cellcounter).monspd = mon_spd;
end

% Now, computing cone weights as a function of lens density
lensat400nm = linspace(0,2,21);
for i = 1:length(data)
    i
    cc_original = data(i).data(:,[2:4]).*repmat(data(i).data(:,5),1,3);
    Loog = data(i).data(:,end);
    allconeweights = [];
    for j = 1:length(lensat400nm)
        lenstransmittance = 1./(10.^(den_lens_ss/den_lens_ss(5).*lensat400nm(j)));
      %  lenstransmittance = 1./(10.^(den_lens_ss.*lensat400nm(j)));
        j
        tmpfund = fundnopig.*repmat(lenstransmittance,1,3);
        tmpfund = tmpfund./repmat(max(tmpfund),81,1);
        M = tmpfund'*data(i).monspd;
        scaled = ConvertConeContrastBasis(data(i).M, M, data(i).bkgndrgb, cc_original);
        [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(scaled, Loog);
        coneweights = (xformmat*planeparams)';
        if (sign(coneweights(2)) == -1)
            coneweights = -coneweights;
        end
        allconeweights = [allconeweights; coneweights];
    end
    data(i).coneweights = allconeweights;
end

% Doing the plotting
% At what lens density do we have zero S-cone weight?

figure('Units','inches','PaperPosition',[0 0 8.5 11],'Position',[0 0 8.5 11],'PaperOrientation','Portrait','DefaultAxesUnits','inches');
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultTextFontName','Helvetica','DefaultTextFontSize',12);
colordef white;

for whichplot = 1:2  % S-cone or L-M
    axes('position',[2 6.5-5*(whichplot-1) 3 3]); hold on;

    % Computing means and standard deviations
    allconeweights = cat(3,data(:).coneweights);
    mn = []; stdev = [];
    for j = 1:length(lensat400nm)
        tmpconeweights = squeeze(allconeweights(j,:,:));
        normconeweights = tmpconeweights./repmat(sum(abs(tmpconeweights)),3,1);
        if whichplot ==1
            mn(j) = mean(normconeweights(3,:));
            stdev(j) = std(normconeweights(3,:));
        else
            if (any(any(tmpconeweights([1 2],:) < 0)))
                mn(j) = nan;
                stdev(j) = nan;
            else
                mn(j) = mean(log10(tmpconeweights(1,:)./tmpconeweights(2,:)));
                stdev(j) = std(log10(tmpconeweights(1,:)./tmpconeweights(2,:)));
            end
        end
    end

    tcrit = tinv(.025,size(normconeweights,2)-1);
    CI_lower = tcrit*stdev./sqrt(size(tmpconeweights,2));
    %h = errorbar(lensat400nm,mn, stdev./sqrt(size(normconeweights,2)));
    %set(h,'LineWidth',2,'Color','black');
    %h = plot(lensat400nm,mn,'ko','LineWidth',1,'MarkerFaceColor','black');
    x = linspace(lensat400nm(1),lensat400nm(end),length(mn));
    Lposwts = ~isnan(mn);
    x(x<min(lensat400nm(Lposwts))) = []; 
    x(x>max(lensat400nm(Lposwts))) = [];

    % Plotting a smooth curve for the mean
    if (whichplot == 1)
        plot(x,interp1(lensat400nm(Lposwts),mn(Lposwts),x,'spline'),'k-','LineWidth',2);
        densityatzeroweight = interp1(mn,lensat400nm,0,'spline');
        plot([lensat400nm(1) densityatzeroweight],[0 0],'k-');
        plot([densityatzeroweight densityatzeroweight],[0 -.05],'k-');
    else
        plot(x,interp1(lensat400nm(Lposwts),10.^mn(Lposwts),x,'spline'),'k-','LineWidth',2);
        densityatzeroweight = interp1(10.^mn(Lposwts),lensat400nm(Lposwts),1,'spline');
        plot([lensat400nm(1) densityatzeroweight],[1 1],'k-');
        plot([densityatzeroweight densityatzeroweight],[1 .25],'k-');
    end
    set(gca,'XLim',[min(lensat400nm) max(lensat400nm)]);

    % Plotting confidence intervals
    if (whichplot == 1)
        % Plotting a smooth curve for the confidence interval
        plot(x,interp1(lensat400nm,mn+CI_lower,x,'spline'),'k-','LineWidth',1);
        plot(x,interp1(lensat400nm,mn-CI_lower,x,'spline'),'k-','LineWidth',1);
        densityatCIweight(1) = interp1(mn-CI_lower,lensat400nm,0,'spline');
        densityatCIweight(2) = interp1(mn+CI_lower,lensat400nm,0,'spline');
        for j = 1:2
            plot([lensat400nm(j) densityatCIweight(j)],[0 0],'k:');
            plot([densityatCIweight(j) densityatCIweight(j)],[0 -.05],'k:');
        end
        set(gca,'YLim',[-.05 0.04]);
        ylabel('S-cone weight','FontSize',12);
    else % L/M  cone weight ratio
        plot(x,interp1(lensat400nm(Lposwts),10.^(mn(Lposwts) + CI_lower(Lposwts)),x,'spline'),'k-','LineWidth',1);
        plot(x,interp1(lensat400nm(Lposwts),10.^(mn(Lposwts) - CI_lower(Lposwts)),x,'spline'),'k-','LineWidth',1);
        densityatCIweight(1) = interp1(10.^(mn(Lposwts) + CI_lower(Lposwts)),lensat400nm(Lposwts), 1,'spline');
        densityatCIweight(2) = interp1(10.^(mn(Lposwts) - CI_lower(Lposwts)),lensat400nm(Lposwts), 1,'spline');
        for j = 1:2
            plot([lensat400nm(j) densityatCIweight(j)],[1 1],'k:');
            plot([densityatCIweight(j) densityatCIweight(j)],[1 .025],'k:');
        end
        set(gca,'Ylim',[.25 4]);
        set(gca,'Yscale','log');
        ylabel({'L-weight/M-weight'},'FontSize',14);
        set(gca,'Ytick',[.5 1 2],'Yticklabel',[.5 1 2])
    end
    xlabel('Lens density at 400 nm','FontSize',14);
end

%%
% Section 5
% S-cone weights (converted to macular pigment density) as a function of
% eccentricity.

load ('den_mac_ss');
macpig = (den_mac_ss./max(den_mac_ss))';

% using a new (synthetic) set of cone fundamentals
load ('T_cones_synthgh2');
fundnopig = T_cones_synthgh2;

% First gathering all the data
data = [];
if (ispc)
    [fnames, spikeIdx] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\NT\NTlum2.txt');
else
    [fnames, spikeIdx] = fnamesFromTxt2('/Volumes/NO BACKUP/NexFiles/nexfilelists/Greg/NT/NTlum2.txt');
end

for cellcounter = 1:size(fnames,1)
    filename = findfile(char(fnames{cellcounter}));
    paradigmID = getparadigmID(filename);
    if paradigmID == 103
        NT = nex2stro(filename);
    end
    out = NTpreprocess(NT,0,Inf);
    fundamentals = NT.sum.exptParams.fundamentals;
    fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
    mon_spd = NT.sum.exptParams.mon_spd;
    mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
    mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
    M = fundamentals'*mon_spd;
    Loog = logical(out(:,7));
    ecc = sqrt((NT.sum.exptParams.rf_x/10).^2 + (NT.sum.exptParams.rf_y/10).^2);
    data(cellcounter).M = M;
    data(cellcounter).data = out;
    data(cellcounter).bkgndrgb = NT.sum.exptParams.bkgndrgb;
    data(cellcounter).ecc = ecc;
    data(cellcounter).monspd = mon_spd;
    fn = char(fnames{cellcounter});
    data(cellcounter).monk = fn(1);
end

% Now, computing cone weights as a function of macular pigment density
% This takes a while to run
macpigpeaks = linspace(-.2,.5,10);
for i = 1:length(data)
    cc_original = data(i).data(:,[2:4]).*repmat(data(i).data(:,5),1,3);
    Loog = data(i).data(:,end);
    allconeweights = [];
    for j = 1:length(macpigpeaks)
        macpigtransmittance = 1./(10.^(macpig.*macpigpeaks(j)));
        tmpfunds = fundnopig.*repmat(macpigtransmittance',1,3);
        tmpfunds = tmpfunds./repmat(max(tmpfunds),81,1);
        M = tmpfunds'*mon_spd;
        scaled = ConvertConeContrastBasis(data(i).M, M, data(i).bkgndrgb, cc_original);
        [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(scaled, Loog);
        coneweights = (xformmat*planeparams)';
        if (sign(coneweights(2)) == -1)
            coneweights = -coneweights;
        end
        allconeweights = [allconeweights; coneweights];
    end
    data(i).coneweights = bsxfun(@rdivide, allconeweights, sum(abs(allconeweights),2));
end

predictedmacpig = [];
for i = 1:length(data)
    predictedmacpig(i) = interp1(data(i).coneweights(:,3),macpigpeaks,0,'linear')
end
iskali = logical([data.monk]=='K');
ecc = [data.ecc];

figure('Units','inches','PaperPosition',[0 0 8.5 11],'Position',[0 0 8.5 11],'PaperOrientation','Portrait','DefaultAxesUnits','inches');
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultTextFontName','Helvetica','DefaultTextFontSize',12);
colordef white;
axes('position',[2 2 4 4]); hold on;

h = [];
h(1) = plot(ecc(iskali), predictedmacpig(iskali),'ko','MarkerFaceColor','black','MarkerSize',7);
h(2) = plot(ecc(~iskali), predictedmacpig(~iskali),'ro','MarkerFaceColor','red','MarkerSize',7);

xlabel('Eccentricity (deg)','FontSize',12);
ylabel('Predicted peak macular pigment density (460 nm)','FontSize',12);
L = ~isnan(predictedmacpig) & predictedmacpig > min(macpigpeaks) & predictedmacpig < max(macpigpeaks);
cordata = [[data.ecc]' predictedmacpig'];
cordata = cordata(L,:);
[r,p] = corr(cordata,'type','Spearman');
title(['r = ',num2str(r(1,2),2),' p = ',num2str(p(1,2),2)]);

[b,bint,r,rint,stats] = regress(predictedmacpig(L)', [ecc(L)' ones(sum(L),1)]);
plot([min(ecc) max(ecc)],b(1)*[min(ecc) max(ecc)]+b(2),'k-','LineWidth',2,'LineStyle',':');
%plot([2 7],[1 1]*macpigmult2deg*max(macpig),':')
%plot([2 7],[1 1]*macpigmult10deg*max(macpig),':')
%plot([2 7],[0 0],':')
set(gca,'YLim',[-.2 .4]);
legend(h,{'Monkey K','Monkey S'},'location','SouthWest');

[b,bint,r,rint,stats] = regress(predictedmacpig(L)', [iskali(L)' ecc(L)' ones(sum(L),1)]); % effect of ecc with monkey covariate
bint
[h,p] = ttest2(predictedmacpig(iskali), predictedmacpig(~iskali))

%%
% Section 5.1 
% S-cone weight as a function of retinal eccentricity and monkey (not
% converting to nay kind of pigment density)

% First gathering all the data
load ('T_cones_myss10.mat');

data = [];
if (ispc)
    [fnames, spikeIdx] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\NT\NTlum2.txt');
else
    [fnames, spikeIdx] = fnamesFromTxt2('/Volumes/NO BACKUP/NexFiles/nexfilelists/Greg/NT/NTlum2.txt');
end

for cellcounter = 1:size(fnames,1)
    filename = findfile(char(fnames{cellcounter}));
    paradigmID = getparadigmID(filename);
    if paradigmID == 103
        NT = nex2stro(filename);
    end
    out = NTpreprocess(NT,0,Inf);
    fundamentals = NT.sum.exptParams.fundamentals;
    fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
    mon_spd = NT.sum.exptParams.mon_spd;
    mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
    mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
    Loog = logical(out(:,7));
    ecc = sqrt((NT.sum.exptParams.rf_x/10).^2 + (NT.sum.exptParams.rf_y/10).^2);  
    
    % Converting to SS 10
    scaled = ConvertConeContrastBasis(fundamentals'*mon_spd, T_cones_ss10'*mon_spd, NT.sum.exptParams.bkgndrgb, out(:,[2:4]).*repmat(out(:,5),1,3));
    Loog = out(:,end);
    [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(scaled, Loog);
    coneweights = (xformmat*planeparams)';
    if (sign(coneweights(2)) == -1)
        coneweights = -coneweights;
    end
    data(cellcounter).coneweights = coneweights./sum(abs(coneweights));
    data(cellcounter).ecc = ecc;
    data(cellcounter).monspd = mon_spd;
    fn = char(fnames{cellcounter});
    data(cellcounter).monk = fn(1);
end
iskali = logical([data.monk]=='K');
ecc = [data.ecc];
coneweights = reshape([data.coneweights],3,length(data))';


figure('Units','inches','PaperPosition',[0 0 8.5 11],'Position',[0 0 8.5 11],'PaperOrientation','Portrait','DefaultAxesUnits','inches');
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultTextFontName','Helvetica','DefaultTextFontSize',12);
colordef white;
axes('position',[2 2 3 3]); hold on;

h = [];
h(1) = plot(ecc(iskali), coneweights(iskali,3),'ko','MarkerFaceColor','black','MarkerSize',7); % Black is Kali
h(2) = plot(ecc(~iskali), coneweights(~iskali,3),'ro','MarkerFaceColor','red','MarkerSize',7); % Red is Sedna

[b,bint,r,rint,stats] = regress(coneweights(:,3), [iskali' ecc' ones(length(data),1)]); % effect of ecc with monkey covariate
stats
p = stats(end)
ylabel('Normalized S-cone weight','FontSize',14);
xlabel('Eccentrcity (DVA)','FontSize',14);
plot([2 7],b(1) + b(3) + b(2).*[2 7],'k-');
plot([2 7],b(3) + b(2).*[2 7],'r-');

% Including an interaction term
[b,bint,r,rint,stats] = regress(coneweights(:,3), [iskali' ecc' ecc'.*iskali' ones(length(data),1)]); % effect of ecc with monkey covariate
legend(h,{'Monkey K','Monkey S'});
%%
% Section 6
% Example Ellipsoidal Neurothresh neuron
% Based on section 3

PLOTSURF = 1;
PLOTDATA = 1;
PLOTAXIS = 1;
filename = 'S073010002.nex';  % Pancolor Good example
%filename = 'S031610003.nex';  % Good one.
filename = 'S120310003.nex';

XLIM = .15; YLIM = .7;
AXESWIDTH = 2;
AXESMARGIN = .5;

load den_lens_ss;
load T_cones_myss10;

lenstransmittance = 1./(10.^(den_lens_ss));
fundnopig = T_cones_ss10./repmat(lenstransmittance,1,3);
fundnopig = fundnopig./repmat(max(fundnopig),size(fundnopig,1),1);
lensat400nm = [0 1 2];

% Getting data ready
stro = nex2stro(findfile(filename));
lmsidxs = [find(strcmp(stro.sum.trialFields(1,:),'lcont'))...
    find(strcmp(stro.sum.trialFields(1,:),'mcont'))...
    find(strcmp(stro.sum.trialFields(1,:),'scont'))];
% -------------------------------
% Converting cone contrasts in nex file to 10 deg fundmentals.
% Must go through excitations first - can't just transform contrasts.
% -------------------------------
fundamentals = stro.sum.exptParams.fundamentals;
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
mon_spd = stro.sum.exptParams.mon_spd;
mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
Moriginal = fundamentals'*mon_spd;
bkgndrgb = stro.sum.exptParams.bkgndrgb;
out = NTpreprocess(stro,0,Inf);
cc_original = out(:,[2 3 4]).*repmat(out(:,5),1,3);
Loog = logical(out(:,7));
wls = [380:5:780];

figure('Units','inches','PaperPosition',[0 0 8.5 11],'Position',[0 0 8.5 11],'PaperOrientation','Portrait','DefaultAxesUnits','inches');
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultTextFontName','Helvetica','DefaultTextFontSize',12);
colordef white;

for i = 1:length(lensat400nm)     
    % making new fundamentals
    lenstransmittance = 1./(10.^(den_lens_ss.*lensat400nm(i)));
    tmpfund = fundnopig.*repmat(lenstransmittance,1,3);
    tmpfund = tmpfund./repmat(max(tmpfund),81,1);
    M = tmpfund'*mon_spd;
    scaled = ConvertConeContrastBasis(Moriginal,M, bkgndrgb, cc_original);
    
    D = [scaled.^2,...
        2.*scaled(:,1) .* scaled(:,2),...
        2.*scaled(:,1) .* scaled(:,3),...
        2.*scaled(:,2) .* scaled(:,3)];
    lssoln = (D' * D) \(D' * ones(size(scaled,1),1));
    options = optimset('MaxFunEvals',50000,'MaxIter',50000,'TolFun',10^-6,'TolX',10^-6,'Display','off');    
    [quadparams, quadSSE, exitflag] = fminsearch(@(x) surfacefiterr4(scaled, x, Loog),lssoln,options);
    
    A = [quadparams(1) quadparams(4) quadparams(5);...
        quadparams(4) quadparams(2) quadparams(6);...
        quadparams(5) quadparams(6) quadparams(3)];
    cond(A)
    [evecs, evals] = eig(A);
    [evals,idxs] = sort(diag(evals));
    evecs = evecs(:,idxs);
    %coneweights = evecs(:,1)./sum(abs(evecs(:,1)))  % Small eigenvalue = long axis
    coneweights = evecs(:,1);
    if (coneweights(1) < 0)  % Positive L-cone weights
        coneweights = -coneweights;
    end

    % Setting x axis to the linear combination of l,m coneweights
    % that makes the maximum projection of long axis visible
    % y is orthogonal to x.
    lmcoefs = mkbasis(coneweights([1 2]));
    if (max(lmcoefs) ~= max(abs(lmcoefs)))
        lmcoefs = -lmcoefs;
    end
    % testing: looking at L+M vs S
    % lmcoefs = [1./sqrt(2); 1./sqrt(2)]
    
    xmat = [lmcoefs' 0; [lmcoefs(2) -lmcoefs(1) 0]; 0 0 1];
    xyscaled = scaled*xmat';
    B = xmat*A*xmat';
    coneweights = xmat*coneweights;
     
    axes('position',[1+(i-1)*(AXESWIDTH+AXESMARGIN) 3 AXESWIDTH AXESWIDTH])
    title(['Lens density @ 400 nm: ',num2str(lensat400nm(i))]);
    hold on;
    axis square;
    % Plotting the data points
    if (PLOTDATA)
        h = plot3(xyscaled(~Loog,1),xyscaled(~Loog,2),xyscaled(~Loog,3),'ko');
        set(h,'Markersize',3,'Markerfacecolor','black');
        h = plot3(-xyscaled(~Loog,1),-xyscaled(~Loog,2),-xyscaled(~Loog,3),'ko');
        set(h,'Markersize',3,'Markerfacecolor','black');
        if (lmcoefs(2) < 0)
            signsym = '-';
        else
            signsym = '+';
        end
        xlabel([num2str(lmcoefs(1),2),'L',signsym,num2str(abs(lmcoefs(2)),2),'M']);
        set(get(gca,'XLabel'),'Color',[0 0 0]);
        if (i == 1)
            zlabel('S'); set(get(gca,'YLabel'),'Color',[0 0 0]);
        else
            set(gca,'YTick',[]);
        end
    end    
     % Adjusting the axes
     LIM = max(abs(scaled(:,3)));
     set(gca,'XLim',LIM*[-1 1],'Ylim',LIM*[-1 1])
 
    % Plot surfaces
    if (PLOTSURF)
        [xx, yy, zz] = meshgrid(linspace(-LIM,LIM,50),...
            linspace(-LIM,XLIM,50),...
            linspace(-LIM,LIM,50));
        xformedxyz = [xx(:) yy(:) zz(:)];
        
        variables = [xformedxyz(:,1).^2 xformedxyz(:,2).^2 xformedxyz(:,3).^2 2*xformedxyz(:,1).*xformedxyz(:,2) 2*xformedxyz(:,1).*xformedxyz(:,3) 2*xformedxyz(:,2).*xformedxyz(:,3)];
        coefficients = [B(1,1) B(2,2) B(3,3) B(1,2) B(1,3) B(2,3)]';
        fr = variables*coefficients;
        fv = isosurface(xx,yy,zz,reshape(fr,size(xx)), 1);
        h = patch(fv);
        set(h,'FaceAlpha',0.5,'EdgeAlpha',0,'FaceColor','green','EdgeColor','none');
     %   set(h,'FaceVertexCData',repmat([0 .75 0],size(fv.vertices,1),1))
     %   set(h,'CDataMapping','direct');
     %   set(h,'FaceColor','interp','EdgeColor','none');
        lighting gouraud
        material metal
        camlight headlight 
    end

    if (PLOTAXIS)
        % Plotting the major axis
        plot3(coneweights(1)*LIM*[-1 1],coneweights(2)*LIM*[-1 1],coneweights(3)*LIM*[-1 1],'k-','LineWidth',4);
        text(1.3,.2,num2str(coneweights(3),4),'units','inches','fontsize',12);
    end
    set(gca,'View',[0 0]);
end
text(1,1,0,filename);
%if (PLOTSURF)
% Import to Illustrator and make the green parts transparent
%    set(gcf,'Renderer','openGL');
%    set(gcf,'InvertHardCopy','off')
%    print -dtiff -r600 -cmyk junk1
%else
    set(gcf,'Renderer','painters');  
%end

%%
% Section 7
% The long axis orientation of ellipsoidal isoresponse surfaces as a
% function of lens density.

load ('T_cones_myss10.mat');
wls = [380:5:780];
load den_lens_ss;
lenstransmittance = 1./(10.^(den_lens_ss));
fundnopig = T_cones_ss10./repmat(lenstransmittance,1,3);
fundnopig = fundnopig./repmat(max(fundnopig),size(fundnopig,1),1);

% First gathering all the data
data = [];
if (ispc)
    [fnames, spikeIdx] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\NT\pancolor2.txt');
else
    [fnames, spikeIdx] = fnamesFromTxt2('/Volumes/NO BACKUP/NexFiles/nexfilelists/Greg/NT/pancolor2.txt');
end
for cellcounter = 1:size(fnames,1)
    filename = findfile(char(fnames{cellcounter}));
    paradigmID = getparadigmID(filename);
    if paradigmID == 103
        NT = nex2stro(filename);
    end
    out = NTpreprocess(NT,0,Inf);
    fundamentals = NT.sum.exptParams.fundamentals;
    fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
    mon_spd = NT.sum.exptParams.mon_spd;
    mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
    mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
    M = fundamentals'*mon_spd;
    Loog = logical(out(:,7));
    ecc = sqrt((NT.sum.exptParams.rf_x/10).^2 + (NT.sum.exptParams.rf_y/10).^2);
    data(cellcounter).M = M;
    data(cellcounter).data = out;
    data(cellcounter).bkgndrgb = NT.sum.exptParams.bkgndrgb;
    data(cellcounter).ecc = ecc;
    data(cellcounter).monspd = mon_spd;
end

% Now, computing cone weights as a function of lens density
lensat400nm = linspace(0,2,20);
for i = 1:length(data)
    i
    cc_original = data(i).data(:,[2:4]).*repmat(data(i).data(:,5),1,3);
    Loog = data(i).data(:,end);
    allconeweights = [];
    for j = 1:length(lensat400nm)
        j
        lenstransmittance = 1./(10.^(den_lens_ss.*lensat400nm(j)));
        tmpfund = fundnopig.*repmat(lenstransmittance,1,3);
        tmpfund = tmpfund./repmat(max(tmpfund),81,1);
        M = tmpfund'*mon_spd;
        scaled = ConvertConeContrastBasis(data(i).M, M, data(i).bkgndrgb, cc_original);
        % [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(scaled, Loog);
        
        % trying an alternative means of fitting ellipsoids
        options = optimset('MaxFunEvals',50000,'MaxIter',50000,'TolFun',10^-6,'TolX',10^-6,'Display','off');
        D = [scaled.^2,...
            2.*scaled(:,1) .* scaled(:,2),...
            2.*scaled(:,1) .* scaled(:,3),...
            2.*scaled(:,2) .* scaled(:,3)];
        lssoln = (D' * D) \(D' * ones(size(scaled,1),1));
        
        [quadparams, quadSSE, exitflag] = fminsearch(@(x) surfacefiterr4(scaled, x, Loog),lssoln,options);
        A = [quadparams(1) quadparams(4) quadparams(5);... % adjusting quadparams
            quadparams(4) quadparams(2) quadparams(6);...
            quadparams(5) quadparams(6) quadparams(3)];
        
        [evecs, evals] = eig(A);
        [evals,idxs] = sort(diag(evals));
        evecs = evecs(:,idxs);
        coneweights = evecs(:,1);  % Small eigenvalue = long axis
        if (coneweights(3) < 0)  % Positive S-cone weights
            coneweights = -coneweights;
        end
         allconeweights = [allconeweights; coneweights'];
    end
    data(i).coneweights = allconeweights;
end

% Computing means and standard deviations
allconeweights = cat(3,data(:).coneweights);
mn = []; stdev = [];
for j = 1:length(lensat400nm)
    tmpconeweights = squeeze(allconeweights(j,:,:));
    normconeweights = tmpconeweights./repmat(sum(abs(tmpconeweights)),3,1);
    mn(j) = mean(normconeweights(3,:));
    stdev(j) = std(normconeweights(3,:));
end

% Doing the plotting
figure('Units','inches','PaperPosition',[0 0 8.5 11],'Position',[0 0 8.5 11],'PaperOrientation','Portrait','DefaultAxesUnits','inches');
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultTextFontName','Helvetica','DefaultTextFontSize',12);
colordef white;
axes('position',[2 2 4 4]); hold on;

tcrit = tinv(.025,length(mn)-1);
sem = stdev./sqrt(size(normconeweights,2));
CI_upper = tcrit*sem;
x = linspace(lensat400nm(1),lensat400nm(end),200); 
plot(x,interp1(lensat400nm,mn,x,'spline'),'k-','LineWidth',2);
plot(x,interp1(lensat400nm,mn+sem,x,'spline'),'k-','LineWidth',1);
plot(x,interp1(lensat400nm,mn-sem,x,'spline'),'k-','LineWidth',1);

xlabel('Lens pigment density at 400 nm','FontSize',14);
ylabel('Mean S-cone component','FontSize',14);

%set(gca,'XLim',[-.03 .4],'Ylim',[.8 .9]);
%plot([1 1]*lensat400nm(mn == max(mn)), [.8 max(mn)],'k-')

%%
% Section 8
% Distributions of S-cone weights among L+M cells using 2, 10, 10 w/o MP
% and new fundamentals. 

load ('T_cones_myss2');
load ('T_cones_myss10');
load ('T_cones_synthgh2');
synthfunds = T_cones_synthgh2;
wls = [380:5:780];

load ('den_mac_ss');
macpig = (den_mac_ss./max(den_mac_ss));
maxpigpeak = 0.15;
macpigtransmittance = 1./(10.^(macpig*maxpigpeak));
T_cones_synthgh5_mp = T_cones_synthgh2.*repmat(macpigtransmittance,1,3);
%peak10deg = 0.095; % from SS paper
%peak2deg = 0.35; % from SS paper

allfunds = cat(3,T_cones_ss2,T_cones_ss10,T_cones_synthgh2);

% First gathering all the data
data = [];
if (ispc)
    [fnames, spikeIdx] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\NT\NTlum2.txt');
else
    [fnames, spikeIdx] = fnamesFromTxt2('/Volumes/NO BACKUP/NexFiles/nexfilelists/Greg/NT/NTlum2.txt');
   
end
for cellcounter = 1:size(fnames,1)
    filename = findfile(char(fnames{cellcounter}));
    paradigmID = getparadigmID(filename);
    if paradigmID == 103
        NT = nex2stro(filename);
    end
    out = NTpreprocess(NT,0,Inf);
    fundamentals = NT.sum.exptParams.fundamentals;
    fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
    mon_spd = NT.sum.exptParams.mon_spd;
    mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
    mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
    originalM = fundamentals'*mon_spd;
    Loog = logical(out(:,7));
    bkgndrgb = NT.sum.exptParams.bkgndrgb;
    cc_original = out(:,[2:4]).*repmat(out(:,5),1,3);
    allconeweights = [];
    for j = 1:size(allfunds,3)
        M = squeeze(allfunds(:,:,j))'*mon_spd;
        scaled = ConvertConeContrastBasis(originalM, M, bkgndrgb, cc_original);
        [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(scaled, Loog);
        coneweights = (xformmat*planeparams)';
        if (sign(coneweights(2)) == -1)
            coneweights = -coneweights;
        end
        allconeweights = [allconeweights; coneweights./sum(abs(coneweights))];
    end
    data(cellcounter).coneweights = allconeweights;
end
Sconeweights = [];
for i = 1:size(data,2)
    Sconeweights(i,:) = data(i).coneweights(:,3)';
end
bins = linspace(min(Sconeweights(:)),max(Sconeweights(:)),15);
binwidth = bins(2)-bins(1);

% Plotting
figure('Units','inches','PaperPosition',[0 0 8.5 11],'Position',[0 0 8.5 11],'PaperOrientation','Portrait','DefaultAxesUnits','inches');
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultTextFontName','Helvetica','DefaultTextFontSize',12);
colordef white;

positions = [1 6 2 2; 1 2 2 2; 4 6 2 2; 4 2 2 2];
titles = {'S-cone weights (2 deg)',
          'S-cone weights (10 deg)',
          'S-cone weights (reduced lens, no macular pigment)',
          ['S-cone weights (reduced lens, ',num2str(maxpigpeak),' MP density)']};
for i = 1:4
    axes('position',positions(i,:)); hold on;
    [n,x] = hist(Sconeweights(:,i),bins);
    h = bar(x,n); set(h,'FaceColor','black');
    set(gca,'Xlim',[bins(1)-binwidth/2 bins(end)+binwidth/2]);
    set(gca,'YLim',[0 10]);
    title(titles{i});
    [h,p] = ttest(Sconeweights(:,i));
    [i p]
    xlabel('S-cone weight');
    ylabel('Count');
    h = plot(mean(Sconeweights(:,i)),9,'kv','Markersize',8);
    if (p < 0.01)
        set(h,'MarkerFaceColor','red','MarkerEdgeColor','red');
    end
end

%%
% Section 9 
% DTscot analyses
% Quantifying thresholds in cd/m^2 (photopic)
AXESWIDTH = 1.5;
if(ispc)
    filelistpath = 'N:\NexFiles\nexfilelists\Greg\DTEM';
else
   filelistpath = '/Volumes/NO BACKUP/NexFiles/nexfilelists/Greg/DTEM';
end

listnames = {'ZackDTscot.txt','GregDTscot.txt','LeahDTscot.txt','FreyaDTscot.txt','SednaDTscot.txt' };
subjectnames = {'Human Z','Human G','Human L','Monkey F','Monkey S'};
%listnames = {'GregDTscot.txt'};
%subjectnames = {'Human G'};

load('T_xyz1931')
vlambda = T_xyz1931(2,:);

data = [];
for i = 1:length(listnames)
    [fnames, spikenums] = fnamesFromTxt2([filelistpath,filesep,listnames{i}]);
    for j = 1:length(fnames)
        stro = nex2stro(findfile(char(fnames{j})));
        
        r_idx = strcmp(stro.sum.trialFields(1,:), 'stim_r');
        g_idx = strcmp(stro.sum.trialFields(1,:), 'stim_g');
        b_idx = strcmp(stro.sum.trialFields(1,:), 'stim_b');
        rgbs = stro.trial(:,r_idx|g_idx|b_idx);
        questmode_idx = strcmp(stro.sum.trialFields(1,:), 'quest_mode');
        thresholds = nan*ones(3,1);
        
        for gun = 1:3
            L = logical(rgbs(:,gun));
            if (any(L))
                all_modes = stro.trial(L, questmode_idx);
                subplot(3,1,gun);
                plot(all_modes,'b.-');
                thresholds(gun) = all_modes(end);
            end
        end
        data = [data; thresholds' i];
    end
end

% Creating predictions
load('wrattennd1.mat');
filtertransmittance = wratten_nd1(:,2);
load('Dell4blackbkgnd.mat');
cal = cals{end};
P_device = SplineSpd(cal.S_device, cal.P_device, [380 5 81]);
load('wrattennd1.mat');
NDtransmittance = wratten_nd1(:,2);
new_mon_spd = P_device .* repmat(NDtransmittance .^6, 1, 3);

gunlums = 680*vlambda*new_mon_spd;
load T_rods;

% wls = [398 420 440 459 479 501 520 539 560 579 599 619 640 659 681 700 720];
% s = [-.697 -.436 -.314 -.146 -.017 0 -.151 -.356 -.749 -1.22 -1.755 -2.312 -3.093 -3.743 -4.503 -5.147 -5.657];
% rodactionspectrum = interp1(wls,s,[380:5:780],'spline','extrap')';
% % extrapolating the rod action spectrum in two parts to avoid the spline
% % coming up again at long wavelengths
% rodactionspectrum1 = interp1(wls,s,[380:5:780],'linear','extrap')';
% L = [380:5:780] > 720;
% rodactionspectrum(L) = rodactionspectrum1(L);
% %plot([380:5:780],10.^rodactionspectrum,'g.') % Sanity check

load den_lens_ss;
lens = den_lens_ss;
lenstransmittance = 1./(10.^(den_lens_ss))';
absorptance = T_rods./lenstransmittance; % starting with CIE data instead of Baylor data

lensdensityat400 = 1;
lenstransmittance = 1./(10.^(lens*(lensdensityat400./lens(5))))';
% opticaldensity = .35; % .4 is the estimate from Baylor et al. (1984)  0.2 is Rushton (1972)
% absorptance = 1-10.^(-(10.^rodactionspectrum).*opticaldensity);
% absorptance = absorptance./max(absorptance);
fund = absorptance.*lenstransmittance;
fund = fund./max(fund);
humanpred = (1./(T_rods * new_mon_spd)).*gunlums; % Predicting threshold, not sensitivity
monkeypred = (1./(fund' * new_mon_spd)).*gunlums;

% Doing the plotting
figure('Units','inches','PaperPosition',[0 0 8.5 11],'Position',[0 0 8.5 11],'PaperOrientation','Portrait','DefaultAxesUnits','inches');
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultTextFontName','Helvetica','DefaultTextFontSize',12);
colordef white;

% % Plotting the human scotopic vlambda and the synthetic monkey version
% axes('position',[2 9 AXESWIDTH AXESWIDTH]); hold on;
% plot([380:5:780],vlambda,'k-','LineWidth',1);
% plot([380:5:780],T_rods,'b-','LineWidth',2);
% plot([380:5:780],fund,'m-.','LineWidth',2);
% set(gca,'XLim',[380 780]);
% ylabel('Relative Sensitivity'); xlabel('Wavelength (nm)');
% legend({'V_\lambda','V''_\lambda'});

for i = 1:length(listnames)
    L = data(:,4) == i;
    r = log10(data(L,1)*gunlums(1));
    g = log10(data(L,2)*gunlums(2));
    b = log10(data(L,3)*gunlums(3));
    nr = sum(~isnan(r));
    ng = sum(~isnan(g));
    nb = sum(~isnan(b));
    
    % Plotting data on the theoretical predictions
    axes('position',[1+2.4*mod(i-1,3) 4*(i<=3)+2.5 AXESWIDTH AXESWIDTH]); hold on;
    if (i == 1)
        ylabel('Log_1_0 threshold (cd m^-^2)');        
    else
        %set(gca,'YTick',[]); 
    end
    
    % Axes #1: Plotting means and SSE
    % plot([r,g,b]','ko','MarkerFaceColor','black','MarkerSize',1);
    plot([1 1; 2 2 ; 3 3]',[-1 1]'*([nanstd(r) nanstd(g) nanstd(b)]./sqrt([nr ng nb]))+repmat([nanmean(r),nanmean(g),nanmean(b)],2,1),'k-','LineWidth',2)
    plot([nanmean(r),nanmean(g),nanmean(b)]','ko','MarkerSize',6,'MarkerFaceColor','black')
    set(gca,'XLim',[.5 3.5],'XTickLabel',{'Red','Green','Blue'});
  
    % Plotting theoretical predictions
    humanmeancorrection = mean(nanmean([r g b])-log10(humanpred));
    plot([1 2 3],log10(humanpred)+humanmeancorrection,'b-','LineWidth',2);
    monkeymeancorrection = mean(nanmean([r g b])-log10(monkeypred));
  %  plot([1 2 3],log10(monkeypred)+monkeymeancorrection,'m-.','LineWidth',2);
    ylims = get(gca,'Ylim');
    text(1,(ylims(2)-ylims(1))*.1+ylims(1),subjectnames{i},'FontSize',7);
    
    % ---------------------------
    % Axes #2: Plotting residuals
    axes('position',[1+2.4*mod(i-1,3) 4*(i<=3)+.5 AXESWIDTH AXESWIDTH]); hold on;
   % if (i == 1)
        ylabel('Residual');
   % else
   %     set(gca,'YTick',[]);
   % end
    set(gca,'Xlim',[.5 3.5]);
    resid = [r,g,b]-(repmat(log10(humanpred),length(r),1)+humanmeancorrection);
    
  %  plot(resid','k.');
    plot([1 1; 2 2 ; 3 3]',[-1 1]'*(nanstd(resid)./sqrt([nr ng nb]))+repmat(nanmean(resid),2,1),'k-','LineWidth',2)
    plot(nanmean(resid),'ko','MarkerSize',6,'MarkerFaceColor','black');
    plot([1 2 3],[0 0 0],'b-','LineWidth',2);
      
    %Looking at residuals from a synthetic "monkey" prediction
  %  plot([1 2 3],((log10(monkeypred)+monkeymeancorrection)-(log10(humanpred)+humanmeancorrection)),'m-.','linewidth',2)
    set(gca,'Ylim',[-.2 .2],'Xlim',[.5 3.5],'XTickLabel',{'Red','Green','Blue'});
end


%%
% Section 9.05
% DTscot analyses. Thresholds in normalized r and g units
% Need to run cell above first.
h = [];
norm_r = data(:,1)./sum(data(:,[1 2 3]),2);
norm_g = data(:,2)./sum(data(:,[1 2 3]),2);
norm_b = data(:,3)./sum(data(:,[1 2 3]),2);
normdata = [norm_r, norm_g, norm_b];

% Doing the plotting
figure('Units','inches','PaperPosition',[0 0 8.5 11],'Position',[0 0 8.5 11],'PaperOrientation','Portrait','DefaultAxesUnits','inches');
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultTextFontName','Helvetica','DefaultTextFontSize',12);
colordef white;

% Plotting the human scotopic vlambda and the synthetic monkey version
markers = ['s','d','^','<','o','p'];
axes('position',[2 2 3 3]); hold on;
for i = 1:data(end,end)
    if (strncmp(subjectnames{i},'H',1))
        color = [0 0 1];
    else
        color = [1 1 0];
    end
    L = data(:,end) == i;
    h(i) = plot(norm_g(L),norm_b(L),'ko','Marker',markers(i),'MarkerFaceColor',color,'MarkerEdgeColor','black');
end
xlabel('G/(R+G+B)');
ylabel('B/(R+G+B)');

% Plotting theoretical predictions: for humans
theory = 1./(T_rods * new_mon_spd);
h(length(h)+1) = plot(theory(2)./sum(theory), theory(3)./sum(theory),'m*','MarkerSize',6);

% Plotting theoretical predictions: lowered lens density
lensdensityat400 = 1;
lenstransmittance = 1./(10.^(lens*(lensdensityat400./lens(5))));
wls = [398 420 440 459 479 501 520 539 560 579 599 619 640 659 681 700 720]; % I forget where these number came from...
s = [-.697 -.436 -.314 -.146 -.017 0 -.151 -.356 -.749 -1.22 -1.755 -2.312 -3.093 -3.743 -4.503 -5.147 -5.657];
opticaldensity = .3;
rodactionspectra = interp1(wls,s,[380:5:780],'linear','extrap');
absorptance = 1-10.^(-(10.^rodactionspectra').*opticaldensity);
fund = absorptance.*lenstransmittance;
theory = 1./(fund' * new_mon_spd);
h(length(h)+1) = plot(theory(2)./sum(theory), theory(3)./sum(theory),'r*','MarkerSize',6);

legend(h,cat(2,subjectnames,{'V''_\lambda','lens_4_0_0=1'}));
set(gca,'Xscale','log','YScale','log');
set(gca,'Xlim',[min(normdata(:,2)) max(normdata(:,2))],'Ylim',[min(normdata(:,3)) max(normdata(:,3))])


% Permutation test on normalized thresholds
tmpdata = normdata;
tmpdata = data(:,[1 2 3]);
humanidx = find(strncmp('Human',subjectnames,5));
Lhuman = ismember(data(:,end),humanidx);
subjectidx = data(:,end);
% Excluding Leah
%Lleah = data(:,end) == find(strcmp(subjectnames,'Human L'));
%sum(Lleah)
%tmpdata(Lleah,:) = [];
%Lhuman(Lleah) = [];
%subjectidx(Lleah) = [];
% Excluding Zack
%Lzack = data(:,end) == find(strcmp(subjectnames,'Human Z'));
%sum(Lzack)
%tmpdata(Lzack,:) = [];
%Lhuman(Lzack) = [];
%subjectidx(Lzack) = [];

% mns = log10([mean(tmpdata(Lhuman,:)), mean(tmpdata(~Lhuman,:))]);
% niter = 500;
% stats = zeros(niter,length(mns));
% for i = 1:niter
%     L = Lhuman(randperm(length(Lhuman)));
%     stats(i,:) = log10([mean(tmpdata(L,:)), mean(tmpdata(~L,:))]);
% end
% figure;
% for i = 1:3
%     subplot(3,1,i);
%     stat = (stats(:,i)-stats(:,i+3));
%     hist(stat); hold on;
%     teststat = mns(i)-mns(i+3);
%     plot(teststat,0,'r*');
%     title(['p = ',num2str(sum(stat>teststat)./niter)]);
% end
% xlabel('human-monkey');

%squared distance in the g,b plane
% stat = (stats(:,2)-stats(:,5)).^2+(stats(:,3)-stats(:,6)).^2;
% figure; axes; hold on;
% hist(stat);
% teststat = (mns(2)-mns(5)).^2+(mns(3)-mns(6)).^2;
% plot(teststat,0,'r*');
% title(['p = ',num2str(sum(stat>teststat)./niter)]);

%[h,p] = ttest2(log10(data(Lhuman,3)./data(Lhuman,2)), log10(data(~Lhuman,3)./data(~Lhuman,2)))
%anovan(tmpdata(:,2),{subjectidx,Lhuman},'nested',[0 1; 0 0],'random',1);
anovan(log10(tmpdata(:,3)./tmpdata(:,2)),{subjectidx,Lhuman},'nested',[0 1; 0 0],'random',1);

%%
% Section 9.1
% DTscot analyses. Thresholds in scotopic trolands.
%

if(ispc)
    filelistpath = 'N:\NexFiles\nexfilelists\Greg\DTEM';
else
   filelistpath = '/Volumes/NO BACKUP/NexFiles/nexfilelists/Greg/DTEM';
end
listnames = {'GregDTscot.txt','ZackDTscot.txt','LeahDTscot.txt','SednaDTscot.txt','FreyaDTscot.txt','ApolloDTscot.txt'};
subjectnames = {'Human G','Human Z','Human L','Monkey S','Monkey F','Monkey A'};
%listnames = {'ZackDTscot.txt','ApolloDTscot.txt'};
%subjectnames = {'Human Z', 'Monkey A'};

data = [];
for i = 1:length(listnames)
    [fnames, spikenums] = fnamesFromTxt2([filelistpath,filesep,listnames{i}]);
    for j = 1:length(fnames)
        stro = nex2stro(findfile(char(fnames{j})));
        
        r_idx = strcmp(stro.sum.trialFields(1,:), 'stim_r');
        g_idx = strcmp(stro.sum.trialFields(1,:), 'stim_g');
        b_idx = strcmp(stro.sum.trialFields(1,:), 'stim_b');
        rgbs = stro.trial(:,r_idx|g_idx|b_idx);
        questmode_idx = strcmp(stro.sum.trialFields(1,:), 'quest_mode');
        thresholds = nan*ones(3,1);
        
        for gun = 1:3
            L = logical(rgbs(:,gun));
            if (any(L))
                all_modes = stro.trial(L, questmode_idx);
                subplot(3,1,gun);
                plot(all_modes,'b.-');
                thresholds(gun) = all_modes(end);
            end
        end
        data = [data; thresholds' i];
    end
end

% Creating predictions
load('Dell4blackbkgnd.mat');
cal = cals{end};
P_device = SplineSpd([380 2 201], cal.P_device, [380 5 81]);
% P_device has to be on a 2 nm lattice otherwise scotopic luminance levels
% are wrong.
load('wrattennd1.mat');
NDtransmittance = wratten_nd1(:,2);
new_mon_spd = P_device .* repmat(NDtransmittance .^6, 1, 3);

load T_rods;
% % Rod absorption spectrum based on Baylor electrophysiology predicts too
% % much blue sensitivity. Switching to a CIE V'(lambda) based observer
% wls = [398 420 440 459 479 501 520 539 560 579 599 619 640 659 681 700 720];
% s = [-.697 -.436 -.314 -.146 -.017 0 -.151 -.356 -.749 -1.22 -1.755 -2.312 -3.093 -3.743 -4.503 -5.147 -5.657];
% rodactionspectrum = interp1(wls,s,[380:5:780],'spline','extrap')';
% % extrapolating the rod action spectrum in two parts to avoid the spline
% % coming up again at long wavelengths
% rodactionspectrum1 = interp1(wls,s,[380:5:780],'linear','extrap')';
% L = [380:5:780] > 720;
% rodactionspectrum(L) = rodactionspectrum1(L);
% absorbance = 10.^rodactionspectrum; % linear units. absorptance = 1-10^(-OD*absorbance)
% %plot([380:5:780],10.^rodactionspectrum,'g.') % Sanity check
% lensdensityat400 = 1;
% load den_lens_ss;
% lens = den_lens_ss;
% lenstransmittance = 1./(10.^(lens*(lensdensityat400./lens(5))));
% opticaldensity = .35 % .4 is the estimate from Baylor et al.
% % 0.25 is the value form Werner. Setting the optical density to 0.2 give
% % predicted B/G of 0.76 - even further from data.
% absorptance = 1-10.^(-(10.^rodactionspectrum).*opticaldensity);  % GH 6/28/14
% fund = absorptance.*lenstransmittance;
% fund = fund./repmat(max(fund),81,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% New CIE-based scotopic observer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load den_lens_ss;
load den_mac_ss

%lensdensityat400 = 1.6;
%lenstransmittance = 1./(10.^(den_lens_ss*(lensdensityat400./lens(5))))';
%wls = [380 385 390 395 400 405 410 415 420 425 430 440 450 460 470 480 490 500 520 540 560 580 600 620 640 660 680 700];
%den_lens_nv_raw = [2.61 2.31 2.02 1.73 1.45 1.19 0.96 0.76 0.59 0.47 0.37 0.27 0.23 0.20 0.18 0.16 0.14 0.13 0.10 0.08 0.06 0.04 0.03 0.02 0.01 0 0 0]*1.16;
%den_lens_nv = interp1(wls,den_lens_nv_raw,[380:5:780],'linear',0);
lenstransmittance = 1./(10.^den_lens_ss)';

absorptance = T_rods./lenstransmittance; % starting with CIE data instead of Baylor data
%macpigtransmittance = 1./(10.^(den_mac_ss./max(den_mac_ss)*0.055))'; % might need to add this back?
% Adding just a little bit of macular pigment to see what effect this has
%absorptance = absorptance./macpigtransmittance;
lensdensityat400 = 1;
lenstransmittance = 1./(10.^(den_lens_ss*(lensdensityat400./den_lens_ss(5))))';
fund = absorptance.*lenstransmittance;
fund = fund./max(fund);

% Doing the plotting
figure('Units','inches','PaperPosition',[0 0 8.5 11],'Position',[0 0 8.5 11],'PaperOrientation','Portrait','DefaultAxesUnits','inches');
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultTextFontName','Helvetica','DefaultTextFontSize',12);
colordef white;

% Plotting the human scotopic vlambda and the synthetic monkey version
axes('position',[3.2 8 2 2]); hold on;
plot([380:5:780],T_rods,'b-','LineWidth',2);
plot([380:5:780],fund,'m-.','LineWidth',2);
set(gca,'XLim',[380 780]);
ylabel('Rel. Sens.'); xlabel('Wavelength (nm)');

scotlum = 1700*T_rods*new_mon_spd*5/2;% The factor of 5/2 comes from the fact that V'(lambda) is measured on a 5 nm lattice and 
% measurements from the PR705 are on a 2 nm lattice. I don't entirely get
% this but it works (see Section C).
scotlumatthresh = repmat(scotlum,size(data,1),1).*data(:,[1:3]);
%pupildiam = 5% % an esitmate
%pupilarea = pi*(pupildiam/2)^2; % mm^2
pupilarea = 50; % from W&S
scottrolatthresh = scotlumatthresh*pupilarea;

%scotlum = 1700*fund*new_mon_spd*5/2;
%scotmonklumatthresh = repmat(scotlum,size(data,1),1).*data(:,[1:3]);
%pupilarea = 50; % from W&S
%scotmonktrolatthresh = scotlumatthresh*pupilarea;

for i = 1:length(listnames)
    L = data(:,4) == i;
    r = log10(1./scotlumatthresh(L,1));
    g = log10(1./scotlumatthresh(L,2));
    b = log10(1./scotlumatthresh(L,3));
    nr = sum(~isnan(r));
    ng = sum(~isnan(g));
    nb = sum(~isnan(b));
    
    % Plotting data on the theoretical predictions
    axes('position',[1+2.4*(i-1) 5 2 2]); hold on;
    if (i == 1)
        ylabel('Log sensitivity (scotopic trolands^-^1)');
    else
  %      set(gca,'YTick',[]); 
    end
    
    % Plotting means and 2*SSE
    % plot([r,g,b]','ko','MarkerFaceColor','black','MarkerSize',1);
    plot([1 1; 2 2 ; 3 3]',[-1 1]'*([nanstd(r) std(g) std(b)]./sqrt([nr ng nb]))+repmat([nanmean(r),nanmean(g),nanmean(b)],2,1),'k-','LineWidth',2)
    plot([nanmean(r),nanmean(g),nanmean(b)]','ko','MarkerSize',4,'MarkerFaceColor','black')
    plot([nanmean(r),nanmean(g),nanmean(b)]','k-','LineWidth',2)
    set(gca,'XLim',[.5 3.5],'XTickLabel',{'Red','Green','Blue'});
  
    r = log10(1./scotmonklumatthresh(L,1));
    g = log10(1./scotmonklumatthresh(L,2));
    b = log10(1./scotmonklumatthresh(L,3));
    x = ([1 1; 2 2 ; 3 3]+.05)';
    plot(x,[-2 2]'*([nanstd(r) nanstd(g) nanstd(b)]./sqrt([nr ng nb]))+repmat([nanmean(r),nanmean(g),nanmean(b)],2,1),'b-','LineWidth',2)
    plot(x, [nanmean(r),nanmean(g),nanmean(b)]','bo','MarkerSize',4,'MarkerFaceColor','blue')
    plot(x, [nanmean(r),nanmean(g),nanmean(b)]','b-','LineWidth',2);
    ylims = get(gca,'Ylim');
    text(2.5,(ylims(2)-ylims(1))*.1+ylims(1),subjectnames{i},'FontSize',7);
end

tmpdata = scottrolatthresh;
subjectidx = data(:,end);
humanidx = find(strncmp('Human',subjectnames,5));
Lhuman = ismember(data(:,end),humanidx);

anovan(log10(tmpdata(:,3)./tmpdata(:,2)),{subjectidx,Lhuman},'nested',[0 1; 0 0],'random',1);

%%
% summary statistics
mn_scot_thresh = [];
for i = 1:max(subjectidx)
    mn_scot_thresh = [mn_scot_thresh; nanmean(log10(scottrolatthresh(subjectidx == i,[2 3])))]; % green and blue
end
Lhuman = ismember([1:max(subjectidx)],humanidx);
mean(mn_scot_thresh(Lhuman,2))./mean(mn_scot_thresh(Lhuman,1))
mean(mn_scot_thresh(~Lhuman,2))./mean(mn_scot_thresh(~Lhuman,1))

monk = [mean(mn_scot_thresh(~Lhuman,1)) std(mn_scot_thresh(~Lhuman,1)) mean(mn_scot_thresh(~Lhuman,2)) std(mn_scot_thresh(~Lhuman,2))]  % Monkeys
humans(3)./humans(1) % blue/green threshold ratios humans have lower B/G threshold ratios?
monk(3)./monk(1)
monkeydata = [mn_scot_thresh(Lhuman,1); mn_scot_thresh(Lhuman,2)];
humandata = [mn_scot_thresh(~Lhuman,1); mn_scot_thresh(~Lhuman,2)];

[mean(monkeydata) std(monkeydata)]
[mean(humandata) std(humandata)]



% % ---------------------------------
% % --------  FOR REVISION   --------
% % Getting scotopic B/G for each subject to address Reviewer #2's question
% % about whether variability in photopic B/G decreases when you correct for
% % the lens density implied by scotopic B/G
% 
% % First making a lookuptable between lens densities and scotopic B/G
% lenBGlut = []; % Lens to B/G look up table
% for i = [-1:.05:10]  % scaling lens absorbance function
%     lenstransmittance = 1./(10.^(lens*(i./lens(5))))';
%     % absorbtanceSpectra = 1 - 10.^(-OD * absorbanceSpectra)
%  %   absorbance = -log10(1-absorptance/max(absorptance)*.3)/.3;
%  %   absorptance_tmp = 1-10.^-(absorbance*.35);
%  %   fund = absorptance_tmp.*lenstransmittance;
%     fund = absorptance.*lenstransmittance;
%     scotlum_novel = 1700*fund*new_mon_spd;
%   %  pred = 10.^(log10(scotlum_novel(2)./scotlum_h(2))-log10(scotlum_novel(3)./scotlum_h(3))) % Long way of writing B/G
%   %  pred = (scotlum_h(3)*scotlum_novel(2))/(scotlum_h(2)*scotlum_novel(3)); % simplified version of above
%     pred = scotlum_novel(3)/scotlum_novel(2); 
%    
%     lenBGlut = [lenBGlut; i pred];
% end
% 
% for i = 1:max(subjectidx)
%     L = subjectidx == i;
%     l = data(L,2)./data(L,3); %  green/blue threshold ratio compared with theoretical blue/green sensitivity ratio
%     logse = std(log10(l))/length(l);
%     logmean = mean(log10(l));
%     predlensOD = interp1(lenBGlut(:,2), lenBGlut(:,1),10.^logmean,'linear');
%     SEbounds = interp1(lenBGlut(:,2), lenBGlut(:,1),[10.^(logmean+2*logse) 10.^(logmean-2*logse)],'linear');
%     [subjectnames{i},' ',num2str(geomean(l)),' ',num2str(predlensOD),' SE bounds: ',num2str(SEbounds)]
%     figure; hist(log10(l))
%    % Now trying to find the lqens density implied by the scotopic B/G
% end
% Interesting, the estimated lens density for the monkeys is uniformly > 1
% at 400 nm. I guess I knew that from Figure 6. All of these densities are
% too high. Use something other than Baylor data?
% Estimated lens densities at 400 nm.
% Human G  3.6793
% Human Z 1.0639
% Human L 1.9026
% Monkey S -0.022796
% Monkey F -0.11215
% Monkey A 0.99919 
%
% 25 year old should have lens OD of 1.6722 at 400 nm
% 40 year old should have lens OD of 1.8708 at 400 nm
% --------  END OF FOR REVISION   --------
% ---------------------------------


% n = length(unique(data(:,end)));
% figure;
% Loutlier = zeros(size(data,1),1);
% nstd = 2.5;
% for i = 1:n
%     L = logical(data(:,end) == i);
%     subplot(1,n,i); hold on
%     plot(log10(scotlumatthresh(L,2)),log10(scotlumatthresh(L,3)),'k.')
%     gm = mean(log10(scotlumatthresh(L,2)));
%     gs = std(log10(scotlumatthresh(L,2)));
%     bm = mean(log10(scotlumatthresh(L,3)));
%     bs = std(log10(scotlumatthresh(L,3)));
%     plot(gm,bm,'y*');
%     axis equal 
%     Loutlier(L) = log10(scotlumatthresh(L,2)) > gm+nstd*gs | log10(scotlumatthresh(L,2)) < gm-nstd*gs;
%     Loutlier(L) = log10(scotlumatthresh(L,3)) > bm+nstd*bs | log10(scotlumatthresh(L,3)) < bm-nstd*bs;
% end
% tmpdata = scotlumatthresh;
% anovan(log10(tmpdata(~Loutlier,3)./tmpdata(~Loutlier,2)),{subjectidx(~Loutlier),Lhuman(~Loutlier)},'nested',[0 1; 0 0],'random',1);
% 

%%
% The start of a new figure?
% Continued from above. Makes the point that monkeys 
% have a slightly lower lens density, consistent with '1' at 400 nm
% This code assumes that every data file has thresholds to blue and green.

% 
% lensdensityat400 = 1;
% lenstransmittance = 1./(10.^(den_lens_ss*(lensdensityat400./lens(5))))';
% fund = absorptance.*lenstransmittance;
% fund = fund./max(fund);
% 
% 
% % Doing the plotting
% figure('Units','inches','PaperPosition',[0 0 8.5 11],'Position',[0 0 8.5 11],'PaperOrientation','Portrait','DefaultAxesUnits','inches');
% set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
% set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
% set(gcf,'DefaultTextFontName','Helvetica','DefaultTextFontSize',12);
% colordef white;
% 
% % Plotting the human scotopic vlambda and the synthetic monkey version
% axes('position',[2.2 2 3 3]); hold on;
% 
% for i = 1:max(data(:,end))
%     L = data(:,end) == i;
%     errorbar(i,nanmean(log10(scotlumatthresh(L,3)./scotlumatthresh(L,2))),...
%         nanstd(log10(scotlumatthresh(L,3)./scotlumatthresh(L,2)))./sqrt(sum(L)),'k-');
%     h = plot(i,nanmean(log10(scotlumatthresh(L,3)./scotlumatthresh(L,2))),'ko','MarkerFaceColor','white','MarkerEdgeColor','black');
%     if (ismember(i,humanidx))
%         set(h,'marker','s','MarkerFaceColor','black');
%     end
%     %bar(i, mean(log10(scottrolatthresh(L,3)./scottrolatthresh(L,2))),.5)
% end
% set(gca,'Xtick',[1:max(data(:,end))],'Xlim',[0 max(data(:,end))]+0.5)
% set(gca,'XTickLabel',subjectnames);
% set(gca,'Ytick',log10([.8 .9 1 1.1]),'YtickLabel',[.8 .9 1 1.1]); 
% 
% % Finding two stimuli that are equal for human scotopic luminous efficiency
% % function and see what their ratio is for the monkey.
% scotlum_h = 1700*T_rods*new_mon_spd;
% scotlum_m = 1700*fund'*new_mon_spd;
% pred_m = log10(scotlum_m(2)./scotlum_h(2))-log10(scotlum_m(3)./scotlum_h(3));
% plot([1 max(data(:,end))],[0 0],'k:');
% plot([1 max(data(:,end))],[pred_m pred_m],'k:');
% 
% %[p,t,stats,terms] = anovan(log10(scottrolatthresh(:,3)./scottrolatthresh(:,2)),{data(:,end)<=3, data(:,end)},'nested',[0 0; 1 0],'random',2);
% %anovan(log10(scottrolatthresh(:,3)./scottrolatthresh(:,2)),{data(:,end)<=3, data(:,end)},'nested',[0 0; 1 0]);
% 
% ylabel({'Blue threshold / Green threshold'},'FontSize',14);
% xlabel('Subject','FontSize',14);
% 
% % Another panel showing a scatterplot?
% 
% axes('position',[2.2 7 3 3]); hold on;
% 
% for i = 1:max(data(:,end))
%     L = data(:,end) == i;
%     h = plot(scotlumatthresh(L,3), scotlumatthresh(L,2),'ko','MarkerFaceColor','white','MarkerEdgeColor','black');
%     if (i <=3 )
%         set(h,'marker','s','MarkerFaceColor','black');
%     end
%     %bar(i, mean(log10(scottrolatthresh(L,3)./scottrolatthresh(L,2))),.5)
% end
% set(gca,'Xscale','log','Yscale','log');
% lims = [min([scotlumatthresh(:,2);scotlumatthresh(:,3)]) max([scotlumatthresh(:,2);scotlumatthresh(:,3)])];
% set(gca,'Xlim',lims,'Ylim',lims);
% xlabel('Blue threshold (scotopic luminance)','FontSize',14);
% ylabel('Green threshold (scotopic luminance)','FontSize',14);

%%
% Section 9.2 Scotopic analysis and figure for J Neurophys resubmission

if(ispc)
    filelistpath = 'N:\NexFiles\nexfilelists\Greg\DTEM';
else
   filelistpath = '/Volumes/NO BACKUP/NexFiles/nexfilelists/Greg/DTEM';
end
listnames = {'GregDTscot.txt','ZackDTscot.txt','LeahDTscot.txt','SednaDTscot.txt','FreyaDTscot.txt','ApolloDTscot.txt'};
subjectnames = {'Human G','Human Z','Human L','Monkey S','Monkey F','Monkey A'};
humanidx = find(strncmp('Human',subjectnames,5));

data = [];
for i = 1:length(listnames)
    [fnames, spikenums] = fnamesFromTxt2([filelistpath,filesep,listnames{i}]);
    for j = 1:length(fnames)
        stro = nex2stro(findfile(char(fnames{j})));
        
        r_idx = strcmp(stro.sum.trialFields(1,:), 'stim_r');
        g_idx = strcmp(stro.sum.trialFields(1,:), 'stim_g');
        b_idx = strcmp(stro.sum.trialFields(1,:), 'stim_b');
        rgbs = stro.trial(:,r_idx|g_idx|b_idx);
        questmode_idx = strcmp(stro.sum.trialFields(1,:), 'quest_mode');
        thresholds = nan*ones(3,1);
        
        for gun = 1:3
            L = logical(rgbs(:,gun));
            if (any(L))
                all_modes = stro.trial(L, questmode_idx);
                subplot(3,1,gun);
                plot(all_modes,'b.-');
                thresholds(gun) = all_modes(end);
            end
        end
        data = [data; thresholds' i];
    end
end
subjectidx = data(:,end);

% Creating predictions
load('Dell4blackbkgnd.mat');
cal = cals{end};
P_device = SplineSpd([380 2 201], cal.P_device, [380 5 81]);
% P_device has to be on a 2 nm lattice otherwise scotopic luminance levels
% are wrong.
load('wrattennd1.mat');
NDtransmittance = wratten_nd1(:,2);
new_mon_spd = P_device .* repmat(NDtransmittance .^6, 1, 3);

load T_rods;
load den_lens_ss;
load den_mac_ss;
wls = [380 385 390 395 400 405 410 415 420 425 430 440 450 460 470 480 490 500 520 540 560 580 600 620 640 660 680 700];
den_lens_nv_raw = [2.61 2.31 2.02 1.73 1.45 1.19 0.96 0.76 0.59 0.47 0.37 0.27 0.23 0.20 0.18 0.16 0.14 0.13 0.10 0.08 0.06 0.04 0.03 0.02 0.01 0 0 0]*1.16;
den_lens_nv = interp1(wls,den_lens_nv_raw,[380:5:780],'linear',0);
lenstransmittance = 1./(10.^den_lens_nv);
absorptance = T_rods./lenstransmittance;
%macpigtransmittance = 1./(10.^(den_mac_ss./max(den_mac_ss)*0.055))'; % might need to add this back?
% Adding just a little bit of macular pigment to see what effect this has
%absorptance = absorptance./macpigtransmittance;
lensdensityat400 = 1;
lenstransmittance = 1./(10.^(den_lens_ss*(lensdensityat400./den_lens_ss(5))))';
fund = absorptance.*lenstransmittance;
fund = fund./max(fund);

scotlum = 1700*T_rods*new_mon_spd*5/2;% The factor of 5/2 comes from the fact that V'(lambda) is measured on a 5 nm lattice and 
% measurements from the PR705 are on a 2 nm lattice. I don't entirely get
% this but it works (see Section C).
scotlumatthresh = repmat(scotlum,size(data,1),1).*data(:,[1:3]); % scotopic cd/m^2
%pupildiam = 5% % an esitmate
%pupilarea = pi*(pupildiam/2)^2; % mm^2
pupilarea = 50; % from W&S
scottrolatthresh = scotlumatthresh*pupilarea;
% looks like 10^-4 scotopic trolands is what we're getting

% Doing the plotting
figure('Units','inches','PaperPosition',[0 0 8.5 11],'Position',[0 0 8.5 11],'PaperOrientation','Portrait','DefaultAxesUnits','inches');
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultTextFontName','Helvetica','DefaultTextFontSize',12);
colordef white;

% "barplot" with errorbars
axes('position',[2.2 2 3 3]); hold on;
for i = 1:max(data(:,end))
    L = data(:,end) == i;
    errorbar(i,nanmean(log10(scotlumatthresh(L,3)./scotlumatthresh(L,2))),...
        nanstd(log10(scotlumatthresh(L,3)./scotlumatthresh(L,2)))./sqrt(sum(L)),'k-');
    h = plot(i,nanmean(log10(scotlumatthresh(L,3)./scotlumatthresh(L,2))),'ko','MarkerFaceColor','white','MarkerEdgeColor','black');
    if (ismember(i,humanidx))
        set(h,'marker','s','MarkerFaceColor','black');
    end
end
set(gca,'Xtick',[1:max(data(:,end))],'Xlim',[0 max(data(:,end))]+0.5)
set(gca,'XTickLabel',subjectnames);
set(gca,'Ytick',log10([.7 .8 .9 1 1.1 1.2]),'YtickLabel',[.7 .8 .9 1 1.1 1.2]);
ylims = get(gca,'Ylim');
ylabel('Scotopic B/G');

% Making a lookup table between lens densities and scotopic B/G
lenBGlut = []; % Lens to B/G look up table
for i = [-2:.05:10]  % scaling lens absorbance function
    lenstransmittance = 1./(10.^(den_lens_ss*(i./den_lens_ss(5))))';
    fund = absorptance.*lenstransmittance;
    scotlum_novel = 1700*fund*new_mon_spd; % scotlum
    lenBGlut = [lenBGlut; i scotlum_novel(2)./scotlum(2) scotlum_novel(3)./scotlum(3)]; % Sensitivities
end
yt = interp1(lenBGlut(:,1),lenBGlut(:,2)./lenBGlut(:,3),linspace(0,4,5));
hax = axes('position',[2.2 2 3 3]);
set(hax,'YAxisLocation','right','Xtick',[],'Ylim',ylims,'Color','none');
set(hax,'YTick',log10(yt),'YTicklabel',linspace(0,4,5))
ylabel('Inferred lens OD @ 400 nm');

for i = 1:max(subjectidx)
    L = subjectidx == i;
    l = scotlumatthresh(L,3)./scotlumatthresh(L,2); % Thresholds
    logse = std(log10(l))/sqrt(sum(L));
    logmean = mean(log10(l));
    predlensOD = interp1(lenBGlut(:,2)./lenBGlut(:,3), lenBGlut(:,1),10.^logmean,'linear');
    SEbounds = interp1(lenBGlut(:,2)./lenBGlut(:,3), lenBGlut(:,1),[10.^(logmean+logse) 10.^(logmean-logse)],'linear');
    [subjectnames{i},' ',num2str(predlensOD),' SE bounds: ',num2str(SEbounds)]
end

%%
% Section 9.3
% Comparing lens density estimates from photopic and scotopic data
% for J Neurophys resubmission.

if(ispc)
    filelistpath = 'N:\NexFiles\nexfilelists\Greg\DTEM';
else
    filelistpath = '/Volumes/NO BACKUP/NexFiles/nexfilelists/Greg/DTEM';
end

%%%%%%%%%%%%%%%%%%
% Scotopic stuff
listnames = {'GregDTscot.txt','ZackDTscot.txt','LeahDTscot.txt','SednaDTscot.txt','FreyaDTscot.txt','ApolloDTscot.txt'};
subjectnames = {'Human G','Human Z','Human L','Monkey S','Monkey F','Monkey A'};
humanidx = find(strncmp('Human',subjectnames,5));

scotdata = [];
for i = 1:length(listnames)
    [fnames, spikenums] = fnamesFromTxt2([filelistpath,filesep,listnames{i}]);
    for j = 1:length(fnames)
        stro = nex2stro(findfile(char(fnames{j})));
        
        r_idx = strcmp(stro.sum.trialFields(1,:), 'stim_r');
        g_idx = strcmp(stro.sum.trialFields(1,:), 'stim_g');
        b_idx = strcmp(stro.sum.trialFields(1,:), 'stim_b');
        rgbs = stro.trial(:,r_idx|g_idx|b_idx);
        questmode_idx = strcmp(stro.sum.trialFields(1,:), 'quest_mode');
        thresholds = nan*ones(3,1);
        
        for gun = 1:3
            L = logical(rgbs(:,gun));
            if (any(L))
                all_modes = stro.trial(L, questmode_idx);
                thresholds(gun) = all_modes(end);  % gun intensity units
            end
        end
        scotdata = [scotdata; thresholds' i];
    end
end

%%%%%%%%%%%%%%%%%%
% Photopic stuff
listnames = {'GregMacPig.txt','ZackMacPig.txt','LeahMacPig.txt','SednaMacPig.txt','FreyaMacPig.txt','ApolloMacPig.txt'};
load Vlambda_star
Vlambda = Vlambda_star(:,2)';
photdata = [];
for i = 1:length(listnames)
    filenames = fnamesFromTxt2([filelistpath,filesep,char(listnames{i})]);
    for fileidx = 1:size(filenames,1)
        stro = nex2stro(findfile(char(filenames{fileidx,:})));
        [thresh, colorDirs, sfs] = DTquestUnpack(stro, 'mode'); close(gcf);
        normcolordirs = colorDirs./repmat(sqrt(sum(colorDirs.^2,2)),1,3);
        guns = reshape(stro.sum.exptParams.mon_spect,81,3);
        bkgndrgb = [stro.sum.exptParams.bkgnd_r stro.sum.exptParams.bkgnd_g stro.sum.exptParams.bkgnd_b];
        M = reshape(stro.sum.exptParams.m_mtx,3,3);
        bkgndlms = M*bkgndrgb';
        bkgndlum = bkgndrgb*guns'*Vlambda';
        
        threshrgb = inv(M)*(repmat(bkgndlms',2,1).*(1+repmat(thresh/100,[1 3]).*normcolordirs([1 2],:)))';
        deltagun = threshrgb-repmat(bkgndrgb',1,2);
        % Quantifying thresholds in 2 deg luminance contrast units
        thresholds = [nan; nan];
        for j = 1:size(deltagun,2)
            thresholds(j) = (guns*deltagun(:,j))'*Vlambda';
        %    thresholds(j) = thresholds(j)./bkgndlum;
        end
        ecc = stro.sum.exptParams.rf_x;
        photdata = [photdata; ecc/10 thresholds' i];
    end
end
% Trimming the ecentricities that are too close to the fovea
greg = photdata;
L = photdata(:,1) >= 5;
photdata(~L,:) = [];

%%%%%%%%%%%%%%%%%%
% Creating scotopic lookup tables/models
load('Dell4blackbkgnd.mat');
cal = cals{end};
P_device = SplineSpd([380 2 201], cal.P_device, [380 5 81]);
% P_device has to be on a 2 nm lattice otherwise scotopic luminance levels
% are wrong.
load('wrattennd1.mat');
NDtransmittance = wratten_nd1(:,2);
new_mon_spd = P_device .* repmat(NDtransmittance .^6, 1, 3);

load T_rods;
load den_lens_ss
%wls = [380 385 390 395 400 405 410 415 420 425 430 440 450 460 470 480 490 500 520 540 560 580 600 620 640 660 680 700];
%den_lens_nv_raw = [2.61 2.31 2.02 1.73 1.45 1.19 0.96 0.76 0.59 0.47 0.37 0.27 0.23 0.20 0.18 0.16 0.14 0.13 0.10 0.08 0.06 0.04 0.03 0.02 0.01 0 0 0]*1.16;
%den_lens_nv = interp1(wls,den_lens_nv_raw,[380:5:780],'linear',0);
lenstransmittance = 1./(10.^den_lens_ss');
rodabsorptance = T_rods./lenstransmittance; % removing lens density

% Making a lookup table between lens densities and scotopic B/G
lenscotBGlut = []; % Lens to B/G look up table
for i = [-2:.05:10]  % scaling lens absorbance function
    lenstransmittance = 1./(10.^(den_lens_ss*(i./den_lens_ss(5))))';
    fund = rodabsorptance.*lenstransmittance;
    scotlum_novel = 1700*fund*new_mon_spd; % scotlum
    lenscotBGlut = [lenscotBGlut; i scotlum_novel(2) scotlum_novel(3)]; % Sensitivities
end

% Getting mean and SE of lens density from scotopic data
scotsumstats = [];
for i = 1:length(listnames)
   L = scotdata(:,end) == i;
   logmnBG = mean(log10(scotdata(L,3)./scotdata(L,2)));
   logsemBG = std(log10(scotdata(L,3)./scotdata(L,2)))./sqrt(sum(L));
   mn = interp1(lenscotBGlut(:,2)./lenscotBGlut(:,3),[-2:.05:10],10.^logmnBG);
   bound1 = interp1(lenscotBGlut(:,2)./lenscotBGlut(:,3),[-2:.05:10],10.^(logmnBG+logsemBG));
   bound2 = interp1(lenscotBGlut(:,2)./lenscotBGlut(:,3),[-2:.05:10],10.^(logmnBG-logsemBG));
   scotsumstats = [scotsumstats; bound2 mn bound1];
end

%%%%%%%%%%%%%%%%%%
% Creating photopic lookup tables/models
load('Dell4BitsCal.mat');
cal = cals{8};  % 2011 calibration structure (small differences don't matter)
P_device = SplineSpd([380 4 101], cal.P_device, [380 5 81]);
glum = Vlambda*P_device(:,2);
blum = Vlambda*P_device(:,3);

load('T_cones_myss10')
load den_lens_ss;
load den_mac_ss;
lenstransmittance = 1./(10.^den_lens_ss);
coneabsorptance = T_cones_ss10./repmat(lenstransmittance,1,3);
coneabsorptance = coneabsorptance./repmat(max(coneabsorptance),size(coneabsorptance,1),1);

photsumstats = [];
% Make sure this calculation is done correctly
for i = 1:length(listnames)
    if (i == 2 || i == 3)
        LMratio = [10 1];
    else
        LMratio = [1 1];        
    end
    lenphotBGlut = [];
    for j = [-2:.05:10]  % scaling lens absorbance function
        lenstransmittance = 1./(10.^(den_lens_ss*(j./den_lens_ss(5))))';
        fund = coneabsorptance'.*repmat(lenstransmittance,3,1);
        fund = fund./repmat(max(fund'),size(fund,2),1)';
        M = fund*P_device;
        gpred = LMratio*M([1 2],2);
        bpred = LMratio*M([1 2],3);
        lenphotBGlut = [lenphotBGlut; j gpred./glum bpred./blum]; % Sensitivities
    end
    L = photdata(:,end) == i;
    logmnBG = mean(log10(photdata(L,3)./photdata(L,2)));
    logsemBG = std(log10(photdata(L,3)./photdata(L,2)))./sqrt(sum(L));
    mn = interp1(lenphotBGlut(:,2)./lenphotBGlut(:,3),[-2:.05:10],10.^logmnBG);
    bound1 = interp1(lenphotBGlut(:,2)./lenphotBGlut(:,3),[-2:.05:10],10.^(logmnBG+logsemBG));
    bound2 = interp1(lenphotBGlut(:,2)./lenphotBGlut(:,3),[-2:.05:10],10.^(logmnBG-logsemBG));
    photsumstats = [photsumstats; bound2 mn bound1];
end

figure; axes; hold on; axis square;
for i = 1:length(listnames)
    plot([photsumstats(i,1) photsumstats(i,3)],[scotsumstats(i,2) scotsumstats(i,2)],'k-');
    plot([photsumstats(i,2) photsumstats(i,2)],[scotsumstats(i,1) scotsumstats(i,3)],'k-'); 
    h = plot(photsumstats(i,2),scotsumstats(i,2),'ks','MarkerSize',10);
    if (ismember(i,humanidx))
        set(h,'MarkerFaceColor','black');
    else
        set(h,'MarkerFaceColor','white');
    end
end
set(gca,'Xlim',[-.5 5],'Xtick',[0 1 2 3 4 5],'Ylim',[-.5 5],'Ytick',[0 1 2 3 4 5]);
plot([0 5],[0 5],'k--');
xlabel('Lens OD (@400 nm) photopic');
ylabel('Lens OD (@400 nm) scotopic');

[r,p] = prettycorr([photsumstats(:,3)  scotsumstats(:,3)])
%%
% Section 10
% Comparing estimates of monkey photopic spectral sensitivity from the
% literature to sums of L and M cone fundamentals.
% Cone fundamentals
load T_cones_synthgh2;
load T_cones_myss10;
load T_cones_myss2;

% Data sets
%load schrier_1966;
load de_valois_1974;
load van_norren_1971;
%load sidley_1965;
load jacobs_1997.csv; % First column is from paper, Second is from erratum published 1999
h = []; g = [];

%schrier1966(:,2) = 10.^schrier1966(:,2);
van_norren_1971(:,2) = 10.^(van_norren_1971(:,2)-max(van_norren_1971(:,2)));
de_valois_1974(:,2) = 10.^de_valois_1974(:,2);
%sidley1965(:,2) = 10.^(sidley1965(:,2)-max(sidley1965(:,2)));
%sidley1965(:,2) = 10.^sidley1965(:,2);
%jacobs_1997(:,2) = 10.^jacobs_1997(:,2)./jacobs_1997(:,1);
jacobs_1997(:,2) = QuantaToEnergy(jacobs_1997(:,1),10.^jacobs_1997(:,3))
%jacobs_1997(:,2) = 10.^jacobs_1997(:,2)

figure; axes; hold on;
% 2 deg funds 2:1 (L:M)
load Vlambda_star
pred_ss2_2 = Vlambda_star(:,2);
h(1)=plot([380:5:780],pred_ss2_2,'r-','LineWidth',2);

% 10 deg funds 1:1 (L:M)
pred_ss10_1 = T_cones_ss10(:,1)+T_cones_ss10(:,2);
pred_ss10_1 = pred_ss10_1./max(pred_ss10_1);
h(2)=plot([380:5:780],pred_ss10_1,'c-','LineWidth',2);

% Lowered lens pigmentation
pred_synth2 = T_cones_synthgh2(:,1)+T_cones_synthgh2(:,2);
pred_synth2 = pred_synth2./max(pred_synth2);
h(3)=plot([380:5:780],pred_synth2,'k-','LineWidth',2);

% % Boettner lens pigmentation
% pred_synth1 = T_cones_synthgh1(:,1)+T_cones_synthgh1(:,2);
% pred_synth1 = pred_synth1./max(pred_synth1);
% plot([380:5:780],pred_synth1,'m-');

set(gca,'Ylim',[10^-3 2]);
set(gca,'YScale','log');

wls = [380:5:780];
wls_thresh = 500;
Lwhich_wls = logical(wls>=wls_thresh);

%mymodel = @(scale,x) (scale+interp1(wls(Lwhich_wls),log10(pred_ss10_1(Lwhich_wls)),log10(x),'linear'));
%mymodel = @(scale,x) (scale*interp1(wls(Lwhich_wls),pred_ss10_1(Lwhich_wls),x,'linear'));

L = de_valois_1974(:,1) >= wls_thresh;
%scale = nlinfit(de_valois_1974(L,1),de_valois_1974(L,2),mymodel,1);
%g(3)=plot(de_valois_1974(:,1),de_valois_1974(:,2)./scale,'c^','MarkerFaceColor','black');
mn(1) = mean(log10(de_valois_1974(L,2)));
mn(2) = mean(interp1(wls(Lwhich_wls),log10(pred_ss10_1(Lwhich_wls)),de_valois_1974(L,1),'linear'));
shift = mn(1)-mn(2);
g(1)=plot(de_valois_1974(:,1),de_valois_1974(:,2)./(10.^shift),'k^','MarkerFaceColor','black');

L = van_norren_1971(:,1) >= wls_thresh;
%scale = nlinfit(van_norren_1971(L,1),van_norren_1971(L,2),mymodel,0);
%g(1)=plot(van_norren_1971(:,1),van_norren_1971(:,2)./scale,'kv','MarkerFaceColor','black');
mn(1) = mean(log10(van_norren_1971(L,2)));
mn(2) = mean(interp1(wls(Lwhich_wls),log10(pred_ss10_1(Lwhich_wls)),van_norren_1971(L,1),'linear'));
shift = mn(1)-mn(2);
g(2)=plot(van_norren_1971(:,1),van_norren_1971(:,2)./(10.^shift),'kv','MarkerFaceColor','black');

L = jacobs_1997(:,1) >= wls_thresh;
%scale = nlinfit(jacobs_1997(L,1),jacobs_1997(L,2),mymodel,1);
%g(5)=plot(jacobs_1997(:,1),jacobs_1997(:,2)./scale,'k+','MarkerFaceColor','black');
mn(1) = mean(log10(jacobs_1997(L,2)));
mn(2) = mean(interp1(wls(Lwhich_wls),log10(pred_ss10_1(Lwhich_wls)),jacobs_1997(L,1),'linear'));
shift = mn(1)-mn(2);
g(3)=plot(jacobs_1997(:,1),jacobs_1997(:,2)./(10.^shift),'ks','MarkerFaceColor','black');

%L = sidley1965(:,1) >= wls_thresh;
%scale = nlinfit(sidley1965(L,1),sidley1965(L,2),mymodel,1);
%g(4)=plot(sidley1965(:,1),scale*sidley1965(:,2),'k+','MarkerFaceColor','black');

%L = schrier1966(:,1) >= wls_thresh;
%scale = nlinfit(schrier1966(L,1),schrier1966(L,2),mymodel,1);
%g(5)=plot(schrier1966(:,1),schrier1966(:,2)./scale,'ko','MarkerFaceColor','black');


xlabel('Wavelength (nm)','FontSize',14);
ylabel('Relative sensitivity','FontSize',14);
legend([h,g],{'V*(\lambda)(2�)','1L:M (10�)','1L:M (10�,lens=1@400nm)','DeValois 1974','Van Norren 1971','Jacobs 1997/1999'},'location',[.42 .25 .2 .2]);
set(gca,'Ylim',[.005 2])
%%
% Section 11: Comparing 15 Hz, 5� DTNT data across observers

CONVERTCC = 0; % Convert cone contrast to synthetic fundamentals
AXISWIDTH = 2;
XMARGIN = .5;
YMARGIN = 1;

LCONELIM = 30; % plotting limits
SCONELIM = 100; % plotting limits
% Getting the gamut of the ViewPixx for plotting
load Dell4BitsCal.mat
load T_cones_smj10
if (CONVERTCC)
    load T_cones_synthgh2
end

PERCGAMUT = 90; % For plotting the surface fits, how near to gamut edge to get (in each direction)
humanlist = {'GregDTNT15Hz.txt','ZackDTNT15Hz.txt','LeahDTNT15Hz.txt'};
monkeylist = {'NutDTNT15Hz.txt','FreyaDTNT15Hz.txt','SednaDTNT15Hz.txt'};
humansubjects = {'Greg','Zack','Leah'};
monkeysubjects = {'Nut','Freya','Sedna'};
if (ismac)
    filelistpath = '/Volumes/NO BACKUP/NexFiles/nexfilelists/Greg/DTNT';
else
    filelistpath = 'N:\NexFiles\nexfilelists\Greg\DTNT';
end

data = [];
for HUMANS = 0:1
    if (HUMANS)
        listnames = humanlist;
    else
        listnames = monkeylist;
    end
    for i = 1:length(listnames)
        filenames = fnamesFromTxt2([filelistpath,filesep,char(listnames{i})]);
        for fileidx = 1:size(filenames,1)
            if (strncmp(filenames{fileidx},'sf',2))
                continue;
            end
            stro = nex2stro(findfile(char(filenames{fileidx})));
            [thresholds, colordirs, sfs, QuestTrajectories] = DTquestUnpackGH(stro, 'mode');
            
            if (CONVERTCC)
                Morig = reshape(stro.sum.exptParams.m_mtx,3,3);
                spds = reshape(stro.sum.exptParams.mon_spect, length(stro.sum.exptParams.mon_spect)/3,3);
                bkgndrgb = [stro.sum.exptParams.bkgnd_r stro.sum.exptParams.bkgnd_g stro.sum.exptParams.bkgnd_b];
                Mnew = T_cones_synthgh2'*spds;
                cc_original = colordirs.*repmat(thresholds,1,size(colordirs,2));
                cc_new = ConvertConeContrastBasis(Morig, Mnew, bkgndrgb, cc_original);
                thresholds = sqrt(sum(cc_new.^2,2));
                colordirs = mkbasis(cc_new')';
            end
            
            threshs = thresholds';
            ncoldir = length(threshs);
            lmsmat = mkbasis(colordirs')';
            Loog = [];
            for j = 1:length(QuestTrajectories)
                Loog(j,1) = QuestTrajectories{j}(end) == max(QuestTrajectories{j});
            end
            data = [data; repmat(HUMANS,ncoldir,1) repmat(i,ncoldir,1) repmat(thresholds,1,3).*colordirs Loog];
        end
    end
end

fundamentals = T_cones_smj10;
fundWavelengthSpacing = S_cones_smj10;
calData = cals{end};
ptb.bkgndRGB = round(255 .* calData.bgColor); %bkgnd voltages discretized b/w 0&255
ptb.bkgndrgb = [calData.gammaTable(ptb.bkgndRGB(1)+1, 1), calData.gammaTable(ptb.bkgndRGB(2)+1, 2), ...
    calData.gammaTable(ptb.bkgndRGB(3)+1, 3)]; %add one to create and index b/w 1 & 256

ptb.monSpd = SplineSpd(calData.S_device, calData.P_device, fundWavelengthSpacing);
ptb.M = fundamentals * ptb.monSpd;
ptb.bkgndlms = ptb.M * ptb.bkgndrgb';

rgb_cube_verts = fullfact([2 2 2]) - 1;
lms_cube_verts = ptb.M * rgb_cube_verts';

cc_cube_verts = bsxfun(@rdivide, bsxfun(@minus, lms_cube_verts, ptb.bkgndlms), ptb.bkgndlms)';
% K = convhull(cc_cube_verts);
% T == NaN -> points outside the gamut
lims = max(abs(data(:,[3 4 5])));
npts = 40;
[xx yy zz] = meshgrid(linspace(-lims(1),lims(1),npts),...
    linspace(-lims(2),lims(2),npts),...
    linspace(-lims(3),lims(3),npts));
xformedxyz = [xx(:) yy(:) zz(:)];
T = tsearchn(cc_cube_verts*PERCGAMUT, delaunayn(cc_cube_verts*PERCGAMUT), xformedxyz);
xformedxyz(isnan(T),:) = nan;
T = tsearchn(cc_cube_verts*PERCGAMUT, delaunayn(cc_cube_verts*PERCGAMUT), -xformedxyz);
xformedxyz(isnan(T),:) = nan;
Loog = logical(data(:,6) == 1);

% Now doing the plotting
figure('Units','inches','PaperPosition',[0 0 8.5 11],'Position',[0 0 8.5 11],'PaperOrientation','Portrait','DefaultAxesUnits','inches');
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultTextFontName','Helvetica','DefaultTextFontSize',12);
colordef white;

for HUMANS = 0:1
    if (HUMANS)
        subjectlist = humansubjects;
    else
        subjectlist = monkeysubjects;
    end
    for i = 1:length(subjectlist)
        axes('position',[1+(i-1)*(AXISWIDTH+XMARGIN),3+HUMANS*(AXISWIDTH+YMARGIN),AXISWIDTH, AXISWIDTH]);
        hold on;
        L = data(:,1) == HUMANS & data(:,2) == i;
        for j = [-1, 1]
            h = plot3(j*data(L&~Loog,3),j*data(L&~Loog,4),j*data(L&~Loog,5),'ko');
            set(h,'MarkerFaceColor','black','MarkerSize',4,'MarkerEdgeColor','none');
            h = plot3([zeros(sum(L&Loog),1) j*data(L&Loog,3)]',[zeros(sum(L&Loog),1) j*data(L&Loog,4)]',[zeros(sum(L&Loog),1) j*data(L&Loog,5)]','k-');
            set(h(:),'Color',[.5 .5 .5]);
        end
        % Fitting data
        [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(data(L,[3 4 5]), Loog(L));
        planeparams=(planeparams'*xformmat'); % adjusting planeparams
        A = [quadparams(1) quadparams(4) quadparams(5);... % adjusting quadparams
            quadparams(4) quadparams(2) quadparams(6);...
            quadparams(5) quadparams(6) quadparams(3)];
        B = xformmat*A*xformmat';
        quadparams = [B(1,1) B(2,2) B(3,3) B(1,2) B(1,3) B(2,3)]';
        variables = [xformedxyz(:,1).^2 xformedxyz(:,2).^2 xformedxyz(:,3).^2 2*xformedxyz(:,1).*xformedxyz(:,2) 2*xformedxyz(:,1).*xformedxyz(:,3) 2*xformedxyz(:,2).*xformedxyz(:,3)];
        [evecs, evals] = eig(B);
        initparams = [2; reshape(evecs*sqrt(evals),9,1)];
        initparams(8:10) = []; % debugging
        
        options = optimset('MaxFunEvals',50000,'MaxIter',50000,'TolFun',10^-5,'Display','none');
        [fpar,fv] = fminsearch(@(x) colefiterr(x,data(L,[3 4 5]),Loog(L),0),initparams, options) 
        %  fr = sum(abs(xformedxyz *reshape(fpar(2:end),3,3)).^fpar(1),2);
        %  % commented out above line, debugging
         fr = sum(abs(xformedxyz *reshape(fpar(2:end),3,2)).^fpar(1),2); % debugging
        surfstruct = isosurface(xx,yy,zz,reshape(fr,size(xx)), 1);
        p = patch(surfstruct);
        set(p,'EdgeColor', 'none', 'FaceAlpha',.5,'FaceColor','green','Edgealpha',0);
        
%         % Plotting contour in the LM plane
%         tmpx = linspace(0,2*pi,100);
%         [th,phi,r] = cart2sph(cos(tmpx)', sin(tmpx)', zeros(length(tmpx),1));
%         wholesum = sum(abs([cos(phi(:)).*cos(th(:)) cos(phi(:)).*sin(th(:)) sin(phi(:))]*reshape(fpar(2:end),3,3)).^fpar(1),2);
%         predr = (1./wholesum).^(1/fpar(1));
%         contour =[predr.*cos(tmpx)'  predr.*sin(tmpx)'];
%         T = tsearchn(cc_cube_verts*PERCGAMUT, delaunayn(cc_cube_verts*PERCGAMUT), [contour zeros(size(contour,1),1)]);
%         contour(isnan(T),:) = nan;
%        plot(contour(:,1),contour(:,2),'LineWidth',3,'Color',[0 .3 0]);
        
        set(gca,'View',[0 90]);
        set(gca,'Xlim',[-LCONELIM LCONELIM],'Ylim',[-LCONELIM LCONELIM],'Zlim',[-SCONELIM SCONELIM]);
        set(gca,'Xtick',[-LCONELIM 0 LCONELIM],'Xticklabel',num2str([-LCONELIM 0 LCONELIM]'/100))
        xlabel('\DeltaL/L');
        set(gca,'Ytick',[-LCONELIM 0 LCONELIM],'Yticklabel',num2str([-LCONELIM 0 LCONELIM]'/100))
        ylabel('\DeltaM/M');
        set(gca,'Ztick',[-SCONELIM 0 SCONELIM],'Zticklabel',num2str([-SCONELIM 0 SCONELIM]'/100))
        zlabel('\DeltaS/S');
        textpos = [10, -28, 0];
         if (HUMANS)
             text(textpos(1), textpos(2), textpos(3),['Human ',subjectlist{i}(1)],'FontSize',10);
         else
             text(textpos(1), textpos(2), textpos(3),['Monkey ',subjectlist{i}(1)],'FontSize',10);
         end
      %  if (i > 1)
      %      set(gca,'Xticklabel',[],'Yticklabel',[],'ZTicklabel',[]); xlabel([]); ylabel([]); zlabel([]);
      %  end
    end
end
set(gcf,'Renderer','painters');

% Use the magic wand tool in Illustrator 
% to select all the green stuff and make it transparent.

% Uncomment this stuff to change the view angle.

h = get(gcf,'Children');
for i = 1:length(h)
    set(h(i),'View',[135 20]); % isoluminant plane
    set(h(i),'View',[0 90]); % LM plane
end



%%
% Section 11.1
% Comparing 15 Hz, 5� DTNT data across observers
% Just like above, but including lines indicating the directions of the
% blue and green phophors. (For Revision). 
ALTERNATIVEROTATION = 1;
AXISWIDTH = 2;
XMARGIN = .5;
YMARGIN = 1;
TWOMECHANISMFIT = 1;

LCONELIM = 30; % plotting limits
SCONELIM = 100; % plotting limits
% Getting the gamut of the ViewPixx for plotting
load Dell4BitsCal.mat
load T_cones_smj10

PERCGAMUT = 90; % For plotting the surface fits, how near to gamut edge to get (in each direction)
humanlist = {'GregDTNT15Hz.txt','ZackDTNT15Hz.txt','LeahDTNT15Hz.txt'};
monkeylist = {'NutDTNT15Hz.txt','FreyaDTNT15Hz.txt','SednaDTNT15Hz.txt'};
humansubjects = {'Greg','Zack','Leah'};
monkeysubjects = {'Nut','Freya','Sedna'};
if (ismac)
    filelistpath = '/Volumes/NO BACKUP/NexFiles/nexfilelists/Greg/DTNT';
else
    filelistpath = 'N:\NexFiles\nexfilelists\Greg\DTNT';
end

data = [];
BGccs = [];
subjectcounter = 1;
for HUMANS = 0:1
    if (HUMANS)
        listnames = humanlist;
    else
        listnames = monkeylist;
    end
    for i = 1:length(listnames)
        filenames = fnamesFromTxt2([filelistpath,filesep,char(listnames{i})]);
        for fileidx = 1:size(filenames,1)
            if (strncmp(filenames{fileidx},'sf',2))
                continue;
            end
            stro = nex2stro(findfile(char(filenames{fileidx})));
            [thresholds, colordirs, sfs, QuestTrajectories] = DTquestUnpackGH(stro, 'mode');
            threshs = thresholds';
            ncoldir = length(threshs);
            lmsmat = mkbasis(colordirs')';
            Loog = [];
            for j = 1:length(QuestTrajectories)
                Loog(j,1) = QuestTrajectories{j}(end) == max(QuestTrajectories{j});
            end
            data = [data; repmat(HUMANS,ncoldir,1) repmat(i,ncoldir,1) repmat(thresholds,1,3).*colordirs Loog];
            M = reshape(stro.sum.exptParams.m_mtx,3,3);
            bkgndrgb = [stro.sum.exptParams.bkgnd_r stro.sum.exptParams.bkgnd_g stro.sum.exptParams.bkgnd_b];
            bkgndlms = M*bkgndrgb(:);
            Ginclms = M*(bkgndrgb.*[1 1.1 1])';
            ccdirG = (Ginclms-bkgndlms)./bkgndlms;
            ccdirG = ccdirG./norm(ccdirG);
            Binclms = M*(bkgndrgb.*[1 1 1.1])';
            ccdirB = (Binclms-bkgndlms)./bkgndlms;
            ccdirB = ccdirB./norm(ccdirB);
            if (subjectcounter > size(BGccs,1))
                BGccs = [BGccs; ccdirG' ccdirB'];
            end
        end
        subjectcounter = subjectcounter+1;
    end
end

fundamentals = T_cones_smj10;
fundWavelengthSpacing = S_cones_smj10;
calData = cals{end};
ptb.bkgndRGB = round(255 .* calData.bgColor); %bkgnd voltages discretized b/w 0&255
ptb.bkgndrgb = [calData.gammaTable(ptb.bkgndRGB(1)+1, 1), calData.gammaTable(ptb.bkgndRGB(2)+1, 2), ...
    calData.gammaTable(ptb.bkgndRGB(3)+1, 3)]; %add one to create and index b/w 1 & 256

ptb.monSpd = SplineSpd(calData.S_device, calData.P_device, fundWavelengthSpacing);
ptb.M = fundamentals * ptb.monSpd;
ptb.bkgndlms = ptb.M * ptb.bkgndrgb';

rgb_cube_verts = fullfact([2 2 2]) - 1;
lms_cube_verts = ptb.M * rgb_cube_verts';

cc_cube_verts = bsxfun(@rdivide, bsxfun(@minus, lms_cube_verts, ptb.bkgndlms), ptb.bkgndlms)';
% K = convhull(cc_cube_verts);
% T == NaN -> points outside the gamut
lims = max(abs(data(:,[3 4 5])));
npts = 40;
[xx yy zz] = meshgrid(linspace(-lims(1),lims(1),npts),...
    linspace(-lims(2),lims(2),npts),...
    linspace(-lims(3),lims(3),npts));
xformedxyz = [xx(:) yy(:) zz(:)];
T = tsearchn(cc_cube_verts*PERCGAMUT, delaunayn(cc_cube_verts*PERCGAMUT), xformedxyz);
xformedxyz(isnan(T),:) = nan;
T = tsearchn(cc_cube_verts*PERCGAMUT, delaunayn(cc_cube_verts*PERCGAMUT), -xformedxyz);
xformedxyz(isnan(T),:) = nan;
Loog = logical(data(:,6) == 1);

% Now doing the plotting
figure('Units','inches','PaperPosition',[0 0 8.5 11],'Position',[0 0 8.5 11],'PaperOrientation','Portrait','DefaultAxesUnits','inches');
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultTextFontName','Helvetica','DefaultTextFontSize',12);
colordef white;

ccdirG = mean(BGccs(:,[1 2 3]));
ccdirB = mean(BGccs(:,[4 5 6]));
scaleG = 8;
scaleB = 20;

for HUMANS = 0:1
    if (HUMANS)
        subjectlist = humansubjects;
    else
        subjectlist = monkeysubjects;
    end
    for i = 1:length(subjectlist)
        axes('position',[1+(i-1)*(AXISWIDTH+XMARGIN),3+HUMANS*(AXISWIDTH+YMARGIN),AXISWIDTH, AXISWIDTH]);
        hold on;
        L = data(:,1) == HUMANS & data(:,2) == i;
        for j = [-1, 1]
            h = plot3(j*data(L&~Loog,3),j*data(L&~Loog,4),j*data(L&~Loog,5),'ko');
            set(h,'MarkerFaceColor','black','MarkerSize',4,'MarkerEdgeColor','none');
            h = plot3([zeros(sum(L&Loog),1) j*data(L&Loog,3)]',[zeros(sum(L&Loog),1) j*data(L&Loog,4)]',[zeros(sum(L&Loog),1) j*data(L&Loog,5)]','k-');
            set(h(:),'Color',[.5 .5 .5]);
            h = plot3(j*data(L&Loog,3),j*data(L&Loog,4),j*data(L&Loog,5),'ko');
            set(h(:),'MarkerFaceColor',[.5 .5 .5],'MarkerSize',4,'MarkerEdgeColor','none');           
        end
        % Fitting data
        [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(data(L,[3 4 5]), Loog(L));
        planeparams=(planeparams'*xformmat'); % adjusting planeparams
        A = [quadparams(1) quadparams(4) quadparams(5);... % adjusting quadparams
            quadparams(4) quadparams(2) quadparams(6);...
            quadparams(5) quadparams(6) quadparams(3)];
        B = xformmat*A*xformmat';
        quadparams = [B(1,1) B(2,2) B(3,3) B(1,2) B(1,3) B(2,3)]';
        variables = [xformedxyz(:,1).^2 xformedxyz(:,2).^2 xformedxyz(:,3).^2 2*xformedxyz(:,1).*xformedxyz(:,2) 2*xformedxyz(:,1).*xformedxyz(:,3) 2*xformedxyz(:,2).*xformedxyz(:,3)];
        [evecs, evals] = eig(B);
        initparams = [2; reshape(evecs*sqrt(evals),9,1)];
        options = optimset('MaxFunEvals',50000,'MaxIter',50000,'TolFun',10^-5,'Display','none');
        [fpar,fv] = fminsearch(@(x) colefiterr(x,data(L,[3 4 5]),Loog(L),0),initparams, options); % three mechanism model
        if (TWOMECHANISMFIT)
            [fpar,fv] = fminsearch(@(x) colefiterr(x,data(L,[3 4 5]),Loog(L),0),fpar(1:7), options); % two mechanism model
            fpar = [fpar; 0;0;0];
        end
        fr = sum(abs(xformedxyz *reshape(fpar(2:end),3,3)).^fpar(1),2);
        surfstruct = isosurface(xx,yy,zz,reshape(fr,size(xx)), 1);
        p = patch(surfstruct);
        set(p,'EdgeColor', 'none', 'FaceAlpha',.5,'FaceColor',[1 .5 .5],'Edgealpha',0);
        
        % Plotting small line segments indicating the color directions of
        % the green and blue phosphors
        [th,phi,r] = cart2sph(ccdirG(1),ccdirG(2),ccdirG(3));
        wholesum = sum(abs([cos(phi).*cos(th) cos(phi).*sin(th) sin(phi)]*reshape(fpar(2:end),3,3)).^fpar(1),2);
        predr = (1./wholesum).^(1/fpar(1));
        [ccsurfaceG(1) ccsurfaceG(2) ccsurfaceG(3)] = sph2cart(th,phi,predr);
        % Need to find the point on the surface in the ccdirB direction
        [th,phi,r] = cart2sph(ccdirB(1),ccdirB(2),ccdirB(3));
        wholesum = sum(abs([cos(phi).*cos(th) cos(phi).*sin(th) sin(phi)]*reshape(fpar(2:end),3,3)).^fpar(1),2);
        predr = (1./wholesum).^(1/fpar(1));
        [ccsurfaceB(1) ccsurfaceB(2) ccsurfaceB(3)] = sph2cart(th,phi,predr);
        
        tmp = (scaleG + norm(ccsurfaceG))./norm(ccsurfaceG);
        h1 = plot3(ccsurfaceG(1)*[1 tmp],ccsurfaceG(2)*[1 tmp],ccsurfaceG(3)*[1 tmp],'k-');
        set(h1,'Color',[0 .5 0],'Linewidth',2);
        h2 = plot3(-ccsurfaceG(1)*[1 tmp],-ccsurfaceG(2)*[1 tmp],-ccsurfaceG(3)*[1 tmp],'k-');
        set(h2,'Color',[0 .5 0],'Linewidth',2);
        
        tmp = (scaleB + norm(ccsurfaceB))./norm(ccsurfaceB);
        h3 = plot3(ccsurfaceB(1)*[1 tmp],ccsurfaceB(2)*[1 tmp],ccsurfaceB(3)*[1 tmp],'k-');
        set(h3,'Color',[0 0 1],'Linewidth',2);
        h4 = plot3(-ccsurfaceB(1)*[1 tmp],-ccsurfaceB(2)*[1 tmp],-ccsurfaceB(3)*[1 tmp],'k-');
        set(h4,'Color',[0 0 1],'Linewidth',2);
        
        % Plotting contour in the LM plane
        if (~ALTERNATIVEROTATION)
            tmpx = linspace(0,2*pi,100);
            [th,phi,r] = cart2sph(cos(tmpx)', sin(tmpx)', zeros(length(tmpx),1));
            wholesum = sum(abs([cos(phi(:)).*cos(th(:)) cos(phi(:)).*sin(th(:)) sin(phi(:))]*reshape(fpar(2:end),3,3)).^fpar(1),2);
            predr = (1./wholesum).^(1/fpar(1));
            contour =[predr.*cos(tmpx)'  predr.*sin(tmpx)'];
            T = tsearchn(cc_cube_verts*PERCGAMUT, delaunayn(cc_cube_verts*PERCGAMUT), [contour zeros(size(contour,1),1)]);
            contour(isnan(T),:) = nan;
            plot(contour(:,1),contour(:,2),'LineWidth',3,'Color',[0 .3 0]);
            set(gca,'View',[0 90]);
        else
            if (subjectlist{i}(1) == 'G')
                set(gca,'View',[48 -8]);
            elseif (subjectlist{i}(1) == 'Z')
                set(gca,'View',[5 -2]);
            elseif (subjectlist{i}(1) == 'L')
                set(gca,'View',[5 -2]);
            elseif (subjectlist{i}(1) == 'N')
                set(gca,'View',[42 -8]);
            elseif (subjectlist{i}(1) == 'F')
                set(gca,'View',[42 -8]);
            elseif (subjectlist{i}(1) == 'S')
                set(gca,'View',[62 -6]);
            end
        end
        set(gca,'Xlim',[-LCONELIM LCONELIM],'Ylim',[-LCONELIM LCONELIM],'Zlim',[-SCONELIM SCONELIM]);
        set(gca,'Xtick',[-LCONELIM 0 LCONELIM],'Xticklabel',num2str([-LCONELIM 0 LCONELIM]'/100))
        xlabel('\DeltaL/L');
        set(gca,'Ytick',[-LCONELIM 0 LCONELIM],'Yticklabel',num2str([-LCONELIM 0 LCONELIM]'/100))
        ylabel('\DeltaM/M');
        set(gca,'Ztick',[-SCONELIM 0 SCONELIM],'Zticklabel',num2str([-SCONELIM 0 SCONELIM]'/100))
        zlabel('\DeltaS/S');
        if (~ALTERNATIVEROTATION)
            textpos = [10, -28, 0];
            if (HUMANS)
                text(textpos(1), textpos(2), textpos(3),['Human ',subjectlist{i}(1)],'FontSize',10);
            else
                text(textpos(1), textpos(2), textpos(3),['Monkey ',subjectlist{i}(1)],'FontSize',10);
            end
        end
      %  if (i > 1)
      %      set(gca,'Xticklabel',[],'Yticklabel',[],'ZTicklabel',[]); xlabel([]); ylabel([]); zlabel([]);
      %  end
    end
    drawnow
end

set(gcf,'Renderer','painters');

%%
% Section 11.2 SEs of beta from bootstrap
% Also estimating the contribution of "nonluminance" mechanisms to the
% detection of the blue and green phosphor
%
% Cone contrast direction for green [0.6306 0.7715 0.0847]
% and blue [ 0.1118 0.2001 0.9734]

% Getting the gamut of the ViewPixx for plotting
%load Dell4BitsCal.mat
%load T_cones_smj10
TWOMECHANISMFIT = 1; % 0 = 3 mechanism fit, 1 = 2 mechanism fit

humanlist = {'GregDTNT15Hz.txt','ZackDTNT15Hz.txt','LeahDTNT15Hz.txt'};
monkeylist = {'NutDTNT15Hz.txt','FreyaDTNT15Hz.txt','SednaDTNT15Hz.txt'};
humansubjects = {'Greg','Zack','Leah'};
monkeysubjects = {'Nut','Freya','Sedna'};
if (ismac)
    filelistpath = '/Volumes/NO BACKUP/NexFiles/nexfilelists/Greg/DTNT';
else
    filelistpath = 'N:\NexFiles\nexfilelists\Greg\DTNT';
end

data = [];
for HUMANS = 0:1
    if (HUMANS)
        listnames = humanlist;
    else
        listnames = monkeylist;
    end
    for i = 1:length(listnames)
        filenames = fnamesFromTxt2([filelistpath,filesep,char(listnames{i})]);
        for fileidx = 1:size(filenames,1)
            if (strncmp(filenames{fileidx},'sf',2))
                continue;
            end
            stro = nex2stro(findfile(char(filenames{fileidx})));
            [thresholds, colordirs, sfs, QuestTrajectories] = DTquestUnpackGH(stro, 'mode');
            threshs = thresholds';
            ncoldir = length(threshs);
            lmsmat = mkbasis(colordirs')';
            Loog = [];
            for j = 1:length(QuestTrajectories)
                Loog(j,1) = QuestTrajectories{j}(end) == max(QuestTrajectories{j});
            end
            data = [data; repmat(HUMANS,ncoldir,1) repmat(i,ncoldir,1) repmat(thresholds,1,3).*colordirs Loog];
        end
    end
end

Loog = logical(data(:,6) == 1);
niter = 2;
beta_se = [];
fpars = [];
for HUMANS = 0:1
    if (HUMANS)
        subjectlist = humansubjects;
    else
        subjectlist = monkeysubjects;
    end
    for i = 1:length(subjectlist)
        L = data(:,1) == HUMANS & data(:,2) == i;
        % Fitting data
        [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(data(L,[3 4 5]), Loog(L));
        planeparams=(planeparams'*xformmat'); % adjusting planeparams
        A = [quadparams(1) quadparams(4) quadparams(5);... % adjusting quadparams
            quadparams(4) quadparams(2) quadparams(6);...
            quadparams(5) quadparams(6) quadparams(3)];
        B = xformmat*A*xformmat';
        quadparams = [B(1,1) B(2,2) B(3,3) B(1,2) B(1,3) B(2,3)]';
        [evecs, evals] = eig(B);
        initparams = [2; reshape(evecs*sqrt(evals),9,1)];
        options = optimset('MaxFunEvals',50000,'MaxIter',50000,'TolFun',10^-5,'Display','none');
        [fpar,fv] = fminsearch(@(x) colefiterr(x,data(L,[3 4 5]),Loog(L),0),initparams, options);
        if (TWOMECHANISMFIT)
            [fpar,fv] = fminsearch(@(x) colefiterr(x,data(L,[3 4 5]),Loog(L),0),fpar(1:7), options); % two mechanism model
            fpar = [fpar; 0;0;0];
        end
        
        which_idxs = find(L);
        bootdata = nan*ones(niter,1);
        for j = 1:niter
            boot_idxs = which_idxs(unidrnd(length(which_idxs),length(which_idxs),1));
            [fpar_boot,~] = fminsearch(@(x) colefiterr(x,data(boot_idxs,[3 4 5]),Loog(boot_idxs),0),fpar, options);
            bootdata(j) = fpar_boot(1);
            bootdata(j)
        end
        figure;
        subplot(2,1,1); hist(bootdata);
        subplot(2,1,2); hist(log10(bootdata));
        beta_se =[beta_se; fpar(1) std(bootdata)]
        fpars = [fpars, fpar];
    end
end

% Trying to predict B/G thresholds on the basis of a single mechanism (the
% dominant one)

Gdir = [0.6306 0.7715 0.0847];
Bdir = [0.1118 0.2001 0.9734];

predrs = [];
for subject = 1:6;
    fpar = fpars(:,subject);
    beta = fpar(1);
    mech1 = fpar(2:4);
    mech2 = fpar(5:7);
    if sign(mech1(1)) ~= sign(mech1(2))
        mech1 = fpar(5:7);
        mech2 = fpar(2:4);
    end
    
    if (TWOMECHANISMFIT)
        mech3 = fpar(8:10);
    end
    
    [th,phi,r] = cart2sph([Gdir(1);Bdir(1)],[Gdir(2);Bdir(2)],[Gdir(3);Bdir(3)]);
    if (TWOMECHANISMFIT)
        wholesum = sum(abs([cos(phi(:)).*cos(th(:)) cos(phi(:)).*sin(th(:)) sin(phi(:))]*[mech1,mech2,mech3]).^beta,2);        
    else
        wholesum = sum(abs([cos(phi(:)).*cos(th(:)) cos(phi(:)).*sin(th(:)) sin(phi(:))]*[mech1,mech2]).^beta,2);
    end
    predr_big = (1./wholesum).^(1/beta); % predictied cone contrast thresholds with 2 mechanisms
    
    wholesum = sum(abs([cos(phi(:)).*cos(th(:)) cos(phi(:)).*sin(th(:)) sin(phi(:))]*mech1).^beta,2);
    predr_small = (1./wholesum).^(1/beta);% predicted cone contrast thresholds with 1 mechanism
    predrs = cat(3,predrs,[predr_small predr_big])
end
BGs = squeeze(predrs(1,:,:)./predrs(2,:,:));
%howmuchBGchanges = abs(1-(BGs(1,:)./BGs(2,:)))  % How much BG changes
howmuchBGchanges = (BGs(1,:)-BGs(2,:))./BGs(1,:)
geomean(howmuchBGchanges(1:3))

howmuchgreenchanges = 1-squeeze(predrs(1,2,:)./predrs(1,1,:));
howmuchbluechanges = 1-squeeze(predrs(2,2,:)./predrs(2,1,:));
geomean([howmuchgreenchanges(1:3); howmuchbluechanges(1:3)]) % monkeys only

%%
% Section 11.3
% Seeing how much B/G changes between 25 Hz DTNT and MacPig data
% (to justify the idea that B/G at 15 Hz are detected by the same
% mechanism)

DTNTsuffix = 'DTNT25Hz.txt';
macpigsuffix = 'MacPig.txt';
subjects = {'Nut','Freya','Sedna'};
if (ismac)
    DTNTfilelistpath = '/Volumes/NO BACKUP/NexFiles/nexfilelists/Greg/DTNT';
    macpigfilelistpath = '/Volumes/NO BACKUP/NexFiles/nexfilelists/Greg/DTEM';
else
    DTNTfilelistpath = 'N:\NexFiles\nexfilelists\Greg\DTNT';
    macpigfilelistpath = 'N:\NexFiles\nexfilelists\Greg\DTEM';
end
Gdir = [0.6306 0.7715 0.0847];
Bdir = [0.1118 0.2001 0.9734];

load Vlambda_star
Vlambda = Vlambda_star(:,2)';
macpigdata = [];
DTNTdata = [];
for i = 1:length(subjects)
    i
    % First getting DTNT data
    tmpdata = [];
    filenames = fnamesFromTxt2([DTNTfilelistpath,filesep,char(subjects{i}),DTNTsuffix]);
    for fileidx = 1:size(filenames,1)
        if (strncmp(filenames{fileidx},'sf',2))
            continue;
        end
        stro = nex2stro(findfile(char(filenames{fileidx})));
        [thresholds, colordirs, sfs, QuestTrajectories] = DTquestUnpackGH(stro, 'mode');
        threshs = thresholds';
        ncoldir = length(threshs);
        lmsmat = mkbasis(colordirs')';
        Loog = [];
        for j = 1:length(QuestTrajectories)
            Loog(j,1) = QuestTrajectories{j}(end) == max(QuestTrajectories{j});
        end
        tmpdata = [tmpdata; repmat(i,ncoldir,1) repmat(thresholds,1,3).*colordirs Loog];
    end
    Loog = logical(tmpdata(:,end) == 1);
    % Fitting data
    [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(tmpdata(:,[2 3 4]), Loog);
    lumvect = planeparams'*xformmat'
    DTNTdata = [DTNTdata; lumvect*Gdir' lumvect*Bdir'];
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Getting macpig data
    tmpdata = [];
    filenames = fnamesFromTxt2([macpigfilelistpath,filesep,char(subjects{i}),macpigsuffix]);
    for fileidx = 1:size(filenames,1)
        stro = nex2stro(findfile(char(filenames{fileidx,:})));
        [thresholds, colordirs, sfs] = DTquestUnpackGH(stro, 'mode');
        threshs = thresholds';
        ncoldir = length(threshs);
        lmsmat = mkbasis(colordirs')';       
        ecc = stro.sum.exptParams.rf_x;
        vectornormscc = sqrt(sum((repmat(thresholds,1,3).*lmsmat).^2,2));
        tmpdata = [tmpdata; ecc/10 vectornormscc'];
    end
    % Trimming the ecentricities that are too close to the fovea
    L = tmpdata(:,1) >= 5;
    tmpdata(~L,:) = [];
    macpigdata = [macpigdata; geomean(tmpdata(:,2)) geomean(tmpdata(:,3))];
end
% macpig data are thresholds, DTNT data are sensitivities
mpBG = macpigdata(:,2)./macpigdata(:,1);
dtntBG = DTNTdata(:,1)./DTNTdata(:,2);
geomean(1-dtntBG./mpBG)

% Nut, Freya, Sedna
% 25 Hz, MacPig, 15 Hz(1-mechanism), 15 Hz(2-mechanisms)
data = [4.7515 4.2661 3.3133; 5.1503 4.5508 3.3751; 4.2801 4.1719 2.9346; 5.5387 4.8526 3.4258]
geomean(1-data(1,:)./data(2,:)) % 25 Hz vs Macpig
geomean(abs(1-data(3,:)./data(2,:)))  % 15 Hz (1-mechanism) vs Macpig
geomean(abs(1-data(4,:)./data(2,:)))  % 15 Hz (2-mechanism) vs Macpig
geomean(abs(1-data(1,:)./data(4,:)))  % 15 Hz (2-mechanism) vs 25 Hz

%%
% Section 12
% Plotting 15 Hz blue and green flicker thresholds on top of 15 Hz DTNT
% isodetection surfaces

AXISWIDTH = 2;
XMARGIN = .25;
YMARGIN = 1;

LCONELIM = 10; % plotting limits
SCONELIM = 100; % plotting limits

% Getting the gamut for plotting sufaces
load Dell4BitsCal.mat
load T_cones_smj10
PERCGAMUT = 90; % For plotting the surface fits, how near to gamut edge to get (in each direction)
humanlist = {'GregDTNT15Hz.txt','ZackDTNT15Hz.txt','LeahDTNT15Hz.txt'};
monkeylist = {'NutDTNT15Hz.txt','FreyaDTNT15Hz.txt','SednaDTNT15Hz.txt'};
humansubjects = {'Zack','Greg','Leah'};
monkeysubjects = {'Nut','Freya','Sedna'};

if (ismac)
    flpdtnt = '/Volumes/NO BACKUP/NexFiles/nexfilelists/Greg/DTNT';
    flpmacpig = '/Volumes/NO BACKUP/NexFiles/nexfilelists/Greg/DTEM';
else
    flpdtnt = 'N:\NexFiles\nexfilelists\Greg\DTNT';
    flpmacpig = 'N:\NexFiles\nexfilelists\Greg\DTEM';
end

% First getting the DTNT data.
% Identical to the cell above.
dtntdata = [];
for HUMANS = 0:1
    if (HUMANS)
        listnames = humanlist;
    else
        listnames = monkeylist;
    end
    for i = 1:length(listnames)
        filenames = fnamesFromTxt2([flpdtnt,filesep,char(listnames{i})]);
        for fileidx = 1:size(filenames,1)
            if (strncmp(filenames{fileidx},'sf',2))
                continue;
            end
            stro = nex2stro(findfile(char(filenames{fileidx})));
            [thresholds, colordirs, sfs, QuestTrajectories] = DTquestUnpackGH(stro, 'mode');
            threshs = thresholds';
            ncoldir = length(threshs);
            lmsmat = mkbasis(colordirs')';
            Loog = [];
            for j = 1:length(QuestTrajectories)
                Loog(j,1) = QuestTrajectories{j}(end) == max(QuestTrajectories{j});
            end
            dtntdata = [dtntdata; repmat(HUMANS,ncoldir,1) repmat(i,ncoldir,1) repmat(thresholds,1,3).*colordirs Loog];
        end
    end
end

% Now getting the MacPig data
macpigdata = [];
for HUMANS = 0:1
    if (HUMANS)
        subjectnames = humansubjects;
    else
        subjectnames = monkeysubjects;
    end
    for i = 1:length(listnames)
        macpigfilename = [char(subjectnames{i}),'MacPig.txt'];
        filenames = fnamesFromTxt2([flpmacpig,filesep,macpigfilename]);
        
        for file = filenames'
            stro = nex2stro(findfile(char(file{:})));
            if stro.sum.exptParams.rf_x == 50 % only get the 5 degree eccentric thresholds
                [threshes, color_dirs] = DTquestUnpackGH(stro, 'mode');
                macpig_scaled = bsxfun(@times, color_dirs, threshes./sqrt(sum(color_dirs.^2, 2))); % make color_dirs unit vectors and multiply by thresholds
                blueidx = find(color_dirs(:,3) == max(color_dirs(:,3))); % blue has the greatest S-cone component
                greenidx = find(color_dirs(:,2) == max(color_dirs(:,2)));% green has the greatest M-cone component
                macpigdata = [macpigdata; repmat([HUMANS i],[1 1 2]) permute([macpig_scaled(greenidx,:); macpig_scaled(blueidx,:)],[3 2 1])];
            end
        end
    end
end

% Stuff for rending surface
fundamentals = T_cones_smj10;
fundWavelengthSpacing = S_cones_smj10;
calData = cals{end};
ptb.bkgndRGB = round(255 .* calData.bgColor); %bkgnd voltages discretized b/w 0&255
ptb.bkgndrgb = [calData.gammaTable(ptb.bkgndRGB(1)+1, 1), calData.gammaTable(ptb.bkgndRGB(2)+1, 2), ...
    calData.gammaTable(ptb.bkgndRGB(3)+1, 3)]; %add one to create and index b/w 1 & 256

ptb.monSpd = SplineSpd(calData.S_device, calData.P_device, fundWavelengthSpacing);
ptb.M = fundamentals * ptb.monSpd;
ptb.bkgndlms = ptb.M * ptb.bkgndrgb';

rgb_cube_verts = fullfact([2 2 2]) - 1;
lms_cube_verts = ptb.M * rgb_cube_verts';

cc_cube_verts = bsxfun(@rdivide, bsxfun(@minus, lms_cube_verts, ptb.bkgndlms), ptb.bkgndlms)';


% Plotting
figure('Units','inches','PaperPosition',[0 0 8.5 11],'Position',[0 0 8.5 11],'PaperOrientation','Portrait','DefaultAxesUnits','inches');
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultTextFontName','Helvetica','DefaultTextFontSize',12);
colordef white;

for HUMANS = 0:1
    if (HUMANS)
        listnames = humanlist;
    else
        listnames = monkeylist;
    end
    for i = 1:length(listnames)
        L = dtntdata(:,1) == HUMANS & dtntdata(:,2) == i;
        [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(dtntdata(L,[3 4 5]), dtntdata(L,6));
        planeparams=(planeparams'*xformmat'); % adjusting planeparams
        A = [quadparams(1) quadparams(4) quadparams(5);... % adjusting quadparams
            quadparams(4) quadparams(2) quadparams(6);...
            quadparams(5) quadparams(6) quadparams(3)];
        B = xformmat*A*xformmat';
        quadparams = [B(1,1) B(2,2) B(3,3) B(1,2) B(1,3) B(2,3)]';
        [evecs, evals] = eig(B);
        initparams = [2; reshape(evecs*sqrt(evals),9,1)];
        options = optimset('MaxFunEvals',50000,'MaxIter',50000,'TolFun',10^-5,'Display','none');
        [fpar,fv] = fminsearch(@(x) colefiterr(x,dtntdata(L,[3 4 5]),dtntdata(L,6),0),initparams, options);
        lims = max(abs(dtntdata(:,[3 4 5])));
        npts = 100;
        [xx, yy, zz] = meshgrid(linspace(-lims(1),lims(1),npts),...
            linspace(-lims(2),lims(2),npts),...
            linspace(-lims(3),lims(3),npts));
        xformedxyz = [xx(:) yy(:) zz(:)];
        T = tsearchn(cc_cube_verts*PERCGAMUT, delaunayn(cc_cube_verts*PERCGAMUT), xformedxyz);
        xformedxyz(isnan(T),:) = nan;
        T = tsearchn(cc_cube_verts*PERCGAMUT, delaunayn(cc_cube_verts*PERCGAMUT), -xformedxyz);
        xformedxyz(isnan(T),:) = nan;
        fr = sum(abs(xformedxyz *reshape(fpar(2:end),3,3)).^fpar(1),2);
        surfstruct = isosurface(xx,yy,zz,reshape(fr,size(xx)), 1);
        
        % Plotting MacPig data
        L = macpigdata(:,1) == HUMANS & macpigdata(:,2) == i;
        [~,s,v] = svd(macpigdata(L,[3 4 5],1));
        %diag(s)
        unitvect_g =  sign(sum(v(:,1)))*v(:,1);
        [~,s,v] = svd(macpigdata(L,[3 4 5],2));
        %diag(s)
        unitvect_b =  sign(sum(v(:,1)))*v(:,1);
        nullvect = null([unitvect_b'; unitvect_g']); % orthogonal to b and g. Should point into the page of the plot.
        
        vert = surfstruct.vertices;
        nvproj = vert*nullvect;
        TOL = .25;
        Lproj = abs(nvproj)<TOL;
        axes('position',[1+(i-1)*(AXISWIDTH+XMARGIN),3+HUMANS*(AXISWIDTH+YMARGIN),AXISWIDTH, AXISWIDTH]); 
        hold on;
        curvexy = vert(Lproj,:)*[unitvect_g unitvect_b];
        whichclust = cluster(linkage(curvexy),'MaxClust',2);
        for j = 1:2 % plotting spline fits
            cs = csaps(curvexy(whichclust == j,2), curvexy(whichclust == j,1),.01);
            xx = linspace(min(curvexy(whichclust == j,2)),max(curvexy(whichclust == j,2)),300)';
            xyz = inv ([unitvect_g unitvect_b nullvect]')*[fnval(cs,xx),xx,zeros(length(xx),1)]';
            plot3(xyz(1,:),xyz(2,:),xyz(3,:),'k-','linewidth',2);
        end
  
      %  plot3(vert(Lproj,1),vert(Lproj,2),vert(Lproj,3),'y.');
        c = {'green','blue'};
        for j = 1:size(macpigdata,3) % looping over green, blue
            h(1) = plot3(macpigdata(L,3,j), macpigdata(L,4,j), macpigdata(L,5,j), 'mp', 'MarkerSize', 10,'MarkerFaceColor',c{j},'MarkerEdgeColor',c{j});
            h(2) = plot3(-macpigdata(L,3,j), -macpigdata(L,4,j), -macpigdata(L,5,j), 'mp', 'MarkerSize', 10,'MarkerFaceColor',c{j},'MarkerEdgeColor',c{j});
        end
        view(nullvect);
        set(gca,'Xlim',[-LCONELIM LCONELIM],'Ylim',[-LCONELIM LCONELIM],'Zlim',[-SCONELIM SCONELIM]);
        set(gca,'Xtick',[-LCONELIM 0 LCONELIM]/2,'Xticklabel',num2str([-LCONELIM 0 LCONELIM]'/200))
        xlabel('\DeltaL/L');
        set(gca,'Ytick',[-LCONELIM 0 LCONELIM]/2,'Yticklabel',num2str([-LCONELIM 0 LCONELIM]'/200))
        ylabel('\DeltaM/M');
        set(gca,'Ztick',[-SCONELIM 0 SCONELIM],'Zticklabel',num2str([-SCONELIM 0 SCONELIM]'/100))
        zlabel('\DeltaS/S');
    end
end

%%
% Section 12.1
% isodeteability contour and blue/green macpig data rendered in a
% color plane spanned by orthogonal blue-gun and green-gun axes.

AXISWIDTH = 2;
XMARGIN = .25;
YMARGIN = 1;
load T_xyz1931;
Vlambda = T_xyz1931(2,:);

humanlist = {'GregDTNT15Hz.txt','ZackDTNT15Hz.txt','LeahDTNT15Hz.txt'};
monkeylist = {'NutDTNT15Hz.txt','FreyaDTNT15Hz.txt','SednaDTNT15Hz.txt'};

if (ismac)
    flpdtnt = '/Volumes/NO BACKUP/NexFiles/nexfilelists/Greg/DTNT';
    flpmacpig = '/Volumes/NO BACKUP/NexFiles/nexfilelists/Greg/DTEM';
else
    flpdtnt = 'N:\NexFiles\nexfilelists\Greg\DTNT';
    flpmacpig = 'N:\NexFiles\nexfilelists\Greg\DTEM';
end

% First getting the DTNT data.
% Identical to the cell above.
dtntdata = [];
for HUMANS = 0:1
    if (HUMANS)
        listnames = humanlist;
    else
        listnames = monkeylist;
    end
    for i = 1:length(listnames)
        filenames = fnamesFromTxt2([flpdtnt,filesep,char(listnames{i})]);
        for fileidx = 1:size(filenames,1)
            if (strncmp(filenames{fileidx},'sf',2))
                continue;
            end
            stro = nex2stro(findfile(char(filenames{fileidx})));
            [thresholds, colordirs, sfs, QuestTrajectories] = DTquestUnpackGH(stro, 'mode');
            threshs = thresholds';
            ncoldir = length(threshs);
            lmsmat = mkbasis(colordirs')';
            Loog = [];
            for j = 1:length(QuestTrajectories)
                Loog(j,1) = QuestTrajectories{j}(end) == max(QuestTrajectories{j});
            end
            dtntdata = [dtntdata; repmat(HUMANS,ncoldir,1) repmat(i,ncoldir,1) repmat(thresholds/100,1,3).*colordirs Loog];
        end
        monitordata(HUMANS+1,i).M = reshape(stro.sum.exptParams.m_mtx,3,3);
        monitordata(HUMANS+1,i).spd = reshape(stro.sum.exptParams.mon_spect,length(stro.sum.exptParams.mon_spect)/3,3);
        monitordata(HUMANS+1,i).bkgndrgb = [stro.sum.exptParams.bkgnd_r stro.sum.exptParams.bkgnd_g stro.sum.exptParams.bkgnd_b];    
    end
end

humansubjects = {};
monkeysubjects = {};
% A bit of preprocessing
suffix = 'DTNT15Hz.txt';
for i = 1:length(humanlist)
    if (~strncmp(fliplr(humanlist{i}),fliplr(suffix),length(suffix)))
        error('filename error');
    else
        humansubjects{i} = humanlist{i}(1:end-length(suffix));
    end
end

for i = 1:length(monkeylist)
    if (~strncmp(fliplr(monkeylist{i}),fliplr(suffix),length(suffix)))
        error('filename error');
    else
        monkeysubjects{i} = monkeylist{i}(1:end-length(suffix));
    end
end

% Now getting the MacPig data
macpigdata = [];
for HUMANS = 0:1
    if (HUMANS)
        subjectnames = humansubjects;
    else
        subjectnames = monkeysubjects;
    end
    for i = 1:length(subjectnames)
        macpigfilename = [char(subjectnames{i}),'MacPig.txt'];
        filenames = fnamesFromTxt2([flpmacpig,filesep,macpigfilename]);
        
        for file = filenames'
            stro = nex2stro(findfile(char(file{:})));
            if stro.sum.exptParams.rf_x == 50 % only get the 5 degree eccentric thresholds
                [thresh, colorDirs, sfs] = DTquestUnpack(stro, 'mode'); close(gcf);
                thresh = thresh./100;
                macpig_scaled = bsxfun(@times, colorDirs, thresh./sqrt(sum(colorDirs.^2, 2))); % make color_dirs unit vectors and multiply by thresholds
                blueidx = find(colorDirs(:,3) == max(colorDirs(:,3))); % blue has the greatest S-cone component
                greenidx = find(colorDirs(:,2) == max(colorDirs(:,2)));% green has the greatest M-cone component
                macpigdata = [macpigdata; repmat([HUMANS i],[1 1 2]) permute([macpig_scaled(greenidx,:); macpig_scaled(blueidx,:)],[3 2 1])];
                % Get rid of factor of 100 above?
                
                
              %  normcolordirs = colorDirs./repmat(sqrt(sum(colorDirs.^2,2)),1,3);
              %  guns = reshape(stro.sum.exptParams.mon_spect,81,3);
              %  bkgndrgb = [stro.sum.exptParams.bkgnd_r stro.sum.exptParams.bkgnd_g stro.sum.exptParams.bkgnd_b];
              %  M = reshape(stro.sum.exptParams.m_mtx,3,3);                
              %  bkgndlms = M*bkgndrgb';
              %  bkgndlum = bkgndrgb*guns'*Vlambda';  
              %  threshrgb = inv(M)*(repmat(bkgndlms',2,1).*(1+repmat(thresh/100,[1 3]).*normcolordirs([1 2],:)))';
              %  deltagun = threshrgb-repmat(bkgndrgb',1,2);
                % Quantifying thresholds in 2 deg luminance contrast units
          %      thresholds = [nan; nan];
          %      for j = 1:size(deltagun,2)
          %          thresholds(j) = (guns*deltagun(:,j))'*Vlambda';
          %          thresholds(j) = thresholds(j)./bkgndlum;
          %      end
          %      macpigdata = [macpigdata; HUMANS i thresholds'.*100];
            end
        end
    end
end


% Now doing the plotting
AXISWIDTH = 2;
XMARGIN = .5;
YMARGIN = 2;

figure('Units','inches','PaperPosition',[0 0 8.5 11],'Position',[0 0 8.5 11],'PaperOrientation','Portrait','DefaultAxesUnits','inches');
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultTextFontName','Helvetica','DefaultTextFontSize',12);
colordef white;

for HUMANS = 0:1
    if (HUMANS)
        subjectnames = humansubjects;
    else
        subjectnames = monkeysubjects;
    end
    for i = 1:length(subjectnames)
        axes('position',[1+(i-1)*(AXISWIDTH+XMARGIN),1+HUMANS*(AXISWIDTH+YMARGIN),AXISWIDTH, AXISWIDTH]);
        hold on;
        % First macpig data
        L = macpigdata(:,1) == HUMANS & macpigdata(:,2) == i;
        normg = sqrt(sum(macpigdata(L,[3 4 5],1).^2,2));
        [u,s,v] = svd(macpigdata(L,[3 4 5],1));
        ccdirg = v(:,1);
        normb = sqrt(sum(macpigdata(L,[3 4 5],2).^2,2));  
        [u,s,v] = svd(macpigdata(L,[3 4 5],2));
        ccdirb = v(:,1);
        
%         plot(log10(normg),0,'kp','markerfacecolor','green','Markersize',12);
 %        plot(-log10(normg),0,'kp','markerfacecolor','green','Markersize',12);
 %        plot(0,log10(normb),'kp','markerfacecolor','blue','Markersize',12);
 %        plot(0,-log10(normb),'kp','markerfacecolor','blue','Markersize',12);

         plot(normg,0,'kp','markerfacecolor','green','Markersize',12);
         plot(-normg,0,'kp','markerfacecolor','green','Markersize',12);
         plot(0,normb,'kp','markerfacecolor','blue','Markersize',12);
         plot(0,-normb,'kp','markerfacecolor','blue','Markersize',12);
        
       % plot(macpigdata(L,3),0,'kp','markerfacecolor','green','Markersize',12);
       % plot(0,macpigdata(L,4),'kp','markerfacecolor','blue','Markersize',12);
       % plot(-macpigdata(L,3),0,'kp','markerfacecolor','green','Markersize',12);
       % plot(0,-macpigdata(L,4),'kp','markerfacecolor','blue','Markersize',12);
       % set(gca,'Xlim',max(macpigdata(:,3))*1.22*[-1 1],'Ylim',max(macpigdata(:,3))*1.22*[-1 1]);
       % ylabel('Blue threshold (% contrast)');
       % xlabel('Green threshold (% contrast)');

        % Now DTNT data
        L = dtntdata(:,1) == HUMANS & dtntdata(:,2) == i;
        [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(dtntdata(L,[3 4 5]), dtntdata(L,6));
        planeparams=(planeparams'*xformmat'); % adjusting planeparams
        A = [quadparams(1) quadparams(4) quadparams(5);... % adjusting quadparams
            quadparams(4) quadparams(2) quadparams(6);...
            quadparams(5) quadparams(6) quadparams(3)];
        B = xformmat*A*xformmat';
        quadparams = [B(1,1) B(2,2) B(3,3) B(1,2) B(1,3) B(2,3)]';
        [evecs, evals] = eig(B);
        initparams = [2; reshape(evecs*sqrt(evals),9,1)];
        options = optimset('MaxFunEvals',50000,'MaxIter',50000,'TolFun',10^-5,'Display','none');
        [fpar,fv] = fminsearch(@(x) colefiterr(x,dtntdata(L,[3 4 5]),dtntdata(L,6),0),initparams, options);
        %spd = monitordata(HUMANS+1,i).spd;
        M = monitordata(HUMANS+1,i).M;
        ndirs = 150;
       % fpar(5:end) = 0; % DEBUGGING
       % fpar(1) = 1; % DUBUGGING
        
        cc = [ccdirg, ccdirb]*[cos(linspace(0,pi*2,ndirs)); sin(linspace(0,pi*2,ndirs))]; % in cc space
        ccnorms = sqrt(sum(cc.^2));
        cc = cc./repmat(ccnorms,3,1);
        wholesum = sum(abs([cc(1,:)', cc(2,:)',cc(3,:)']*reshape(fpar(2:end),3,3)).^fpar(1),2);
        predr = (1./wholesum').^(1/fpar(1));
        % "r" in cone contrast space
        % Don't plot pieces of the curve that are out of gamut
        ingamut = [];
        for j = 1:ndirs
            ingamut(j) = gamutCheck(predr(j)*[cc(1,j),cc(2,j),cc(3,j)], monitordata(HUMANS+1,i).bkgndrgb, M, 'both');
        end
        predr(~ingamut) = nan;
        % plot(cos(linspace(0,2*pi,ndirs)).*predr,sin(linspace(0,2*pi,ndirs)).*predr,'k-')
        x = cos(linspace(0,pi*2,ndirs)).*predr./ccnorms;
        y = sin(linspace(0,pi*2,ndirs)).*predr./ccnorms;
        plot(x,y,'k-','LineWidth',2);
        set(gca,'Xlim',[-1.25 1.25]*max(abs(normg)),'Ylim', [-1.25 1.25]*max(abs(normb)));
        drawnow
    end
end


%%
% Section 12.2
% 15 Hz isodetecability surface and plane defined by the the blue/green
% phosphor axes. Showing intersection. (For Reviewer 2). Single subject
% only

dtntdata = [];
if (ismac)
    flpdtnt = '/Volumes/NO BACKUP/NexFiles/nexfilelists/Greg/DTNT';
else
    flpdtnt = 'N:\NexFiles\nexfilelists\Greg\DTNT';
end

filenames = fnamesFromTxt2([flpdtnt,filesep,'FreyaDTNT15Hz.txt']);
for fileidx = 1:size(filenames,1)
    if (strncmp(filenames{fileidx},'sf',2))
        continue;
    end
    stro = nex2stro(findfile(char(filenames{fileidx})));
    [thresholds, colordirs, sfs, QuestTrajectories] = DTquestUnpackGH(stro, 'mode');
    threshs = thresholds';
    ncoldir = length(threshs);
    lmsmat = mkbasis(colordirs')';
    Loog = [];
    for j = 1:length(QuestTrajectories)
        Loog(j,1) = QuestTrajectories{j}(end) == max(QuestTrajectories{j});
    end
    dtntdata = [dtntdata; repmat(thresholds/100,1,3).*colordirs Loog];
    if (fileidx == 1)
        M = reshape(stro.sum.exptParams.m_mtx,3,3);
        spd = reshape(stro.sum.exptParams.mon_spect,length(stro.sum.exptParams.mon_spect)/3,3);
        bkgndrgb = [stro.sum.exptParams.bkgnd_r stro.sum.exptParams.bkgnd_g stro.sum.exptParams.bkgnd_b]';  
    end
end

[planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(dtntdata(:,[1 2 3]), dtntdata(:,4));
planeparams=(planeparams'*xformmat'); % adjusting planeparams
A = [quadparams(1) quadparams(4) quadparams(5);... % adjusting quadparams
    quadparams(4) quadparams(2) quadparams(6);...
    quadparams(5) quadparams(6) quadparams(3)];
B = xformmat*A*xformmat';
quadparams = [B(1,1) B(2,2) B(3,3) B(1,2) B(1,3) B(2,3)]';
[evecs, evals] = eig(B);
initparams = [2; reshape(evecs*sqrt(evals),9,1)];
options = optimset('MaxFunEvals',50000,'MaxIter',50000,'TolFun',10^-5,'Display','none');
[fpar,fv] = fminsearch(@(x) colefiterr(x,dtntdata(:,[1 2 3]), dtntdata(:,4),0),initparams, options);

% Stuff for rending surface
rgb_cube_verts = fullfact([2 2 2]) - 1;
lms_cube_verts = M * rgb_cube_verts';
bkgndlms = M*bkgndrgb(:);
cc_cube_verts = bsxfun(@rdivide, bsxfun(@minus, lms_cube_verts, bkgndlms), bkgndlms)';
PERCGAMUT = .9;

lims = max(abs(dtntdata(:,[1 2 3])));
npts = 100;
[xx, yy, zz] = meshgrid(linspace(-lims(1),lims(1),npts),...
    linspace(-lims(2),lims(2),npts),...
    linspace(-lims(3),lims(3),npts));
xformedxyz = [xx(:) yy(:) zz(:)];
T = tsearchn(cc_cube_verts*PERCGAMUT, delaunayn(cc_cube_verts*PERCGAMUT), xformedxyz);
xformedxyz(isnan(T),:) = nan;
T = tsearchn(cc_cube_verts*PERCGAMUT, delaunayn(cc_cube_verts*PERCGAMUT), -xformedxyz);
xformedxyz(isnan(T),:) = nan;
fr = sum(abs(xformedxyz *reshape(fpar(2:end),3,3)).^fpar(1),2);
surfstruct = isosurface(xx,yy,zz,reshape(fr,size(xx)), 1);

figure; axes; hold on;
p = patch(surfstruct);
set(p,'EdgeColor', 'none', 'FaceAlpha',.5,'FaceColor','green','Edgealpha',0);

% Now drawing B and G color directions
ccdirG = (M(:,2)-bkgndlms)./bkgndlms;
ccdirB = (M(:,3)-bkgndlms)./bkgndlms;
scale = 1.5;

% Need to find the point on the surface in the ccdirG direction
[th,phi,r] = cart2sph(ccdirG(1),ccdirG(2),ccdirG(3));
wholesum = sum(abs([cos(phi).*cos(th) cos(phi).*sin(th) sin(phi)]*reshape(fpar(2:end),3,3)).^fpar(1),2);
predr = (1./wholesum).^(1/fpar(1));
[ccsurfaceG(1) ccsurfaceG(2) ccsurfaceG(3)] = sph2cart(th,phi,predr);
% Need to find the point on the surface in the ccdirB direction
[th,phi,r] = cart2sph(ccdirB(1),ccdirB(2),ccdirB(3));
wholesum = sum(abs([cos(phi).*cos(th) cos(phi).*sin(th) sin(phi)]*reshape(fpar(2:end),3,3)).^fpar(1),2);
predr = (1./wholesum).^(1/fpar(1));
[ccsurfaceB(1) ccsurfaceB(2) ccsurfaceB(3)] = sph2cart(th,phi,predr);
    
h1 = plot3(ccsurfaceG(1)*[1 scale],ccsurfaceG(2)*[1 scale],ccsurfaceG(3)*[1 scale],'k-');
set(h1,'Color',[0 .5 0],'Linewidth',2);
h2 = plot3(-ccsurfaceG(1)*[1 scale],-ccsurfaceG(2)*[1 scale],-ccsurfaceG(3)*[1 scale],'k-');
set(h2,'Color',[0 .5 0],'Linewidth',2);

h3 = plot3(ccsurfaceB(1)*[1 scale],ccsurfaceB(2)*[1 scale],ccsurfaceB(3)*[1 scale],'k-');
set(h3,'Color',[0 0 1],'Linewidth',2);
h4 = plot3(-ccsurfaceB(1)*[1 scale],-ccsurfaceB(2)*[1 scale],-ccsurfaceB(3)*[1 scale],'k-');
set(h4,'Color',[0 0 1],'Linewidth',2);
xlabel('\DeltaL/L');
ylabel('\DeltaM/M');
zlabel('\DeltaS/S');
set(gca,'View', [42 -74]);
set(gca,'Xlim',[-.3 .3],'Ylim',[-.3 .3],'Zlim',[-1 1]);
set(gca,'XTick',[-.3 0 .3],'YTick',[-.3 0 .3],'ZTick',[-1 0 1]);


%%
% Section 13 
% Orientation of monkey DTNT 25 Hz surfaces using a couple of different
% sets of fundamentals.
load T_xyz1931;
Vlambda = T_xyz1931(2,:);

AXISWIDTH = 2;
XMARGIN = .5;
YMARGIN = 2;

LCONELIM = 30; % plotting limits
SCONELIM = 100; % plotting limits
% Getting the gamut of the ViewPixx for plotting
load Dell4BitsCal.mat
load T_cones_smj10
load T_cones_synthgh2

PERCGAMUT = 90; % For plotting the surface fits, how near to gamut edge to get (in each direction)
dtntlistnames = {'NutDTNT25Hz.txt','FreyaDTNT25Hz.txt','SednaDTNT25Hz.txt'};
macpiglistnames = {'NutMacPig.txt','FreyaMacPig.txt','SednaMacPig.txt'};
subjectnames = {'Nut','Freya','Sedna'};
if (ismac)
    filelistpath = '/Volumes/NO BACKUP/NexFiles/nexfilelists/Greg/DTNT';
else
    filelistpath = 'N:\NexFiles\nexfilelists\Greg\DTNT';
end

data10deg = [];
datasynthfund = [];
for i = 1:length(dtntlistnames)
    filenames = fnamesFromTxt2([filelistpath,filesep,char(dtntlistnames{i})]);
    for fileidx = 1:size(filenames,1)
        if (strncmp(filenames{fileidx},'sf',2))
            continue;
        end
        stro = nex2stro(findfile(char(filenames{fileidx})));
        [thresholds, colordirs, sfs, QuestTrajectories] = DTquestUnpackGH(stro, 'mode');        
        Morig = reshape(stro.sum.exptParams.m_mtx,3,3);
        spds = reshape(stro.sum.exptParams.mon_spect, length(stro.sum.exptParams.mon_spect)/3,3);
        bkgndrgb = [stro.sum.exptParams.bkgnd_r stro.sum.exptParams.bkgnd_g stro.sum.exptParams.bkgnd_b];
        Mnew = T_cones_synthgh2'*spds;
        cc_original = colordirs.*repmat(thresholds,1,size(colordirs,2));
        cc_new = ConvertConeContrastBasis(Morig, Mnew, bkgndrgb, cc_original);
        synththresholds = sqrt(sum(cc_new.^2,2));
        synthcolordirs = mkbasis(cc_new')';
           
        ncoldir = length(thresholds);
        lmsmat = mkbasis(colordirs')';
        Loog = [];
        for j = 1:length(QuestTrajectories)
            Loog(j,1) = QuestTrajectories{j}(end) == max(QuestTrajectories{j});
        end
        data10deg = [data10deg; repmat(i,ncoldir,1) repmat(thresholds,1,3).*colordirs Loog];
        datasynthfund = [datasynthfund; repmat(i,ncoldir,1) repmat(synththresholds,1,3).*synthcolordirs Loog];
    end
end

% Getting some Macpig data
if (ismac)
    filelistpath = '/Volumes/NO BACKUP/NexFiles/nexfilelists/Greg/DTEM';
else
    filelistpath = 'N:\NexFiles\nexfilelists\Greg\DTEM';
end
macpigdata = [];
for i = 1:length(macpiglistnames)
    filenames = fnamesFromTxt2([filelistpath,filesep,char(macpiglistnames{i})]);
    for fileidx = 1:size(filenames,1)
        stro = nex2stro(findfile(char(filenames{fileidx})));
        if stro.sum.exptParams.rf_x == 50 % only get the 5 degree eccentric thresholds
            [thresh, colorDirs, sfs] = DTquestUnpack(stro, 'mode'); close(gcf);
            normcolordirs = colorDirs./repmat(sqrt(sum(colorDirs.^2,2)),1,3);
            guns = reshape(stro.sum.exptParams.mon_spect,81,3);
            bkgndrgb = [stro.sum.exptParams.bkgnd_r stro.sum.exptParams.bkgnd_g stro.sum.exptParams.bkgnd_b];
            M = reshape(stro.sum.exptParams.m_mtx,3,3);
            bkgndlms = M*bkgndrgb';
            bkgndlum = bkgndrgb*guns'*Vlambda';
            threshrgb = inv(M)*(repmat(bkgndlms',2,1).*(1+repmat(thresh/100,[1 3]).*normcolordirs([1 2],:)))';
            deltagun = threshrgb-repmat(bkgndrgb',1,2);
            % Quantifying thresholds in 2 deg luminance contrast units
            thresholds = [nan; nan];
            for j = 1:size(deltagun,2)
                thresholds(j) = (guns*deltagun(:,j))'*Vlambda';
                thresholds(j) = thresholds(j)./bkgndlum;
            end
            macpigdata = [macpigdata; i thresholds'.*100];
        end
    end
end

% Now doing the plotting
figure('Units','inches','PaperPosition',[0 0 8.5 11],'Position',[0 0 8.5 11],'PaperOrientation','Portrait','DefaultAxesUnits','inches');
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultTextFontName','Helvetica','DefaultTextFontSize',12);
colordef white;

Loog = data10deg(:,5);
Sconeweights = [];
for k = 1:2
    if (k ==1)
        data = datasynthfund;
        
    else
        data = data10deg;
    end
    for i = 1:length(subjectnames)
        axes('position',[1+(i-1)*(AXISWIDTH+XMARGIN),3+(k-1)*(AXISWIDTH+YMARGIN),AXISWIDTH, AXISWIDTH]);
        hold on;
        L = data(:,1) == i;
        
        for j = [-1, 1]
            h = plot3(j*data(L&~Loog,2),j*data(L&~Loog,3),j*data(L&~Loog,4),'ko');
            set(h,'MarkerFaceColor','black','MarkerSize',4,'MarkerEdgeColor','none');
            h = plot3([zeros(sum(L&Loog),1) j*data(L&Loog,2)]',[zeros(sum(L&Loog),1) j*data(L&Loog,3)]',[zeros(sum(L&Loog),1) j*data(L&Loog,4)]','k-');
            set(h(:),'Color',[.5 .5 .5]);
        end
        % Fitting data
        [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(data(L,[2 3 4]), Loog(L));
        coneweights = (xformmat*planeparams)';
        if (sign(coneweights(1)) == -1)
            coneweights = -coneweights;
        end
        title(num2str(coneweights,2));
        
        [xx yy zz] = meshgrid(linspace(-20,20,10),...
            linspace(-20,20,10),...
            linspace(-100,100,10));
        xformedxyz = [xx(:) yy(:) zz(:)];
        
        fr = abs(xformedxyz*coneweights');
        fv = isosurface(xx,yy,zz,reshape(fr,size(xx)), 1);
        h = patch(fv);
        
        set(h,'EdgeColor', 'none', 'FaceAlpha',.5,'FaceColor','green','Edgealpha',0);
     %   set(h,'FaceVertexCData',repmat([0 .75 0],size(fv.vertices,1),1))
     %   set(h,'CDataMapping','direct');
      % set(h,'FaceColor','interp','EdgeColor','none');
        lighting gouraud
        material metal
        camlight headlight
        set(gca,'View',[45 22]);
        set(gca,'xlim',[-50 50],'ylim',[-50 50]);
        set(gca,'Xtick',[-LCONELIM 0 LCONELIM],'Xticklabel',num2str([-LCONELIM 0 LCONELIM]'/100))
        xlabel('\DeltaL/L');
        set(gca,'Ytick',[-LCONELIM 0 LCONELIM],'Yticklabel',num2str([-LCONELIM 0 LCONELIM]'/100))
        ylabel('\DeltaM/M');
        set(gca,'Ztick',[-SCONELIM 0 SCONELIM],'Zticklabel',num2str([-SCONELIM 0 SCONELIM]'/100))
        zlabel('\DeltaS/S');
        Sconeweights(i,k) = coneweights(3)./sum(abs(coneweights));
    end
end
set(gcf,'Renderer','painters');
    

%%
% Section 13.1
% Bar plot of normalized S-cone weights (with errorbars?)
% requires the previous section be run
AXISWIDTH = 3;
xMARGIN = .5;
% Bootstrapping for SEs
se = [];
niter = 00;
for k = 1:2
    if (k ==1)
        data = datasynthfund;
    else
        data = data10deg;
    end
    tmpsconeweights = [];
    for i = 1:length(subjectnames)
        tmpsconeweights = [];
        L = find(data(:,1) == i);
        for j = 1:niter
            tmpL = unidrnd(size(L,1),size(L,1),1);
            [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(data(tmpL,[2 3 4]), Loog(tmpL));
            coneweights = (xformmat*planeparams)';
            tmpsconeweights(j) = coneweights(3)./sum(abs(coneweights))
        end
        se(i,k) = std(tmpsconeweights);
    end
end

bg = [];
for i = 1:macpigdata(end,1)
    L = macpigdata(:,1) == i;
    bg(i)=geomean(macpigdata(L,3)./ macpigdata(L,2));
end

% Now doing the plotting
BARWIDTH = .3;
BGOFFSET = .5;
figure('Units','inches','PaperPosition',[0 0 8.5 11],'Position',[0 0 8.5 11],'PaperOrientation','Portrait','DefaultAxesUnits','inches');
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultTextFontName','Helvetica','DefaultTextFontSize',12);
colordef white;

for i = 1:2
    axes('position',[3, 1+(i-1)*(AXISWIDTH+XMARGIN), AXISWIDTH, AXISWIDTH]);
    hold on
    h1 = bar([1:length(subjectnames)]+.2, Sconeweights(:,i)',BARWIDTH,'black');
    set(gca,'Ylim',[-.1 .1],'Xlim',[0 4]);
    %if (i==1)
        axes('position',[3, 1+(i-1)*(AXISWIDTH+XMARGIN), AXISWIDTH, AXISWIDTH]);
        h2 = bar([1:length(subjectnames)]-.2, bg-BGOFFSET,BARWIDTH,'white');
        set(gca,'color','none','Xlim',[0 4], 'Ylim',[-.4 .4],'YAxisLocation','right');
        ytick = get(gca,'Ytick'); set(gca,'Yticklabel',num2str(ytick'+BGOFFSET,3))
    %end
    if (niter > 0)
        errorbar([1:length(subjectnames)],Sconeweights(:,i),se(:,i),'k-','Linestyle','none');
    end
end

%%
% Section 14: Monitor spectra across intensities
figure;
plotcounter = 0;
for calname = {'Dell4BitsCal.mat' 'ViewPixx.mat'}
    cals = load(calname{:});
    cal = cals.cals{end};
    colors = [1 0 0; 0 .5 0; 0 0 .75];
    for channel = 2:3
        spds = reshape(cal.rawdata.mon(:,channel), cal.S_device(3), []);
        spds = spds.*cal.S_device(2);  % Are spds in 1 nm bins?
        [u,s,~] = svd(spds);
        
        subplot(4,4,1+4*plotcounter); hold on
        plot(diag(s(1:4,1:4))./sum(diag(s)), ['o' colors(channel,:)],'MarkerFaceColor',colors(channel,:),'MarkerEdgeColor',colors(channel,:))
        set(gca, 'xtick', [], 'Xlim',[0.5 4.5],'yscale', 'log','ylim',[0.001 1],'Ytick',[.001 .01 .1 1])
        
        subplot(4,4,2+4*plotcounter); hold on
        plot(linspace(380,780,cal.S_device(3)),u(:,2),'Color',colors(channel,:),'LineWidth',2);
        set(gca,'xlim',[380 780],'xtick',[380 580 780],'ylim',[min(u(:,2))-.05 max(u(:,2))+.05],'ytick',[])
       
        subplot(4,4,3+4*plotcounter); hold on
        plot(linspace(380,780,cal.S_device(3)),-u(:,1),'Color',colors(channel,:),'LineWidth',2);
        set(gca,'xlim',[380 780],'xtick',[380 580 780],'ylim',[min(-u(:,1))-.05 max(-u(:,1))+.05],'ytick',[])
 
        subplot(4,4,4+4*plotcounter); hold on
        h = plot(linspace(380,780,cal.S_device(3)),spds./repmat(max(spds),size(spds,1),1),'Color',colors(channel,:));
        set(gca,'xlim',[380 780],'xtick',[380 580 780]);
        
        for i = 1:length(h)
            set(h(i),'Color',colors(channel,:)*i/length(h))
        end
        plotcounter = plotcounter + 1;
    end
end

%%
% Section A: Are monkey B/G ratios related to human B/G ratios by a
% (vertical) scaling?

load T_xyz1931;
Vlambda = T_xyz1931(2,:);
filelistpath = 'C:\NO BACKUP\NexFiles\nexfilelists\Greg\DTEM';
listnames = {'GregMacPig.txt','ZackMacPig.txt','KaliMacPig.txt','SednaMacPig.txt'};
totaldata = {};
for i = 1:length(listnames)
    filenames = fnamesFromTxt2([filelistpath,'/',char(listnames{i})]);
    data = [];
    for fileidx = 1:size(filenames,1)
        stro = nex2stro(findfile(char(filenames{fileidx,:})));
        [thresh, colorDirs, sfs] = DTquestUnpack(stro, 'mode'); close(gcf);
        normcolordirs = colorDirs./repmat(sqrt(sum(colorDirs.^2,2)),1,3);
        guns = reshape(stro.sum.exptParams.mon_spect,81,3);
        bkgndrgb = [stro.sum.exptParams.bkgnd_r stro.sum.exptParams.bkgnd_g stro.sum.exptParams.bkgnd_b];
        M = reshape(stro.sum.exptParams.m_mtx,3,3);
        bkgndlms = M*bkgndrgb';
        bkgndlum = bkgndrgb*guns'*Vlambda';
        threshrgb = inv(M)*(repmat(bkgndlms',2,1).*(1+repmat(thresh/100,[1 3]).*normcolordirs([1 2],:)))';
        thresholds = (threshrgb'*guns'*Vlambda')./bkgndlum;
        ecc = stro.sum.exptParams.rf_x;
        data = [data; ecc/10 thresholds'];
    end
    totaldata{i} = data;
end

cols = [1 0 0; 0 1 0; 0 0 1; 0 0 0];
figure; axes; hold on;
x = linspace(0,8,50);
for i =1:length(totaldata)
    tmp = totaldata{i};
    if (i>=3)
        tmp(:,3) = tmp(:,3)*1.013;
    end
    pp_g  = csaps(tmp(:,1),log10((tmp(:,2)-1)*100),.5);
    pp_b  = csaps(tmp(:,1),log10((tmp(:,3)-1)*100),.5);
    h = plot(x,10.^fnval(x,pp_b)./10.^fnval(x,pp_g),'k-','LineWidth',2)
    set(h,'Color',cols(i,:));
    h = plot(tmp(:,1),(tmp(:,3)-1)./(tmp(:,2)-1),'ko');
    set(h,'MarkerFaceColor',cols(i,:));
end

%%
% Section B
% Looking at how sensitive the L/M cone ratio is when measured with red vs
% green phosphors or blue vs green.
% Bottom line: Neither blue vs green nore red vs green is superior for
% esimating cone weights/isoluminant points.

load Dell4BitsCal.mat
load T_cones_smj10
mon_spd = SplineRaw([380:2:780]', cals{end}.P_device, [380:5:780]');

M = T_cones_smj10*mon_spd;
bkgndrgb = [.5 .5 .5];
bkgndlms = M*bkgndrgb';
% Looking the response of a linear L+M mechanism to different modulation
% depths of out of phase red and green
LMweights = [.2 .2];
level = [0:.1:1];
for i = 1:length(level)
    stim1 = bkgndrgb;
    stim2 = [1 level(i) bkgndrgb(3)];
    stim3 = [bkgndrgb(1) level(i) 1];
    resp1 = LMweights*M([1 2],:)*stim1';
    resp2 = LMweights*M([1 2],:)*stim2';
    resp3 = LMweights*M([1 2],:)*stim3';
    rgmodulation(i) = resp1-resp2;
    bgmodulation(i) = resp1-resp3;    
end
figure; axes; hold on;
plot(level,rgmodulation,'r-');
plot(level,bgmodulation,'b-');
xlabel('green level vs red = 1 or blue = 1');


%%
% Section C
% Computing the luminance of the ViewPixx monitor
% We know from direct measurements that the photopic ViewPixx background is
% ~50 cd/m^2 and the scotopic background on the Dell 4 should be about ~8
% (photopic) cd/m^2

load ViewPixx
load T_xyz1931
cal = cals{end};
RGB = cal.bgColor;
rgb = [];
for i = 1:3
    rgb(i) = cal.gammaTable(round(RGB(i)*256)-1);
end
spect = cal.P_device*rgb';
spect1 = SplineSpd([380:2:780]',spect,[380:5:780]');
lum1 = T_xyz1931(2,:)*spect1*683 % according to Mike at Photoresearch, this should be in cd/m^2

splinedVlambda = SplineRaw([380:5:780]',T_xyz1931(2,:)',[380:2:780]');
lum2 = splinedVlambda'*spect*683
% For some reason I have to multiply these by 5/2 (doesn't matter which I
% spline - spd or vlambda)
lum1*5/2
lum2*5/2
%
load Dell4blackbkgnd
cal = cals{end};
spect = cal.P_ambient;
plot(spect);
lum = splinedVlambda'*spect*683;
lum*5/2

% getting scotopic luminance of scotopic background
load('wrattennd1.mat');
load T_rods;
NDtransmittance = wratten_nd1(:,2);
filteredspect = SplineSpd([380:2:780]',spect,[380:5:780]') .* NDtransmittance .^6;
scotlum = 1700*T_rods*filteredspect*5/2;% The factor of 5/2 comes from the fact that V'(lambda) is measured on a 5 nm lattice and 
% measurements from the PR705 are on a 2 nm lattice. I don't entirely get
% this but it works (see Section C).

%%
% Section D
% den_mac_bone agrees with CVRL database. 390:1:830
% den_lens_ssf agrees with CVRL database. 390:1:830

load den_mac_bone
load den_lens_ssf
load T_cones_ss10
load T_cones_ss2
load T_xyz1931;
Vlambda = [1.98 1 0]*T_cones_ss2;
Vlambda = Vlambda./max(Vlambda);

macpigtransmittance = 1./(10.^(den_mac_bone./max(den_mac_bone)*0.095));
lenstransmittance = 1./(10.^(den_lens_ssf./1.76)*1.63);
lenstransmittance = 1./(10.^den_lens_ssf);

absorptance = T_cones_ss10'./repmat(macpigtransmittance,1,3)./repmat(lenstransmittance,1,3);

% s-cone peaks at 424 normally
shift = 435-424;
interpvals = interp1([390:830],absorptance(:,3),[390-shift:390-1],'linear','extrap')
tmp = [interpvals'; absorptance(1:end-shift,3)];
fundsnopig = [absorptance(:,[1 2]) tmp];
funds = fundsnopig.*repmat(macpigtransmittance,1,3).*repmat(lenstransmittance,1,3);
funds = funds./repmat(max(funds),size(funds,1),1);

figure; axes; hold on;
plot([390:830],funds);
plot([390:830],T_cones_ss10','k-');

load Dell4BitsCal
cal = cals{end};
spds = [];
for i = 1:size(cal.P_device,2)
    spds(:,i) = SplineRaw([380:2:780]',cal.P_device(:,i),[390:1:830]')
end

% predictions with new fundamentals
M = funds'*spds;
glum = Vlambda*spds(:,2);
blum = Vlambda*spds(:,3);
coneweights = [1 1 .18];
gpred = coneweights./sum(abs(coneweights))*M(:,2);
bpred = coneweights./sum(abs(coneweights))*M(:,3);
newpred = (gpred./bpred)./(glum./blum)

coneweights./sum(abs(coneweights))

% Prediction with old fundamentals
M = T_cones_ss10*spds;
glum = Vlambda*spds(:,2);
blum = Vlambda*spds(:,3);
coneweights = [1 1 .19];
gpred = coneweights./sum(abs(coneweights))*M(:,2);
bpred = coneweights./sum(abs(coneweights))*M(:,3);
oldpred = (gpred./bpred)./(glum./blum)
% Trying to hit .5554
coneweights./sum(abs(coneweights))




%%
% Section D cont.d
% OK, now reanalyzing the S-cone weight to luminance-tuned V1 neurons using
% this new set of fundamentals (which got from 390:830, remember)
% Pretty much all of this code was take from section 4.1.
% Splining the monitor SPDs instead of the fundamentals.

load den_lens_ssf;
macpigtransmittance = 1./(10.^(den_mac_bone./max(den_mac_bone)*0.095));
lenstransmittance = 1./(10.^(den_lens_ssf));
% fundsnopig taken from above
fundsnolens = fundsnopig.*repmat(macpigtransmittance,1,3).*repmat(lenstransmittance,1,3);
fundsnolens = fundsnolens./repmat(max(fundsnolens),size(fundsnolens,1),1);

% First gathering all the data
data = [];
if(ispc)
    filelistpath = 'N:\NexFiles\nexfilelists\Greg\NT';
else
   filelistpath = '/Volumes/NO BACKUP/NexFiles/nexfilelists/Greg/NT';
end
[fnames, spikeIdx] = fnamesFromTxt2([filelistpath, filesep, 'NTlum2.txt']);
for cellcounter = 1:size(fnames,1)
    filename = findfile(char(fnames{cellcounter}));
    paradigmID = getparadigmID(filename);
    if paradigmID == 103
        NT = nex2stro(filename);
    end
    out = NTpreprocess(NT,0,Inf);
    fundamentals = NT.sum.exptParams.fundamentals;
    fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);

    mon_spd = NT.sum.exptParams.mon_spd;
    mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
    mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
    M = fundamentals'*mon_spd;
    Loog = logical(out(:,7));
    data(cellcounter).M = M;
    data(cellcounter).data = out;
    data(cellcounter).bkgndrgb = NT.sum.exptParams.bkgndrgb;
    data(cellcounter).monspd = mon_spd;
end

% Now, computing cone weights as a function of lens density
lensat400nm = linspace(0,2,20);
for i = 1:length(data)
    i
    cc_original = data(i).data(:,[2:4]).*repmat(data(i).data(:,5),1,3);
    Loog = data(i).data(:,end);
    allconeweights = [];
    for j = 1:length(lensat400nm)
        j
        idx = find(390:830 == 400);
        lenstransmittance = 1./(10.^(den_lens_ssf./den_lens_ssf(idx).*lensat400nm(j)));
        tmpfund = fundsnolens.*repmat(lenstransmittance,1,3);
        tmpfund = tmpfund./repmat(max(tmpfund),size(tmpfund,1),1);
        mon_spd = SplineRaw([380:5:780]', data(cellcounter).monspd, [390:1:830]');
        M = tmpfund'*mon_spd;
        scaled = ConvertConeContrastBasis(data(i).M, M, data(i).bkgndrgb, cc_original);
        [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(scaled, Loog);
        coneweights = (xformmat*planeparams)';
        if (sign(coneweights(2)) == -1)
            coneweights = -coneweights;
        end
        allconeweights = [allconeweights; coneweights];
    end
    data(i).coneweights = allconeweights;
end

% Doing the plotting
% At what lens density do we have zero S-cone weight?

figure('Units','inches','PaperPosition',[0 0 8.5 11],'Position',[0 0 8.5 11],'PaperOrientation','Portrait','DefaultAxesUnits','inches');
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultTextFontName','Helvetica','DefaultTextFontSize',12);
colordef white;

for whichplot = 1:2  % S-cone or L-M
    axes('position',[2 6.5-5*(whichplot-1) 3 3]); hold on;

    % Computing means and standard deviations
    allconeweights = cat(3,data(:).coneweights);
    mn = []; stdev = [];
    for j = 1:length(lensat400nm)
        tmpconeweights = squeeze(allconeweights(j,:,:));
        normconeweights = tmpconeweights./repmat(sum(abs(tmpconeweights)),3,1);
        if whichplot ==1
            mn(j) = mean(normconeweights(3,:));
            stdev(j) = std(normconeweights(3,:));
        else
            if (any(any(tmpconeweights([1 2],:) < 0)))
                mn(j) = nan;
                stdev(j) = nan;
            else
                mn(j) = mean(log10(tmpconeweights(1,:)./tmpconeweights(2,:)));
                stdev(j) = std(log10(tmpconeweights(1,:)./tmpconeweights(2,:)));
            end
        end
    end

    tcrit = tinv(.025,size(normconeweights,2)-1);
    CI_lower = tcrit*stdev./sqrt(size(tmpconeweights,2));
    %h = errorbar(lensat400nm,mn, stdev./sqrt(size(normconeweights,2)));
    %set(h,'LineWidth',2,'Color','black');
    %h = plot(lensat400nm,mn,'ko','LineWidth',1,'MarkerFaceColor','black');
    x = linspace(lensat400nm(1),lensat400nm(end),length(mn));
    Lposwts = ~isnan(mn);
    x(x<min(lensat400nm(Lposwts))) = []; 
    x(x>max(lensat400nm(Lposwts))) = [];

    % Plotting a smooth curve for the mean
    if (whichplot == 1)
        plot(x,interp1(lensat400nm(Lposwts),mn(Lposwts),x,'spline'),'k-','LineWidth',2);
        densityatzeroweight = interp1(mn,lensat400nm,0,'spline');
        plot([lensat400nm(1) densityatzeroweight],[0 0],'k-');
        plot([densityatzeroweight densityatzeroweight],[0 -.05],'k-');
    else
        plot(x,interp1(lensat400nm(Lposwts),10.^mn(Lposwts),x,'spline'),'k-','LineWidth',2);
        densityatzeroweight = interp1(10.^mn(Lposwts),lensat400nm(Lposwts),1,'spline');
        plot([lensat400nm(1) densityatzeroweight],[1 1],'k-');
        plot([densityatzeroweight densityatzeroweight],[1 .25],'k-');
    end
    set(gca,'XLim',[min(lensat400nm) max(lensat400nm)]);

    % Plotting confidence intervals
    if (whichplot == 1)
        % Plotting a smooth curve for the confidence interval
        plot(x,interp1(lensat400nm,mn+CI_lower,x,'spline'),'k-','LineWidth',1);
        plot(x,interp1(lensat400nm,mn-CI_lower,x,'spline'),'k-','LineWidth',1);
        densityatCIweight(1) = interp1(mn-CI_lower,lensat400nm,0,'spline');
        densityatCIweight(2) = interp1(mn+CI_lower,lensat400nm,0,'spline');
        for j = 1:2
            plot([lensat400nm(j) densityatCIweight(j)],[0 0],'k:');
            plot([densityatCIweight(j) densityatCIweight(j)],[0 -.05],'k:');
        end
        set(gca,'YLim',[-.05 0.04]);
        ylabel('S-cone weight','FontSize',12);
    else % L/M  cone weight ratio
        plot(x,interp1(lensat400nm(Lposwts),10.^(mn(Lposwts) + CI_lower(Lposwts)),x,'spline'),'k-','LineWidth',1);
        plot(x,interp1(lensat400nm(Lposwts),10.^(mn(Lposwts) - CI_lower(Lposwts)),x,'spline'),'k-','LineWidth',1);
        densityatCIweight(1) = interp1(10.^(mn(Lposwts) + CI_lower(Lposwts)),lensat400nm(Lposwts), 1,'spline');
        densityatCIweight(2) = interp1(10.^(mn(Lposwts) - CI_lower(Lposwts)),lensat400nm(Lposwts), 1,'spline');
        for j = 1:2
            plot([lensat400nm(j) densityatCIweight(j)],[1 1],'k:');
            plot([densityatCIweight(j) densityatCIweight(j)],[1 .025],'k:');
        end
        set(gca,'Ylim',[.25 4]);
        set(gca,'Yscale','log');
        ylabel({'L-weight/M-weight'},'FontSize',14);
        set(gca,'Ytick',[.5 1 2],'Yticklabel',[.5 1 2])
    end
    xlabel('Lens density at 400 nm','FontSize',14);
end

% Sanity check
mon_spd = SplineRaw([380:5:780]', data(cellcounter).monspd, [390:1:830]');
M = fundsnolens'*mon_spd;
glum = Vlambda'*spds(:,2);
blum = Vlambda'*spds(:,3);
gpred = [1 1 .14]*M(:,2);
bpred = [1 1 .14]*M(:,3);
(gpred./bpred)./(glum./blum)


%%
% Section E
% What lens density is consistent with Wald's 1945 data (@ 400 nm)?
% Log sensitivity measurements
x = [365 405 436 492 546 578 621 691 713 750]; % wavelegnth
y = [-2.042 0.427 1.675 2.295 2.095 1.375 0.038 -3.635 -4.787 -5.890]; % rod sensitivity with lens
z = [-2.62 -1.05 -1.48 0 -1.78 -1.12 -3.87 -5.39 -6.52 -7.68]; % rod sensitivity without lens

% I think this z table has some values misplaced..

% Well that wasn't very useful

plot(x,10.^(z-y))

%%
% Section F
% Finding an L:M cone ratio and lens density that reconcile Freya andGH
% Apollo's photopic and scotopic B/G ratio data.

% For Apollo (photopic) blue was somewhat brighter than it was for Freya.

load T_xyz1931;
Vlambda = T_xyz1931(2,:);
load wrattennd1;
load T_rods;

if (ispc)
    filelistpath = 'C:\NO BACKUP\NexFiles\nexfilelists\Greg\DTEM';
else
    filelistpath = '/Volumes/NO BACKUP/NexFiles/nexfilelists/Greg/DTEM';
end
data = [];
monspect = [];
for PHOTOPIC = 0:1
    if (PHOTOPIC)
        listnames = {'GregMacPig.txt','LeahMacPig.txt','ZackMacPig.txt','SednaMacPig.txt','ApolloMacPig.txt','FreyaMacPig.txt'};
    else
        listnames = {'GregDTscot.txt','LeahDTscot.txt','ZackDTscot.txt','SednaDTscot.txt','ApolloDTscot.txt','FreyaDTscot.txt'};
    end
    for i = 1:length(listnames)
        filenames = fnamesFromTxt2([filelistpath,filesep,char(listnames{i})]);
        for fileidx = 1:size(filenames,1)
            stro = nex2stro(findfile(char(filenames{fileidx})));     
            if (PHOTOPIC)
                if stro.sum.exptParams.rf_x == 50 % only get the 5 degree eccentric thresholds
                    [thresh, colorDirs, sfs] = DTquestUnpack(stro, 'mode'); close(gcf);
                    normcolordirs = colorDirs./repmat(sqrt(sum(colorDirs.^2,2)),1,3);
                    guns = reshape(stro.sum.exptParams.mon_spect,81,3);
                    bkgndrgb = [stro.sum.exptParams.bkgnd_r stro.sum.exptParams.bkgnd_g stro.sum.exptParams.bkgnd_b];
                    M = reshape(stro.sum.exptParams.m_mtx,3,3);
                    bkgndlms = M*bkgndrgb';
                    bkgndlum = bkgndrgb*guns'*Vlambda';
                    threshrgb = inv(M)*(repmat(bkgndlms',2,1).*(1+repmat(thresh/100,[1 3]).*normcolordirs([1 2],:)))';
                    deltagun = threshrgb-repmat(bkgndrgb',1,2);
                    % Quantifying thresholds in 2 deg luminance contrast units
                    thresholds = [nan; nan];
                    for j = 1:size(deltagun,2)
                        thresholds(j) = (guns*deltagun(:,j))'*Vlambda';
                        thresholds(j) = thresholds(j)./bkgndlum;
                    end
                    %thresholds = (threshrgb'*guns'*Vlambda')./bkgndlum;
                    data = [data; PHOTOPIC i thresholds'];
                    if length(stro.sum.exptParams.mon_spect) == 303
                        mon_spd = SplineSpd([380:4:780]', reshape(stro.sum.exptParams.mon_spect,101,3), [380:5:780]');
                    elseif length(stro.sum.exptParams.mon_spect) == 603
                        mon_spd = SplineSpd([380:2:780]', reshape(stro.sum.exptParams.mon_spect,201,3), [380:5:780]');
                    else
                       mon_spd =  reshape(stro.sum.exptParams.mon_spect,81,3);
                    end
                    monspect = [monspect, mon_spd(:)];
                end
                
            else % SCOTOPIC, not PHOTOPIC
                r_idx = strcmp(stro.sum.trialFields(1,:), 'stim_r');
                g_idx = strcmp(stro.sum.trialFields(1,:), 'stim_g');
                b_idx = strcmp(stro.sum.trialFields(1,:), 'stim_b');
                rgbs = stro.trial(:,r_idx|g_idx|b_idx);
                questmode_idx = strcmp(stro.sum.trialFields(1,:), 'quest_mode');
                thresholds = nan*ones(2,1);
                if length(stro.sum.exptParams.mon_spd) == 303
                    mon_spd = SplineRaw([380:4:780]', reshape(stro.sum.exptParams.mon_spd,101,3), [380:5:780]');
                else
                    mon_spd = SplineRaw([380:2:780]', reshape(stro.sum.exptParams.mon_spd,201,3), [380:5:780]');
                end
                new_mon_spd = mon_spd.*repmat(wratten_nd1(:,2).^6,1,3);
                scotlum = 1700*T_rods*new_mon_spd*5/2;% The factor of 5/2 comes from the fact that V'(lambda) is measured on a 5 nm lattice and
                monspect = [monspect, new_mon_spd(:)];
 
                for gun = 2:3
                    L = logical(rgbs(:,gun));
                    if (any(L))
                        all_modes = stro.trial(L, questmode_idx);
                        thresholds(gun) = scotlum(gun).*all_modes(end);
                    end
                end
                data = [data; PHOTOPIC i thresholds(2:3)'];
              
            end
        end
    end
end

% Now getting the means
figure; axes; hold on;
symbols = {'ks','k^','ko','bs','b^','bo'}; % Apollo: black squares, Freya: magenta triangles, Sedna: blue circles
for monkey = 1:data(end,2)
    Lscot = data(:,1) == 0;
    Lphot = data(:,1) == 1;
    Lmonk = data(:,2) == monkey;
    
    tmp = [10.^mean(log10(data(Lmonk&Lscot,4)./data(Lmonk&Lscot,3))),...
        10.^mean(log10(data(Lmonk&Lphot,4)./data(Lmonk&Lphot,3)))]
    
    plot(tmp(1),tmp(2),symbols{monkey});
    xlabel('Scotopic B/G');
    ylabel('Photopic B/G');  
end

% stuff for photopic predictions
lensdensities = logspace(log10(1/2),log10(5),10);
LvMs = logspace(log10(1/100),log10(100),20);
load('T_cones_myss10'); load den_lens_ss;  load den_mac_ss;
macpigtransmittance = 1./(10.^(den_mac_ss./max(den_mac_ss)*0.095)); % might need to add this back?
lenstransmittance = 1./(10.^(den_lens_ss));
absorptance = T_cones_ss10./repmat(macpigtransmittance,1,3)./repmat(lenstransmittance,1,3);
synthfunds = absorptance./repmat(max(absorptance),81,1);

% stuff for scotopic predictions
wls = [398 420 440 459 479 501 520 539 560 579 599 619 640 659 681 700 720];
s = [-.697 -.436 -.314 -.146 -.017 0 -.151 -.356 -.749 -1.22 -1.755 -2.312 -3.093 -3.743 -4.503 -5.147 -5.657];
opticaldensity = .3;
rodactionspectra = interp1(wls,s,[380:5:780],'linear','extrap');

synthpred = [];
for i = 1:length(lensdensities)
    for j = 1:length(LvMs)
        lenstransmittance = 1./(10.^(den_lens_ss./den_lens_ss(5)*lensdensities(i)));
        tmpfunds = synthfunds.*repmat(lenstransmittance,1,3);
        tmpfunds = tmpfunds./repmat(max(tmpfunds),81,1);
        [u,s,v] = svd(monspect(:,Lphot));
        spd = sign(max(u(:,1)))*reshape(u(:,1),81,3);
        M = tmpfunds'*spd;
        LVLAMBDA = LvMs(j);
        gpred = [LVLAMBDA 1]*M([1 2],2);
        bpred = [LVLAMBDA 1]*M([1 2],3);
        glum = Vlambda*spd(:,2);
        blum = Vlambda*spd(:,3);
        synthpred(i,j,1) = (gpred./bpred)./(glum./blum);

        rodfund = 1-10.^(-(10.^rodactionspectra').*opticaldensity);
        rodfund = rodfund.*lenstransmittance.*macpigtransmittance;
        rodfund = rodfund./repmat(max(rodfund),81,1);
        [u,s,v] = svd(monspect(:,Lscot));
        spd = sign(max(u(:,1)))*reshape(u(:,1),81,3);
   
        scotlum_h = 1700*T_rods*spd(:,[2 3]);
        scotlum_m = 1700*rodfund'*spd(:,[2 3]);
        synthpred(i,j,2) = (scotlum_m(1)./scotlum_m(2))./(scotlum_h(1)./scotlum_h(2));
    end
end
plot(squeeze(synthpred(:,:,2)),squeeze(synthpred(:,:,1)),'.');



%%
% Section G
% Comparing monitor phosphor spectra across data files.

if (ispc)
    filelistpath = 'C:\NO BACKUP\NexFiles\nexfilelists\Greg\DTEM';
else
    filelistpath = '/Volumes/NO BACKUP/NexFiles/nexfilelists/Greg/DTEM';
end
listnames = {'GregMacPig.txt','ApolloMacPig.txt','FreyaMacPig.txt','SednaMacPig.txt'};
listnames = {'GregDTscot.txt','ApolloDTscot.txt','FreyaDTscot.txt','SednaDTscot.txt'};

monspect = [];
data = [];
for i = 1:length(listnames)
    filenames = fnamesFromTxt2([filelistpath,filesep,char(listnames{i})]);
    for fileidx = 1:size(filenames,1)
        stro = nex2stro(findfile(char(filenames{fileidx,:})));
        
        if (isfield(stro.sum.exptParams,'mon_spect'))
            tmpspect = stro.sum.exptParams.mon_spect;
        elseif (isfield(stro.sum.exptParams,'mon_spd'))
            tmpspect = stro.sum.exptParams.mon_spd;
        else
            error('Cannot find phosphor spectra');
        end
        
        if length(tmpspect) == 303
            tmpspect = SplineSpd([380:4:780]', reshape(tmpspect,101,3), [380:5:780]');
        elseif length(tmpspect) == 603
            tmpspect = SplineSpd([380:2:780]', reshape(tmpspect,201,3), [380:5:780]');
        else
            tmpspect =  reshape(tmpspect,81,3);
        end
        monspect = [monspect, tmpspect(:)];
        data = [data; i fileidx];
    end
end

figure; axes; hold on;
colors = colormap(lines(6))
for i = 1:size(data,1)
    plot(monspect(:,i)+data(i,1)/10000,'Color',colors(data(i,1),:),'linewidth',2);
end

% The DTscotopic monitor spectra all look about the same except that the
% monitor was recalibrated once in the middle of data collection from Greg.


%%
% Section H
% Does estimating the lens density separately for each observer (Leah,
% Zack, Greg, Apollo, Sedna, Freya) reduce the variance in photopic B/G?
% Code largely taken from section 1.1. estimated lens densities taken from
% section 9.1?

if (ispc)
    filelistpath = 'C:\NO BACKUP\NexFiles\nexfilelists\Greg\DTEM';
else
    filelistpath = '/Volumes/NO BACKUP/NexFiles/nexfilelists/Greg/DTEM';
end
load Vlambda_star;
Vlambda_star = Vlambda_star(:,2);

% making a version of Vlambda star that does not have any lens density
load den_lens_ss;
lens = den_lens_ss;
lenstransmittance = 1./(10.^(den_lens_ss));
fundnopig = T_cones_ss2./repmat(lenstransmittance,1,3);
Vlambda_star_nolens = 1.98*fundnopig(:,1)+fundnopig(:,2);  % 1.98L+ 1M=vlambda* GH 2/17/14 
Vlambda_star_nolens = Vlambda_star_nolens./max(Vlambda_star_nolens);

lenscorr400 = [5.2197 2.5451 3.4103 1.4023 1.3063 2.4779];
% lenscorr400 = [den_lens_ss(5) den_lens_ss(5) den_lens_ss(5)
% den_lens_ss(5) den_lens_ss(5) den_lens_ss(5)]; % for debugging

listnames = {'GregMacPig.txt','ZackMacPig.txt','LeahMacPig.txt','ApolloMacPig.txt','FreyaMacPig.txt','SednaMacPig.txt'};
bgdata = [];
for i = 1:length(listnames)
    filenames = fnamesFromTxt2([filelistpath,filesep,char(listnames{i})]);        
    lenstransmittance = 1./(10.^(den_lens_ss./den_lens_ss(5).*lenscorr400(i)));
    Vlambda_corrected = Vlambda_star_nolens.*lenstransmittance;
    Vlambda_corrected = Vlambda_corrected./max(Vlambda_corrected);
    data = [];
    for fileidx = 1:size(filenames,1)
        stro = nex2stro(findfile(char(filenames{fileidx,:})));
        [thresh, colorDirs, sfs] = DTquestUnpack(stro, 'mode'); close(gcf);
        normcolordirs = colorDirs./repmat(sqrt(sum(colorDirs.^2,2)),1,3);
        guns = reshape(stro.sum.exptParams.mon_spect,81,3);
        bkgndrgb = [stro.sum.exptParams.bkgnd_r stro.sum.exptParams.bkgnd_g stro.sum.exptParams.bkgnd_b];
        M = reshape(stro.sum.exptParams.m_mtx,3,3);
        bkgndlms = M*bkgndrgb';
        threshrgb = inv(M)*(repmat(bkgndlms',2,1).*(1+repmat(thresh/100,[1 3]).*normcolordirs([1 2],:)))';
        deltagun = threshrgb-repmat(bkgndrgb',1,2);
        ecc = stro.sum.exptParams.rf_x;
        
        % First calculating thresholds wrt vlambda star
        thresholds = [nan; nan; nan; nan];
        bkgndlum = bkgndrgb*guns'*Vlambda_star;
        for j = 1:size(deltagun,2)
            thresholds(j) = ((guns*deltagun(:,j))'*Vlambda_star)./bkgndlum;
            thresholds(j+2) = ((guns*deltagun(:,j))'*Vlambda_corrected)./bkgndlum;
        end
        data = [data; ecc/10 thresholds'];
    end
    
    % Collecting data for stats
    bg = data(data(:,1) == 7,[2 3 4 5]); % only looking at 7 deg eccentric
    bgdata = [bgdata; mean(log10(bg(:,2)./bg(:,1))) mean(log10(bg(:,4)./bg(:,3)))]
end

std(bgdata(:,2))./std(bgdata(:,1))
% looks like using customized lens densities does reduce the observer to
% observer variability in B/G
bins = [.5:.025:1];
figure;
for i = 1:2
    subplot(2,1,i); hold on;
    [n,x] = hist(10.^bgdata(:,i),bins);
    bar(x,n,'facecolor','black');
    [n,x] = hist(10.^bgdata([1 2 3],i),bins);
    bar(x,n,'facecolor','red');
end
std(10.^bgdata(:,2))./std(10.^bgdata(:,1))

bins = [-.3:.05:0];
figure;
titles = {'V*(\lambda)','Individually lens-adjusted luminous efficiency function'};
for i = 1:2
    subplot(2,1,i); hold on;
    [n,x] = hist(bgdata(:,i),bins);
    bar(x,n,'facecolor','black');
    [n,x] = hist(bgdata([1 2 3],i),bins);
    bar(x,n,'facecolor','red');
    title(titles(i));
    xlabel('log_1_0(B/G)'); ylabel('count');
end
legend('Monkeys','Humans');
std(bgdata(:,2))./std(bgdata(:,1))

%%
% Section I
% A big model looking at a variety of factors that could
% contribute to the photopic (and scotopic?) B/G with the approximately
% measured values. Factors are:
%     1) L:M ratio
%     2) S-cone contribution
%     3) Rod contribution
%     4) Lens optical density
%     5) photopigment optical density
%     6) Lambda max of S-cone (430 nm) (Not dealing with this)
% Code largely lifted from TenDegFundStuff, section 7.2
% 3/12/13 GDLH Fixing a bug. Using the Stockman and Sharpe nomogram instead
% of the "action spectra" that I was computing in TenDegFundStuff which was
% wrong.

load Dell4BitsCal
cal = cals{end};
spds = SplineSpd(cal.S_device, cal.P_device, [380 5 81]);

load den_lens_ss; 
load den_mac_ss;
load T_cones_myss10;
load Vlambda_star;
Vlambda = Vlambda_star(:,2)';

% Preparing cone actionspectra
macpigtransmittance = 1./(10.^(den_mac_ss./max(den_mac_ss)*0.095));
lenstransmittance = 1./(10.^(den_lens_ss));
absorptance = T_cones_ss10./repmat(macpigtransmittance,1,3)./repmat(lenstransmittance,1,3);
absorptance = absorptance./repmat(max(absorptance),81,1);
%plot([380:5:780],absorptance);
opticaldensity = 0.38;
%actionspectra = -log10(1-absorptance*(1-10^-opticaldensity));  <-- WRONG!
coneactionspectra = EnergyToQuanta([380 5 81],StockmanSharpeNomogram([380 5 81],[558.9 530.3 420.7])');
coneactionspectra = coneactionspectra./repmat(max(coneactionspectra),size(coneactionspectra,1),1);
%plot([380:5:780],coneactionspectra,'.');

% Preparing a rod fundamental
wls = [398 420 440 459 479 501 520 539 560 579 599 619 640 659 681 700 720];
s = [-.697 -.436 -.314 -.146 -.017 0 -.151 -.356 -.749 -1.22 -1.755 -2.312 -3.093 -3.743 -4.503 -5.147 -5.657];
rodactionspectrum = interp1(wls,s,[380:5:780],'spline','extrap')';
rodactionspectrum1 = interp1(wls,s,[380:5:780],'linear','extrap')';
L = [380:5:780] > 720;
rodactionspectrum(L) = rodactionspectrum1(L);
opticaldensity = .35; % .4 is the estimate from Baylor et al. (1984)  0.2 is Rushton (1972)
rodabsorptance = 1-10.^(-(10.^rodactionspectrum).*opticaldensity);
rodabsorptance = rodabsorptance./max(rodabsorptance);

% Preparing parameters to test
nlevels = 15;
titles = {'mac pig','PP OD','lens','log(L:M)','S','rod'};
mp_s = linspace(0,.15,nlevels); % 1
ppod_s = linspace(.1,.7,nlevels); % 2
len_s = linspace(0,4,nlevels); % 3
lm_s =  log10(logspace(log10(1/4), log10(4),nlevels)); % 4
%scone_s = 0;
%rod_s = 0;
scone_s = linspace(0,.1,nlevels); % 5
rod_s = linspace(0,.1,nlevels); % 6

ns = [length(mp_s) length(ppod_s) length(len_s) length(lm_s) length(scone_s) length(rod_s)]
designmatrix = fullfact(ns);
glum = Vlambda*spds(:,2);
blum = Vlambda*spds(:,3);

BGdata = nan*ones(size(designmatrix,1),1);
for i = 1:size(designmatrix,1)
    i
    mp = mp_s(designmatrix(i,1));
    ppod = ppod_s(designmatrix(i,2));
    l = len_s(designmatrix(i,3));
    lm = lm_s(designmatrix(i,4));
    s = scone_s(designmatrix(i,5));
    r = rod_s(designmatrix(i,6));
    
    macpigtransmittance = 1./(10.^(den_mac_ss./max(den_mac_ss)*mp));
    lenstransmittance = 1./(10.^(den_lens_ss./den_lens_ss(5)*l));
    
    funds = 1-10.^(-coneactionspectra.*ppod); % no need to normalize here
    funds = funds.*repmat(macpigtransmittance,1,3).*repmat(lenstransmittance,1,3);
    funds = funds./repmat(max(funds),81,1);
    sens = 10.^lm*funds(:,1)+funds(:,2);
    sens = sens./max(sens); % So everything peaks at 1 before we scale and add S-cone and rod funds.
    sens = sens+s*funds(:,3) + r*rodabsorptance;
    sens = sens./max(sens);
    gpred = sens'*spds(:,2);
    bpred = sens'*spds(:,3);
    
    BGdata(i) = (gpred./bpred)./(glum./blum);
end

L = BGdata > .59 & BGdata < .65;
acc = designmatrix(L,:);
[mp_s(acc(:,1))' ppod_s(acc(:,2))' len_s(acc(:,3))' 10.^lm_s(acc(:,4))' scone_s(1,acc(:,5))' rod_s(1,acc(:,6))'];
sum(L)./size(designmatrix,1)  % percent 

acc = designmatrix(L,:);
goodmodels = [mp_s(acc(:,1))' ppod_s(acc(:,2))' len_s(acc(:,3))' lm_s(acc(:,4))' scone_s(1,acc(:,5))' rod_s(1,acc(:,6))'];

% Plotting
figure('Units','inches','PaperPosition',[0 0 8.5 11],'Position',[0 0 8.5 11],'PaperOrientation','Portrait','DefaultAxesUnits','inches');
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',10);
set(gcf,'DefaultTextFontName','Helvetica','DefaultTextFontSize',12);
colordef white;
AXESWIDTH = 1.75;
MARGIN =.5;
XAXSTART = 1;
hax = []; bins = {};
for i = 1:6
    hax(i) = axes('position',[XAXSTART+mod(i-1,3)*(AXESWIDTH+MARGIN) 6+(AXESWIDTH+1.5*MARGIN)*(ceil(i/3)-1) AXESWIDTH AXESWIDTH]);
    bins{i} = unique(goodmodels(:,i));
    binwidth = (bins{i}(2)- bins{i}(1));
    [n,x] = hist(goodmodels(:,i), bins{i});
    if (~isempty(n))
        h = bar(x+binwidth*.2,n./sum(n),.5); % Normalized density
        set(h,'FaceColor','black','EdgeColor','none');
        xlabel(titles{i});
        xlims = [min(goodmodels(:,i)) max(goodmodels(:,i))];
        xlims = xlims + [-1 1]*binwidth/2;
    end
    set(gca,'xlim',xlims,'ylim',[0 max(n./sum(n))*1.1],'Ytick',[]);
end

% Correlation matrix. This stuff taken from prettycorr.m.
axes('position',[XAXSTART 1.5 2*AXESWIDTH 2*AXESWIDTH])
[r,p] = corr(goodmodels,'type','Spearman');
him = image((r+1)*32); colormap(gray(64));
for i = 1:size(goodmodels,2)
    for j = 1:size(goodmodels,2)
        h = text(i,j,[num2str(r(i,j),2)],'HorizontalAlignment','center','VerticalAlignment','middle','fontsize',7);
        if (r(i,j)<0)
            set(h,'color',[1 1 1]);
        end
        set(gca,'XTick',1:size(goodmodels,2),'XTickLabel',titles);
        set(gca,'YTick',1:size(goodmodels,2),'YTickLabel',titles);
        axis tight
    end
end

L = L & designmatrix(:,1) < 6; % only reasonable MP
sum(L)./size(designmatrix,1)

L = L & designmatrix(:,5) < 4; % only reasonable S-cone contributions
sum(L)./size(designmatrix,1)

L = L & designmatrix(:,6) < 4; % only reasonable rod contributions
sum(L)./size(designmatrix,1)

L = L & designmatrix(:,5) < 2; % zero S-cone contributions
sum(L)./size(designmatrix,1)

L = L & designmatrix(:,6) < 2; % zero rod contributions
sum(L)./size(designmatrix,1)

acc = designmatrix(L,:);
goodmodels = [mp_s(acc(:,1))' ppod_s(acc(:,2))' len_s(acc(:,3))' lm_s(acc(:,4))' scone_s(1,acc(:,5))' rod_s(1,acc(:,6))'];

for i = 1:6
    axes(hax(i)); hold on;
    [n,x] = hist(goodmodels(:,i),bins{i});
    binwidth = (bins{i}(2)- bins{i}(1));
    h = bar(x-binwidth*.2,n./sum(n),.5); % Normalized density
    set(h,'FaceColor','red','EdgeColor','none');
    if (n(1) ~= sum(L))
        set(gca,'ylim',[0 max(n./sum(n))*1.1]);
    end
end

% Small correlation matrix
% This stuff taken from prettycorr.m
ROWSIZE = 2*AXESWIDTH/size(goodmodels,2);
k = find(var(goodmodels) > 0);
[r,p] = corr(goodmodels(:,k),'type','Spearman');
axes('position',[5 1.5+2*ROWSIZE length(k)*ROWSIZE*[1 1]])
him = image((r+1)*32);
for i = k
    for j = k
        h = text(i,j,[num2str(r(i,j),2)],'HorizontalAlignment','center','VerticalAlignment','middle','fontsize',7);
        if (r(i,j)<0)
            set(h,'color',[1 1 1]);
        end
        set(gca,'XTick',1:size(goodmodels,2),'XTickLabel',titles);
        set(gca,'YTick',[]);
        axis tight
    end
end
set(gca,'YColor',[1 0 0],'XColor',[1 0 0]);

% Getting mean and SD of lens distribution
mean(len_s(acc(:,3)'))
std(len_s(acc(:,3)'))
% PCA
% goodmodels_z = (goodmodels-repmat(mean(goodmodels),[size(goodmodels,1),1]))./repmat(std(goodmodels),[size(goodmodels,1),1]);
% goodmodels_z(:,isnan(goodmodels_z(1,:))) = [];
% [coeff, score, latent, tsquare] = princomp(goodmodels_z);
% 
% figure;
% subplot(3,1,1);
% bar(coeff(:,1));
% set(gca,'XTickLabel',titles)
% subplot(3,1,2);
% bar(coeff(:,2));
% set(gca,'XTickLabel',titles)
% subplot(3,1,3);
% plot(latent,'ko','markerSize',10,'markerFaceColor','black');
% 


%%
% Section K
% Comparing equal energy white to illuminant D65
load spd_D65

lms_D65 = spd_D65'* T_cones_ss2;
lms_D65(1)./lms_D65(2)   % L/M with D65

lms_eew1 = sum(T_cones_ss2);
lms_eew1(1)/lms_eew1(2)   % L/M with *real* EEW

load dell4BitsCal
cal = cals{8};
RGB =  cal.bgColor;
rgb = [];
for i = 1:3
     rgb(i) = cal.gammaTable(round(RGB(i)*255)+1,i);
end
bkgnd = SplineSpd(cal.S_device, cal.P_device* rgb', [380 5 81]);
lms_bkgnd = bkgnd'* T_cones_ss2;
lms_bkgnd(1)./lms_bkgnd(2)   % L/M with our EE background
% About a 5% differences in L/M excitation with D65 and EEW

%%
% Section J: trying to get estimates of L:M from B/G
% Also trying to get lens density estimates?

if (ispc)
    filelistpath = 'C:\NO BACKUP\NexFiles\nexfilelists\Greg\DTEM';
else
    filelistpath = '/Volumes/NO BACKUP/NexFiles/nexfilelists/Greg/DTEM';
end

LMweights = [];
for HUMANS = 0:1
    if (HUMANS)
        listnames = {'GregMacPig.txt','ZackMacPig.txt','LeahMacPig.txt'};
    else
        listnames = {'ApolloMacPig.txt','FreyaMacPig.txt','SednaMacPig.txt','KaliMacPig.txt','NutMacPig.txt',};
    end
    
    for i = 1:length(listnames)
        filenames = fnamesFromTxt2([filelistpath,filesep,char(listnames{i})]);
        
        data = [];
        for fileidx = 1:size(filenames,1)
            stro = nex2stro(findfile(char(filenames{fileidx,:})));
            [thresh, colorDirs, sfs] = DTquestUnpack(stro, 'mode'); close(gcf);
            normcolordirs = colorDirs./repmat(sqrt(sum(colorDirs.^2,2)),1,3);
            threshlms = repmat(thresh/100,[1 3]).*normcolordirs([1 2],:); % Is this correct? In cc units?
            ecc = stro.sum.exptParams.rf_x;
            data = [data; ecc/10 threshlms(1,:) threshlms(2,:)];
        end
        L = data(:,1)>=5;
        LccG = data(L,2);
        MccG = data(L,3);
        LccB = data(L,5);
        MccB = data(L,6);
        figure; axes; hold on;
        plot(LccG,MccG,'o');
        plot(LccB,MccB,'o');
        plot([LccG'; LccB'],[MccG'; MccB'],'b-')
        pt1 = [geomean(LccG) geomean(MccG)];
        pt2 = [geomean(LccB) geomean(MccB)];
        plot(pt1(1),pt1(2),'*');
        plot(pt2(1),pt2(2),'*');
%        greg = [LccG MccG; LccB MccB];
%        LM = ones(size(greg,1),1)\greg;
%        LMweights = [LMweights; LM./norm(LM)]
       LMweights = [LMweights, inv([pt1; pt2])*[1;1]];

    end
end

%%
% Section K
% Looking at the Stockman nomogram fit to the S-cone photopigment
% absorbance spectrum.

a = -188862.970810906644;
b = 90228.966712600282;
c = -2483.531554344362;
d = -6675.007923501414;
e = 1813.525992411163;
f = -215.177888526334;
g = 12.487558618387;
h = -0.289541500599;
nm = [390:.1:500];
x = log10(nm)-log10(420.7/558);
y = a + b*x.^2 + c*x.^4 + d*x.^6 + e*x.^8 +f*x.^10 + g*x.^12 + h*x.^14;
%plot(x,y);
plot(nm,10.^y);

% Getting the S-cone photopigment absorbance function
load T_cones_myss10
load den_lens_ss
load den_mac_ss
%lenstransmittance = 1./(10.^den_lens_ss./den_lens_ss(5)*1.63);
macpigtransmittance = 1./(10.^(den_mac_ss./max(den_mac_ss)*0.095));
lenstransmittance = 1./(10.^den_lens_ss);

absorptance = T_cones_ss10./repmat(macpigtransmittance,1,3)./repmat(lenstransmittance,1,3);
absorptance = absorptance(:,3); % S-cone only
absorbance = -log10(1-absorptance/max(absorptance)*.2);  % I'm not sure about this caluculation
absorbance = absorbance./[380:5:780]' % Converting to *quantal* units, appropriate for absorbance
absorbance = absorbance./max(absorbance);

nm = [380:5:780];
figure; axes; hold on;
plot(nm,absorbance);
% 
% err = [];
% lmaxes = [415:434];
% for lmax = lmaxes
%     x = log10(nm)-log10(lmax/558);
%     y = a + b*x.^2 + c*x.^4 + d*x.^6 + e*x.^8 +f*x.^10 + g*x.^12 + h*x.^14;
%     plot(nm,10.^y,'g-');
%     tmp = (y-log10(absorbance')).^2;
%   %  tmp = (10.^y-absorbance').^2;
%    % tmp(16:end) = [];  % only going up to 455?
%     err = [err; sum(tmp)];
% end
% figure;
% plot(lmaxes,err,'.')

%%
% Section J
% Comparing isodet surface to isoresponse surface

load ('T_cones_smj10');
if (ismac)
    DTNTfilelistpath = '/Volumes/NO BACKUP/NexFiles/nexfilelists/Greg/DTNT';
else
    DTNTfilelistpath = 'N:\NexFiles\nexfilelists\Greg\DTNT';
end
DTNTdata = [];
filenames = fnamesFromTxt2([DTNTfilelistpath,filesep,'SednaDTNT25Hz.txt']);
for fileidx = 1:size(filenames,1)
    if (strncmp(filenames{fileidx},'sf',2))
        continue;
    end
    stro = nex2stro(findfile(char(filenames{fileidx})));
    [thresholds, colordirs, sfs, QuestTrajectories] = DTquestUnpackGH(stro, 'mode');
    threshs = thresholds';
    ncoldir = length(threshs);
    lmsmat = mkbasis(colordirs')';
    Loog = [];
    for j = 1:length(QuestTrajectories)
        Loog(j,1) = QuestTrajectories{j}(end) == max(QuestTrajectories{j});
    end
    DTNTdata = [DTNTdata; repmat(thresholds,1,3).*colordirs Loog];
end

DTNT_monspd = reshape(stro.sum.exptParams.mon_spect,81,3);
DTNT_M = reshape(stro.sum.exptParams.m_mtx,3,3);
DTNT_bkgndrgb = [stro.sum.exptParams.bkgnd_r  stro.sum.exptParams.bkgnd_g  stro.sum.exptParams.bkgnd_b];

% Everything forced into SMJ10
DTNTdata(:,[1 2 3]) = ConvertConeContrastBasis(DTNT_M, T_cones_smj10*DTNT_monspd, DTNT_bkgndrgb, DTNTdata(:,[1 2 3]));

Loog = logical(DTNTdata(:,end) == 1);
figure; axes; hold on;
plot3(DTNTdata(~Loog,1),DTNTdata(~Loog,2),DTNTdata(~Loog,3),'m.');
plot3(-DTNTdata(~Loog,1),-DTNTdata(~Loog,2),-DTNTdata(~Loog,3),'m.');
set(gca,'View',[-134 2]);

% Fitting DTNT data
[planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(DTNTdata(:,[1 2 3]), Loog);
DTNTweights = planeparams'*xformmat';
normDTNTweights = DTNTweights./sum(abs(DTNTweights));

% First gathering the neurothresh data
data = [];
[fnames, spikeIdx] = fnamesFromTxt2('N:\NexFiles\nexfilelists\Greg\NT\NTlum2.txt');
for cellcounter = 1:size(fnames,1)
    filename = findfile(char(fnames{cellcounter}));
    paradigmID = getparadigmID(filename);
    if paradigmID == 103
        NT = nex2stro(filename);
    end
    out = NTpreprocess(NT,0,Inf);
    fundamentals = NT.sum.exptParams.fundamentals;
    fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]);
    mon_spd = NT.sum.exptParams.mon_spd;
    mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);
    mon_spd = SplineRaw([380:4:780]', mon_spd, [380:5:780]');
    M = fundamentals'*mon_spd;
    data(cellcounter).data = out(:,[2 3 4]).*repmat(out(:,5),1,3);
    % Forcing everything into SMJ10
    data(cellcounter).data = ConvertConeContrastBasis(M, T_cones_smj10*DTNT_monspd, NT.sum.exptParams.bkgndrgb, data(cellcounter).data);
    
    % data(cellcounter).M = M;
    data(cellcounter).M = T_cones_smj10*DTNT_monspd;
    data(cellcounter).bkgndrgb = NT.sum.exptParams.bkgndrgb;
    data(cellcounter).monspd = mon_spd;
    data(cellcounter).Loog = logical(out(:,7));
end

figure; axes; hold on;
coneweights = [];
for i = 1:length(data)
    i
    tmp = data(i).data;
    Loog = data(i).Loog;
    plot3(tmp(~Loog,1),tmp(~Loog,2),tmp(~Loog,3),'m.');
    plot3(-tmp(~Loog,1),-tmp(~Loog,2),-tmp(~Loog,3),'m.');
    plot3(tmp(Loog,1),tmp(Loog,2),tmp(Loog,3),'y.');
    plot3(-tmp(Loog,1),-tmp(Loog,2),-tmp(Loog,3),'y.');
    [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(tmp(:,[1 2 3]), Loog);
    coneweights(i,:) = planeparams'*xformmat';
end
L = logical(coneweights(:,2) < 0);
coneweights(L,:) = -coneweights(L,:);
set(gca,'View',[-134 2],'Xlim',[-1 1],'Ylim',[-1 1],'Zlim',[-2 2]);
normconeweights = coneweights./repmat(sum(abs(coneweights),2),1,3)

figure; axes; hold on;
hist(normconeweights(:,3));
plot(normDTNTweights(3),0,'y*')

% % Which set of cone fundamentals was used in each NT experiment?
% for i = 1:length(data)
%    [T_cones_smj10*data(i).monspd  T_cones_smj*data(i).monspd data(i).M]
%    pause
% end

% 
% % --------------------------------------------
% % Switching to other sets of fundamentals
% % --------------------------------------------
load ('den_lens_ss');
load ('den_mac_ss');
lensat400nm = den_lens_ss(5);

macpigtransmittance = 1./(10.^(den_mac_ss./max(den_mac_ss)*0.095));
macpigtransmittance = ones(size(macpigtransmittance,1),1)
lenstransmittance = 1./(10.^(den_lens_ss));
fundnopig = T_cones_ss10'./repmat(macpigtransmittance,1,3)./repmat(lenstransmittance,1,3);
fundnopig = fundnopig./repmat(max(fundnopig),81,1);

lenstransmittance = 1./(10.^(den_lens_ss/den_lens_ss(5).*lensat400nm));
tmpfunds = fundnopig.*repmat(lenstransmittance,1,3);
tmpfunds = tmpfunds./repmat(max(tmpfunds),81,1);
Mnew = tmpfunds'*DTNT_monspd;
tmpscaled = ConvertConeContrastBasis(DTNT_M, Mnew, DTNT_bkgndrgb, DTNTdata(:,[1 2 3]));

Loog = DTNTdata(:,4);
[planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(tmpscaled, Loog);
DTNTweights = planeparams'*xformmat';
normDTNTweights = DTNTweights./sum(abs(DTNTweights))

figure; axes; hold on;
plot3(tmpscaled(~Loog,1),tmpscaled(~Loog,2),tmpscaled(~Loog,3),'m.');
plot3(-tmpscaled(~Loog,1),-tmpscaled(~Loog,2),-tmpscaled(~Loog,3),'m.');
set(gca,'View',[-134 2]);

% NT
coneweights = [];
for i = 1:length(data)
    i
    Mnew = tmpfunds'*data(i).monspd;
    tmpscaled = ConvertConeContrastBasis(data(i).M, Mnew, data(i).bkgndrgb, data(i).data);
    [planeparams, planeSSE, quadparams, quadSSE, xformmat] = NTsurfacefit(tmpscaled, data(i).Loog);
    coneweights(i,:) = planeparams'*xformmat';
end
L = logical(coneweights(:,2) < 0);
coneweights(L,:) = -coneweights(L,:);
normconeweights = coneweights./repmat(sum(abs(coneweights),2),1,3);
mean(normconeweights)

figure; axes; hold on;
hist(normconeweights(:,3));
plot(normDTNTweights(3),0,'y*')

