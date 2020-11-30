% As of now the script is pretty much like the daylight_props_5 but need to
% see what new new things can I add to this. I believe there is a relation
% between the interaction of STA and PC1 and how signals are combined
% between 2 different non-overlapping subunits. I am going to advance this
% script with mainanalysis.m 
% Author - Abhishek De, 11/17

% Creating different daylight spectra from the CIE daylight functions
close all; clearvars;
plot_counter = 1;
global illuminants T_xyz fundamentals
wave = 400:10:720; % Taking this bizzare range to match the wavelength ranges of the natural images 
dayBasis = ieReadSpectra('cieDaylightBasis',wave); % Daylight spectra basis functions from isetbio
num_spectras = 100;
x = linspace(0.25,0.40,num_spectras);
y = 2.870*x - 3.000*(x.*x) - 0.275;
coeff1 = (-1.3515-1.7703*x+5.9114*y)./(0.0241+0.2562*x-0.7341*y);
coeff2 = (0.0300-31.4424*x+30.0717*y)./(0.0241+0.2562*x-0.7341*y);
coeffs = cat(2,ones(num_spectras,1),coeff1',coeff2'); % Limiting the coefficients between 0 and 1
illuminants = coeffs * dayBasis';

load T_xyz1964.mat; % 1964 CIE color matching functions
T_xyz = T_xyz1964;
T_xyY = T_xyz./(repmat(sum(T_xyz),[3 1])+0.0001);
rgb = reshape([T_xyY(1,1:end-20)',T_xyY(2,1:end-20)',T_xyY(3,1:end-20)'],[size(T_xyY(:,1:end-20),2) 1 3]);
figure(plot_counter),subplot(231),patch(T_xyY(1,1:end-20),T_xyY(2,1:end-20),rgb); hold on; plot(0.33,0.34,'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0]);
plot(x,y,'k','Linewidth',2),xlabel('x'), ylabel('y'); title('CIE daylight');  hold off;
subplot(232),plot(wave,illuminants'), xlabel('Wavelength'), ylabel('Energy'); set(gca,'Xlim',[390 780]); title('Daylight spectra');
subplot(233),plot(coeff1,coeff2,'Linewidth',2), xlabel('M1'),ylabel('M2'); title('Coefficients');

% Loading the cone action spectra and moni tor spectral distributions
load fundamentals.mat;
load mon_spd;
fundamentals = reshape(fundamentals,[length(fundamentals)/3,3]); %1st column - L, 2nd- M, 3rd- S
mon_spd = reshape(mon_spd,[length(mon_spd)/3, 3]);  % MONITOR SPECTRAL DISTRIBUTION IN R,G,B
mon_spd = spline([380:4:780], mon_spd', [380:5:780]); % fitting a cubic spline
lo = find(380:5:780==400);
hi = find(380:5:780==720);
fundamentals = fundamentals(lo:2:hi,:); % Starting the fundamentals from 390 nm
mon_spd = mon_spd(:,lo:2:hi);
M = fundamentals'*mon_spd'; % matrix that converts RGB phosphor intensites to L,M,S cone excitations
figure(plot_counter),subplot(234), plot(wave,fundamentals','Linewidth',2), set(gca,'Xlim',[390 780]); xlabel('Wavelength'); title('fundamentals');
subplot(235), plot(wave,mon_spd','Linewidth',2), set(gca,'Xlim',[390 780]); xlabel('Wavelength'); title('Mon spd');
subplot(236), plot(380:5:780,T_xyz','Linewidth',2), set(gca,'Xlim',[380 780]); xlabel('Wavelength'); title('CIE color matching functions');
plot_counter = plot_counter + 1;


%% Beginning to use a natural hyperspectral image from Nascimento and Foster repository. Giving this a try, don't know as of now where I want to go with this
% but I haven't dealt with natural images before, so might be interesting
load('ref_cyflower1bb_reg1.mat'); % Loading the image file, size is 1019x1337x33
% load('ref_ribeira1bbb_reg1.mat');
% load('ref_ruivaes1bb_reg1.mat');
[R,C,D] = size(reflectances);
reflectances = RGB2XWFormat(reflectances);
idx = [1 25 50 75 100];
im = cell(1,numel(idx));
cexc = cell(1,numel(idx));
imXYZ = cell(1,numel(idx));
imxyY = cell(1,numel(idx));
bkgndxyY = [];
bkgndLMSexc = [];
spatialpts = randi(size(reflectances,1),[10000 2]); % Selecting 100 pairs of pts
Lconeexc = cell(1,numel(idx)); % Essential for storing the cone excitations
Mconeexc = cell(1,numel(idx));
Sconeexc = cell(1,numel(idx));
Lconecon = cell(1,numel(idx)); % Essential for storing the cone contrasts
Mconecon = cell(1,numel(idx));
Sconecon = cell(1,numel(idx));
Lconeratio = cell(1,numel(idx)); % Essential for storing cone ratios
Mconeratio = cell(1,numel(idx));
Sconeratio = cell(1,numel(idx));
LMconecon = cell(1,numel(idx)); % for storing L-M cone contrast of the selected pixels
Lumconecon = cell(1,numel(idx));
angulardiffinchromaticities = cell(1,numel(idx));
for ii = 1:numel(idx)
    [imtmp,cexctmp,imXYZtmp] = calcRGB(reflectances,idx(ii),lo,hi,R,C);
    im{ii} = imtmp;
    cexc{ii} = cexctmp;
    imXYZ{ii} = imXYZtmp;
    imxyYtmp = transpose(XYZToxyY(imXYZtmp'));
    imxyY{ii} = imxyYtmp;
    bkgndLMSexc = [bkgndLMSexc mean(cexctmp',2)];
    bkgndxyY = [bkgndxyY XYZToxyY(mean(imXYZtmp',2))];
    Lconeratio{ii} = cexctmp(spatialpts(:,1),1)./ cexctmp(spatialpts(:,2),1);
    Mconeratio{ii} = cexctmp(spatialpts(:,1),2)./ cexctmp(spatialpts(:,2),2);
    Sconeratio{ii} = cexctmp(spatialpts(:,1),3)./ cexctmp(spatialpts(:,2),3);
    Lconeexc{ii} = [cexctmp(spatialpts(:,1),1) cexctmp(spatialpts(:,2),1)]; % corresponds to photoreceptor adaptation
    Mconeexc{ii} = [cexctmp(spatialpts(:,1),2) cexctmp(spatialpts(:,2),2)];
    Sconeexc{ii} = [cexctmp(spatialpts(:,1),3) cexctmp(spatialpts(:,2),3)];
    Lconecon{ii} = (Lconeexc{ii}-repmat(bkgndLMSexc(1,end),size(spatialpts)))./repmat(bkgndLMSexc(1,end),size(spatialpts));
    Mconecon{ii} = (Mconeexc{ii}-repmat(bkgndLMSexc(2,end),size(spatialpts)))./repmat(bkgndLMSexc(2,end),size(spatialpts));
    Sconecon{ii} = (Sconeexc{ii}-repmat(bkgndLMSexc(3,end),size(spatialpts)))./repmat(bkgndLMSexc(3,end),size(spatialpts));
    LMconecon{ii} = Lconecon{ii} - Mconecon{ii}; % Useful for analysis later on, storing the L-M cone contrast of the selected pixels
    Lumconecon{ii} = Lconecon{ii} + Mconecon{ii};
    vec1 = imxyYtmp(spatialpts(:,1),1:2) - repmat(bkgndxyY(1:2,end)',[size(spatialpts,1) 1]);
    vec2 = imxyYtmp(spatialpts(:,2),1:2) - repmat(bkgndxyY(1:2,end)',[size(spatialpts,1) 1]);
    angulardiffinchromaticities{ii} = (180/pi)*acos(sum(vec1.*vec2,2)./(sqrt(sum(vec1.^2,2)).*sqrt(sum(vec2.^2,2))));
end
figure(plot_counter),subplot(231),image(im{1}.^0.23), set(gca,'XTick',[],'YTick',[]);
subplot(232),image(im{2}.^0.23), set(gca,'XTick',[],'YTick',[]);
subplot(233),image(im{3}.^0.23), set(gca,'XTick',[],'YTick',[]);
subplot(234),image(im{4}.^0.23), set(gca,'XTick',[],'YTick',[]);
subplot(235),image(im{5}.^0.23), set(gca,'XTick',[],'YTick',[]);
subplot(236),plot(wave,illuminants(idx,:)','Linewidth',2); xlabel('Wavelength');ylabel('Energy'); title('Illuminant spectras');
plot_counter = plot_counter + 1;

%%
% Plotting the cone ratios between 2 spatial regions in the image under 2 different illuminants
figure(plot_counter),subplot(331),plot(Lconeratio{1},Lconeratio{2},'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0]); hold on; lsline; xlabel('ill1');ylabel('ill2');title('L-cone ratio'); hold off;
subplot(332);plot(Mconeratio{1},Mconeratio{2},'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[0 1 0]); hold on; lsline; xlabel('ill1');ylabel('ill2');title('M-cone ratio'); hold off;
subplot(333);plot(Sconeratio{1},Sconeratio{2},'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1]); hold on; lsline; xlabel('ill1');ylabel('ill2');title('S-cone ratio'); hold off;
subplot(334);plot(Lconeratio{1},Lconeratio{3},'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0]); hold on; lsline; xlabel('ill1');ylabel('ill3');title('L-cone ratio'); hold off;
subplot(335);plot(Mconeratio{1},Mconeratio{3},'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[0 1 0]); hold on; lsline; xlabel('ill1');ylabel('ill3');title('M-cone ratio'); hold off;
subplot(336);plot(Sconeratio{1},Sconeratio{3},'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1]); hold on; lsline; xlabel('ill1');ylabel('ill3');title('S-cone ratio'); hold off;
subplot(337);plot(Lconeratio{1},Lconeratio{4},'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0]); hold on; lsline; xlabel('ill1');ylabel('ill4');title('L-cone ratio'); hold off;
subplot(338);plot(Mconeratio{1},Mconeratio{4},'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[0 1 0]); hold on; lsline; xlabel('ill1');ylabel('ill4');title('M-cone ratio'); hold off;
subplot(339);plot(Sconeratio{1},Sconeratio{4},'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1]); hold on; lsline; xlabel('ill1');ylabel('ill4');title('S-cone ratio'); hold off;
plot_counter = plot_counter + 1;

% plotting the differences in the chromaticities between 2 spatial regions under different illuminants
figure(plot_counter),subplot(241), plot(angulardiffinchromaticities{1},angulardiffinchromaticities{2},'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0]); hold on; lsline; xlabel('ill1'); ylabel('ill2');
subplot(242), plot(angulardiffinchromaticities{1},angulardiffinchromaticities{3},'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0]); hold on; lsline; xlabel('ill1'); ylabel('ill3');
subplot(243), plot(angulardiffinchromaticities{1},angulardiffinchromaticities{4},'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0]); hold on; lsline; xlabel('ill1'); ylabel('ill4');
subplot(244), plot(angulardiffinchromaticities{1},angulardiffinchromaticities{5},'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0]); hold on; lsline; xlabel('ill1'); ylabel('ill5');
subplot(245), hist(angulardiffinchromaticities{1}-angulardiffinchromaticities{2}); hold on; plot(mean(angulardiffinchromaticities{1}-angulardiffinchromaticities{2}),0,'kv','MarkerFacecolor','g'); title('ill1 - ill2'); xlabel('Angular diff'); set(gca,'Xlim',[-180 180]); hold off;
subplot(246), hist(angulardiffinchromaticities{1}-angulardiffinchromaticities{3}); hold on; plot(mean(angulardiffinchromaticities{1}-angulardiffinchromaticities{3}),0,'kv','MarkerFacecolor','g'); title('ill1 - ill3'); xlabel('Angular diff'); set(gca,'Xlim',[-180 180]); hold off;
subplot(247), hist(angulardiffinchromaticities{1}-angulardiffinchromaticities{3}); hold on; plot(mean(angulardiffinchromaticities{1}-angulardiffinchromaticities{4}),0,'kv','MarkerFacecolor','g'); title('ill1 - ill4'); xlabel('Angular diff'); set(gca,'Xlim',[-180 180]); hold off;
subplot(248), hist(angulardiffinchromaticities{1}-angulardiffinchromaticities{5}); hold on; plot(mean(angulardiffinchromaticities{1}-angulardiffinchromaticities{5}),0,'kv','MarkerFacecolor','g'); title('ill1 - ill5'); xlabel('Angular diff'); set(gca,'Xlim',[-180 180]); hold off;
plot_counter = plot_counter + 1;

% Plotting the chromaticities of the background under different illuminants
% Quantifying the change in the variance as the iluminant spectra shifts from having max energy at low wavelenght to high wavelength
stdofthehists = [std(angulardiffinchromaticities{1}-angulardiffinchromaticities{2}); std(angulardiffinchromaticities{1}-angulardiffinchromaticities{3});std(angulardiffinchromaticities{1}-angulardiffinchromaticities{4});std(angulardiffinchromaticities{1}-angulardiffinchromaticities{5})];
figure(plot_counter),subplot(121); patch(T_xyY(1,1:end-20),T_xyY(2,1:end-20),rgb); hold on; plot(bkgndxyY(1,:),bkgndxyY(2,:),'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0]);
plot(x,y,'k','Linewidth',2),xlabel('x'), ylabel('y'); title('Chromaticities of the background');  hold off;
subplot(122); plot(stdofthehists,'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0]); hold on; lsline; ylabel('st deviation'); set(gca,'XTick',[1 2 3 4],'XTickLabel',{'ill:1-2','ill:1-3','ill:1-4','ill:1-5'}); hold off;
plot_counter = plot_counter + 1;

% Now plotting to see if the differences in cone contrast between 2 spatial
% regions remains preserved as I change the illuminant. If the cone
% contrast differences are preserved, then I would expect a line, See the
% plot for more clarity
figure(plot_counter), subplot(331), plot(diff(Lconecon{1}'),diff(Lconecon{2}'),'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0]);  hold on; lsline; xlabel('ill1');ylabel('ill2');title('L-cone contrast diff'); hold off;
subplot(332), plot(diff(Mconecon{1}'),diff(Mconecon{2}'),'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[0 1 0]);  hold on; lsline; xlabel('ill1');ylabel('ill2');title('M-cone contrast diff'); hold off;
subplot(333), plot(diff(Sconecon{1}'),diff(Sconecon{2}'),'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1]);  hold on; lsline; xlabel('ill1');ylabel('ill2');title('S-cone contrast diff'); hold off;
subplot(334), plot(diff(Lconecon{1}'),diff(Lconecon{3}'),'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0]);  hold on; lsline; xlabel('ill1');ylabel('ill3');title('L-cone contrast diff'); hold off;
subplot(335), plot(diff(Mconecon{1}'),diff(Mconecon{3}'),'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[0 1 0]);  hold on; lsline; xlabel('ill1');ylabel('ill3');title('M-cone contrast diff'); hold off;
subplot(336), plot(diff(Sconecon{1}'),diff(Sconecon{3}'),'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1]);  hold on; lsline; xlabel('ill1');ylabel('ill3');title('S-cone contrast diff'); hold off;
subplot(337), plot(diff(Lconecon{1}'),diff(Lconecon{4}'),'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0]);  hold on; lsline; xlabel('ill1');ylabel('ill4');title('L-cone contrast diff'); hold off;
subplot(338), plot(diff(Mconecon{1}'),diff(Mconecon{4}'),'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[0 1 0]);  hold on; lsline; xlabel('ill1');ylabel('ill4');title('M-cone contrast diff'); hold off;
subplot(339), plot(diff(Sconecon{1}'),diff(Sconecon{4}'),'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1]);  hold on; lsline; xlabel('ill1');ylabel('ill4');title('S-cone contrast diff'); hold off;
plot_counter = plot_counter + 1;

% Quantifying Luminance and L-M cone contrasts present in the image under
% different illuminants using histograms
figure(plot_counter)
for ii = 1:5
    Lum = Lumconecon{ii};
    LM = LMconecon{ii};
    subplot(2,5,ii), hist(Lum(:),30), hold on; plot(mean(Lum(:)),0,'kv','MarkerFacecolor','g'); title('Lum var'); hold off;
    subplot(2,5,5+ii), hist(LM(:),30), hold on; plot(mean(LM(:)),0,'kv','MarkerFacecolor','r'); title('L-M var'); hold off;
end
plot_counter = plot_counter + 1;

% Wanted to see the L-M cone contrast differences between 2 spatial regions
% is preserved under different illuminants. Extending it to Luminance
% contrast also
figure(plot_counter), subplot(441), plot(diff(LMconecon{1}'),diff(LMconecon{2}'),'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0]);  hold on; lsline; xlabel('ill1');ylabel('ill2');title('L-M cone contrast diff'); hold off;
subplot(442), plot(diff(LMconecon{1}'),diff(LMconecon{3}'),'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[0 1 0]);  hold on; lsline; xlabel('ill1');ylabel('ill3');title('L-M cone contrast diff'); hold off;
subplot(443), plot(diff(LMconecon{1}'),diff(LMconecon{4}'),'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1]);  hold on; lsline; xlabel('ill1');ylabel('ill4');title('L-M cone contrast diff'); hold off;
subplot(444), plot(diff(LMconecon{1}'),diff(LMconecon{5}'),'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[0 1 1],'MarkerEdgeColor',[0 1 1]);  hold on; lsline; xlabel('ill1');ylabel('ill5');title('L-M cone contrast diff'); hold off;
subplot(445), hist(diff(LMconecon{1}')-diff(LMconecon{2}'),100); hold on; plot(mean(diff(LMconecon{1}')-diff(LMconecon{2}')),0,'kv','MarkerFacecolor','g'); title('ill1 - ill2'); xlabel('L-M cone contrast diff'); set(gca,'Xlim',[-0.2 0.2]); hold off;
subplot(446), hist(diff(LMconecon{1}')-diff(LMconecon{3}'),100); hold on; plot(mean(diff(LMconecon{1}')-diff(LMconecon{3}')),0,'kv','MarkerFacecolor','g'); title('ill1 - ill3'); xlabel('L-M cone contrast diff'); set(gca,'Xlim',[-0.2 0.2]); hold off;
subplot(447), hist(diff(LMconecon{1}')-diff(LMconecon{4}'),100); hold on; plot(mean(diff(LMconecon{1}')-diff(LMconecon{4}')),0,'kv','MarkerFacecolor','g'); title('ill1 - ill4'); xlabel('L-M cone contrast diff'); set(gca,'Xlim',[-0.2 0.2]); hold off;
subplot(448), hist(diff(LMconecon{1}')-diff(LMconecon{5}'),100); hold on; plot(mean(diff(LMconecon{1}')-diff(LMconecon{5}')),0,'kv','MarkerFacecolor','g'); title('ill1 - ill5'); xlabel('L-M cone contrast diff'); set(gca,'Xlim',[-0.2 0.2]); hold off;
subplot(449), plot(diff(Lumconecon{1}'),diff(Lumconecon{2}'),'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0]);  hold on; lsline; xlabel('ill1');ylabel('ill2');title('Lum cone contrast diff'); hold off;
subplot(4,4,10), plot(diff(Lumconecon{1}'),diff(Lumconecon{3}'),'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[0 1 0]);  hold on; lsline; xlabel('ill1');ylabel('ill3');title('Lum cone contrast diff'); hold off;
subplot(4,4,11), plot(diff(Lumconecon{1}'),diff(Lumconecon{4}'),'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1]);  hold on; lsline; xlabel('ill1');ylabel('ill4');title('Lum cone contrast diff'); hold off;
subplot(4,4,12), plot(diff(Lumconecon{1}'),diff(Lumconecon{5}'),'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[0 1 1],'MarkerEdgeColor',[0 1 1]);  hold on; lsline; xlabel('ill1');ylabel('ill5');title('Lum cone contrast diff'); hold off;
subplot(4,4,13), hist(diff(Lumconecon{1}')-diff(Lumconecon{2}'),100); hold on; plot(mean(diff(LMconecon{1}')-diff(LMconecon{2}')),0,'kv','MarkerFacecolor','g'); title('ill1 - ill2'); xlabel('L+M cone contrast diff'); set(gca,'Xlim',[-0.2 0.2]); hold off;
subplot(4,4,14), hist(diff(Lumconecon{1}')-diff(Lumconecon{3}'),100); hold on; plot(mean(diff(LMconecon{1}')-diff(LMconecon{3}')),0,'kv','MarkerFacecolor','g'); title('ill1 - ill3'); xlabel('L+M cone contrast diff'); set(gca,'Xlim',[-0.2 0.2]); hold off;
subplot(4,4,15), hist(diff(Lumconecon{1}')-diff(Lumconecon{4}'),100); hold on; plot(mean(diff(LMconecon{1}')-diff(LMconecon{4}')),0,'kv','MarkerFacecolor','g'); title('ill1 - ill4'); xlabel('L+M cone contrast diff'); set(gca,'Xlim',[-0.2 0.2]); hold off;
subplot(4,4,16), hist(diff(Lumconecon{1}')-diff(Lumconecon{5}'),100); hold on; plot(mean(diff(LMconecon{1}')-diff(LMconecon{5}')),0,'kv','MarkerFacecolor','g'); title('ill1 - ill5'); xlabel('L+M cone contrast diff'); set(gca,'Xlim',[-0.2 0.2]); hold off;
plot_counter = plot_counter + 1;

% Quantifying the increasing variance in the L-M differences as the
% illuminant progressively becomes more different; extending the work to L+M
stdofLMhists = [std(diff(LMconecon{1}')-diff(LMconecon{2}')); std(diff(LMconecon{1}')-diff(LMconecon{3}')); std(diff(LMconecon{1}')-diff(LMconecon{4}')); std(diff(LMconecon{1}')-diff(LMconecon{5}'))];
stdofLumhists = [std(diff(Lumconecon{1}')-diff(Lumconecon{2}')); std(diff(Lumconecon{1}')-diff(Lumconecon{3}')); std(diff(Lumconecon{1}')-diff(Lumconecon{4}')); std(diff(Lumconecon{1}')-diff(Lumconecon{5}'))];
figure(plot_counter),subplot(121), plot(stdofLMhists,'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0]); hold on; lsline; ylabel('st deviation'); set(gca,'XTick',[1 2 3 4],'XTickLabel',{'ill:1-2','ill:1-3','ill:1-4','ill:1-5'}); title('L-M'); hold off;
subplot(122), plot(stdofLumhists,'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0]); hold on; lsline; ylabel('st deviation'); set(gca,'XTick',[1 2 3 4],'XTickLabel',{'ill:1-2','ill:1-3','ill:1-4','ill:1-5'}); title('L+M'); hold off;
plot_counter = plot_counter + 1;

% Checking if the difference of cone excitations remains constant as I
% change the illuminant: comparing illuminant 1 and illuminant 5 
l1 = Lconeexc{1}; m1 = Mconeexc{1}; s1 = Sconeexc{1};
l5 = Lconeexc{5}; m5 = Mconeexc{5}; s5 = Sconeexc{5};
figure(plot_counter), subplot(231), plot(l1(:,1)-l1(:,2),l5(:,1)-l5(:,2),'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0]);  hold on; lsline; xlabel('ill1'); ylabel('ill5'); title('L exc diff'); hold off;
subplot(232), plot(m1(:,1)-m1(:,2),m5(:,1)-m5(:,2),'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[0 1 0]);  hold on; lsline; xlabel('ill1'); ylabel('ill5'); title('M exc diff'); hold off;
subplot(233), plot(s1(:,1)-s1(:,2),s5(:,1)-s5(:,2),'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1]);  hold on; lsline; xlabel('ill1'); ylabel('ill5'); title('S exc diff'); hold off;
subplot(234), plot(l1(:),l5(:),'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0]);  hold on; lsline; xlabel('ill1'); ylabel('ill5'); title('L exc '); hold off;
subplot(235), plot(m1(:),m5(:),'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[0 1 0]);  hold on; lsline; xlabel('ill1'); ylabel('ill5'); title('M exc '); hold off;
subplot(236), plot(s1(:),s5(:),'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1]);  hold on; lsline; xlabel('ill1'); ylabel('ill5'); title('S exc '); hold off;
plot_counter = plot_counter + 1;

% Next exercise will be fun. I want to take just two spatial pts from the image and look
% how their CIE chromaticities change as a function of changing illuminants
figure(plot_counter),subplot(221), patch(T_xyY(1,1:end-20),T_xyY(2,1:end-20),rgb); hold on;
for ii = 1:100
    surf1radiance = [];
    surf2radiance = [];
    surf1LMS = [];
    surf2LMS = [];
    surfidxs = randi(size(reflectances,1),[2 1]);
    for kk = 1:size(illuminants,1)
        surf1radiance =  [surf1radiance; illuminants(kk,:).*reflectances(surfidxs(1),:)];
        surf2radiance =  [surf2radiance; illuminants(kk,:).*reflectances(surfidxs(2),:)];
    end
    surf1LMS = surf1radiance*fundamentals;
    surf2LMS = surf2radiance*fundamentals;
    XYZsurf1 = surf1radiance*T_xyz(:,lo:2:hi)';
    XYZsurf2 = surf2radiance*T_xyz(:,lo:2:hi)';
    xyYsurf1 = XYZToxyY(XYZsurf1');
    xyYsurf2 = XYZToxyY(XYZsurf2');
    diffchromaticity = xyYsurf1 - xyYsurf2;
    subplot(221),plot(xyYsurf1(1,:),xyYsurf1(2,:),'o','MarkerSize',1,'LineWidth',0.5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0]); hold on;
    subplot(221),plot(xyYsurf2(1,:),xyYsurf2(2,:),'o','MarkerSize',1,'LineWidth',0.5,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[0 1 0]); hold on;
    subplot(222), plot(diffchromaticity(1,:),diffchromaticity(2,:),'k','Linewidth',2); hold on;
    subplot(223),plot3(surf1LMS(:,1),surf1LMS(:,2),surf1LMS(:,3),'o','MarkerSize',1,'LineWidth',0.5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0]); hold on;
    subplot(223),plot3(surf2LMS(:,1),surf2LMS(:,2),surf2LMS(:,3),'o','MarkerSize',1,'LineWidth',0.5,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[0 1 0]); hold on;
    subplot(224),plot3(surf1LMS(:,1)-surf2LMS(:,1),surf1LMS(:,2)-surf2LMS(:,2),surf1LMS(:,3)-surf2LMS(:,3),'k', 'LineWidth',1); hold on;
end
subplot(221), xlabel('x'), ylabel('y'); title('CIE daylight');  hold off;
subplot(222),xlabel('\delta x'), ylabel('\delta y'); title('Diff chromaticity'); hold off;
subplot(223),xlabel('L exc'), ylabel('M exc'); zlabel('S exc'); hold off;
subplot(224),xlabel('\delta L exc'), ylabel('\delta M exc'); zlabel('\delta S exc'); hold off;
plot_counter = plot_counter + 1;

%% In this section I will be dealing with munsell chips
% Load up the munsell chips, will do a head to head comparison with the
% hyperspectral images 
load munsell380_800_1.mat % rows - wavelength, columns - different reflectances
lo1= find(380:1:780==400);
hi1 = find(380:1:780==720);
munsell = munsell(lo1:10:hi1,:)'; % all the reflectances are column vectors
totmunsells = size(munsell,1);
munsellsurfs = randi(totmunsells,[1000 2]);
idx1 = [1 100];
munsellLconeexc = cell(1,numel(idx1));
munsellMconeexc = cell(1,numel(idx1));
munsellSconeexc = cell(1,numel(idx1));
munsellLconeratio = cell(1,numel(idx1));
munsellMconeratio = cell(1,numel(idx1));
munsellSconeratio = cell(1,numel(idx1));

for ii = 1:numel(idx1)
    radiance = [];
    for jj = 1:totmunsells
        radiance = [radiance; illuminants(idx1(ii),:).*munsell(jj,:)];
    end
    LMS = radiance*fundamentals;
    munsellLconeratio{ii} = LMS(munsellsurfs(:,1),1)./LMS(munsellsurfs(:,2),1);
    munsellMconeratio{ii} = LMS(munsellsurfs(:,1),2)./LMS(munsellsurfs(:,2),2);
    munsellSconeratio{ii} = LMS(munsellsurfs(:,1),3)./LMS(munsellsurfs(:,2),3);
    munsellLconeexc{ii} = [LMS(munsellsurfs(:,1),1) LMS(munsellsurfs(:,2),1)];
    munsellMconeexc{ii} = [LMS(munsellsurfs(:,1),2) LMS(munsellsurfs(:,2),2)];
    munsellSconeexc{ii} = [LMS(munsellsurfs(:,1),3) LMS(munsellsurfs(:,2),3)];
end

% First comaring the cone ratios and then the actual cone excitations
l1 = munsellLconeexc{1}; m1 = munsellMconeexc{1}; s1 = munsellSconeexc{1};
l2 = munsellLconeexc{2}; m2 = munsellMconeexc{2}; s2 = munsellSconeexc{2};
figure(plot_counter); subplot(331),plot(munsellLconeratio{1},munsellLconeratio{2},'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0]);  hold on; lsline; xlabel('ill1'); ylabel('ill5'); title('L cone ratio'); hold off;
subplot(332),plot(munsellMconeratio{1},munsellMconeratio{2},'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[0 1 0]);  hold on; lsline; xlabel('ill1'); ylabel('ill5'); title('M cone ratio'); hold off;
subplot(333),plot(munsellSconeratio{1},munsellSconeratio{2},'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1]);  hold on; lsline; xlabel('ill1'); ylabel('ill5'); title('S cone ratio'); hold off;
subplot(334),plot(l1(:),l2(:),'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0]);  hold on; lsline; xlabel('ill1'); ylabel('ill5'); title('L cone'); hold off;
subplot(335),plot(m1(:),m2(:),'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[0 1 0]);  hold on; lsline; xlabel('ill1'); ylabel('ill5'); title('M cone exc'); hold off;
subplot(336),plot(s1(:),s2(:),'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1]);  hold on; lsline; xlabel('ill1'); ylabel('ill5'); title('S cone exc'); hold off;
subplot(337),plot(l1(:,1)-l1(:,2),l2(:,1)-l2(:,2),'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0]);  hold on; lsline; xlabel('ill1'); ylabel('ill5'); title('L exc diff'); hold off;
subplot(338),plot(m1(:,1)-m1(:,2),m2(:,1)-m2(:,2),'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[0 1 0]);  hold on; lsline; xlabel('ill1'); ylabel('ill5'); title('M exc diff'); hold off;
subplot(339),plot(s1(:,1)-s1(:,2),s2(:,1)-s2(:,2),'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1]);  hold on; lsline; xlabel('ill1'); ylabel('ill5'); title('S exc diff'); hold off;
plot_counter = plot_counter + 1;

% Next exercise will be fun. I want to take just two munsell edges and look
% how their CIE chromaticities change as a function of changing illuminants
figure(plot_counter),subplot(221), patch(T_xyY(1,1:end-20),T_xyY(2,1:end-20),rgb); hold on;
for ii = 1:50
    surf1radiance = [];
    surf2radiance = [];
    surf1LMS = [];
    surf2LMS = [];
    surfidxs = randi(totmunsells,[2 1]);
    for kk = 1:size(illuminants,1)
        surf1radiance =  [surf1radiance; illuminants(kk,:).*munsell(surfidxs(1),:)];
        surf2radiance =  [surf2radiance; illuminants(kk,:).*munsell(surfidxs(2),:)];
    end
    surf1LMS = surf1radiance*fundamentals;
    surf2LMS = surf2radiance*fundamentals;
    XYZsurf1 = surf1radiance*T_xyz(:,lo:2:hi)';
    XYZsurf2 = surf2radiance*T_xyz(:,lo:2:hi)';
    xyYsurf1 = XYZToxyY(XYZsurf1');
    xyYsurf2 = XYZToxyY(XYZsurf2');
    diffchromaticity = xyYsurf1 - xyYsurf2;
    subplot(221),plot(xyYsurf1(1,:),xyYsurf1(2,:),'o','MarkerSize',1,'LineWidth',0.5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0]); hold on;
    subplot(221),plot(xyYsurf2(1,:),xyYsurf2(2,:),'o','MarkerSize',1,'LineWidth',0.5,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[0 1 0]); hold on;
    subplot(222), plot(diffchromaticity(1,:),diffchromaticity(2,:),'k','Linewidth',1); hold on;
    subplot(223),plot3(surf1LMS(:,1),surf1LMS(:,2),surf1LMS(:,3),'o','MarkerSize',1,'LineWidth',0.5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0]); hold on;
    subplot(223),plot3(surf2LMS(:,1),surf2LMS(:,2),surf2LMS(:,3),'o','MarkerSize',1,'LineWidth',0.5,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[0 1 0]); hold on;
    subplot(224),plot3(surf1LMS(:,1)-surf2LMS(:,1),surf1LMS(:,2)-surf2LMS(:,2),surf1LMS(:,3)-surf2LMS(:,3),'k', 'LineWidth',1); hold on;
    
end
subplot(221),xlabel('x'), ylabel('y'); title('CIE daylight');  hold off;
subplot(222),xlabel('\delta x'), ylabel('\delta y'); title('Diff chromaticity'); hold off;
subplot(223),xlabel('L exc'), ylabel('M exc'); zlabel('S exc'); hold off;
subplot(224),xlabel('\delta L exc'), ylabel('\delta M exc'); zlabel('\delta S exc'); hold off;
plot_counter = plot_counter + 1;

%% From the above set of simulations it is clear that the difference in cone excitations across changes in illuminants has a very high correlation 
% suggesting difference in cone excitations is preserved. Next step is to understand if the double opponent cells are calculating a difference 
% in cone excitations or difference in chromaticities 

bkgndRGB = [0.5 0.5 0.5];
RGBS1 = [0.5 -0.5 0]; % Subunit 1 of DO cell
RGBS2 = [-0.5 0.5 0]; % Subunit 2 of DO cell
x = linspace(-1,1,51); % x represents the contrast (RGB amplitudes) of subunit 1  
y = -x + 0.5; % y represents the contrast (RGB amplitudes) of subunit 2
S1RGBvals = RGBS1'*x;
S1RGBvals = S1RGBvals + repmat(bkgndRGB',[1 size(S1RGBvals,2)]);
S2RGBvals = RGBS2'*y;
S2RGBvals = S2RGBvals + repmat(bkgndRGB',[1 size(S2RGBvals,2)]);
S1xyY = XYZToxyY(SRGBPrimaryToXYZ(S1RGBvals));
S2xyY = XYZToxyY(SRGBPrimaryToXYZ(S2RGBvals));
S1coneexc = M * S1RGBvals;
S2coneexc = M * S2RGBvals;
diffchromaticities = S1xyY - S2xyY;
diffconeexcitations = S1coneexc - S2coneexc;
figure(plot_counter),subplot(221),plot(x,y,'Linewidth',2); hold on; 
plot([-1 1],[0 0],'-k','Linewidth',2); hold on; % x- axis 
plot([0 0],[-1 1],'-k','Linewidth',2); % y- axis
xlabel('S1 RGB contrast'); ylabel('S2 RGB contrast'); title('Isoresponse contour'); hold off;
subplot(222), patch(T_xyY(1,1:end-20),T_xyY(2,1:end-20),rgb); hold on;
plot(S1xyY(1,:),S1xyY(2,:),'o','MarkerSize',1,'LineWidth',0.5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0]); hold on;
plot(S2xyY(1,:),S2xyY(2,:),'o','MarkerSize',1,'LineWidth',0.5,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[0 1 0]); hold on;
xlabel('x'), ylabel('y'); title('CIE daylight');  hold off;
subplot(223),plot(diffchromaticities(1,:),diffchromaticities(2,:),'o','MarkerSize',2,'LineWidth',0.5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0]); hold on;
xlabel('\delta x'), ylabel('\delta y'); title('Diff chromaticity model'); hold off;
subplot(224),plot3(diffconeexcitations(1,:),diffchromaticities(2,:), diffchromaticities(3,:),'o','MarkerSize',2,'LineWidth',0.5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0]); hold on;
xlabel('\delta L'), ylabel('\delta M'), zlabel('\delta S'); title('Diff cone excitations model'); hold off;