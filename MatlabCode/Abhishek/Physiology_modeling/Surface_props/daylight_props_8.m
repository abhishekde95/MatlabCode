% This is new script which has been long time due. The goal of this script
% to see which type of double opponent cell i.e S+M vs S+L give rise to
% linear/nonlinear isoresponse contour. Might be relevant for CoSYNE
% Author - Abhishek De, 2/18,
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

%% Here's when I create a lime-magenta and orange-cyan cell
iters = 1;
trials = 100;%size(munsellsurfs,1);
LMDOrespS1 = cell(trials,1); LMDOrespS2 = cell(trials,1);
OCDOrespS1 = cell(trials,1); OCDOrespS2 = cell(trials,1);
AreaOC = zeros(trials,1); AreaLM = zeros(trials,1);
numsubplots = ceil(sqrt(iters));
mean_cone_exc = mean(munsell(:))*(fundamentals'*dayBasis(:,1));
ratioeigLM = zeros(trials,1);
ratioeigOC = zeros(trials,1);
binsizes1 = linspace(95,100,6);
binsizes2 = linspace(95,100,6);
idx = [0 0];
Areadiff = zeros(iters,1); % takes the difference between the LM and OC cdfs
for kk = 1:iters
    tmp = rand(1,2);
    tmp = tmp/norm(tmp);
    %     LMDO_S1 = [0.5 -0.5 0.25];
    %     OCDO_S1 = [0.5 -0.5 -0.25];
    LMDO_S1 = [tmp(1) -tmp(1) tmp(2)];
    OCDO_S1 = [-tmp(1) tmp(1) tmp(2)];
    LMDO_S2 = -LMDO_S1;
    OCDO_S2 = -OCDO_S1;
    for jj = 1:trials
        % Here is where I will select the munsell edges
        LMSexc1 = []; LMSexc2 = [];
        LMScon1 = []; LMScon2 = [];
        %         while idx(1) == idx(2)
        %             idx = randi(totmunsells,[1 2]);
        %         end
        for ii = 1:size(illuminants,1)
            surf1 = munsell(munsellsurfs(jj,1),:); % reflectance spectra of surface 1
            surf2 = munsell(munsellsurfs(jj,2),:); % reflectance spectra of surface 2
            lightfromsurf1 = illuminants(ii,:).*surf1;
            lightfromsurf2 = illuminants(ii,:).*surf2;
            % LMS calculations
            LMSexc1 = [LMSexc1 fundamentals'*lightfromsurf1'];
            LMSexc2 = [LMSexc2 fundamentals'*lightfromsurf2'];
            %         mean_cone_exc = mean([LMSexc1(:,end) LMSexc2(:,end)],2);
            LMScon1 = [LMScon1 (LMSexc1(:,end) - mean_cone_exc)./mean_cone_exc];
            LMScon2 = [LMScon2 (LMSexc2(:,end) - mean_cone_exc)./mean_cone_exc];
        end
        LMDOrespS1{jj} = LMDO_S1 * LMScon1;
        LMDOrespS2{jj} = LMDO_S2 * LMScon2;
        OCDOrespS1{jj} = OCDO_S1 * LMScon1;
        OCDOrespS2{jj} = OCDO_S2 * LMScon2;
        
        % Lime-Magenta
        tmpvar = [LMDOrespS1{jj}; LMDOrespS2{jj}];
        [~,eigval] = eig(cov(tmpvar'));
        ratioeigLM(jj) = 100*max(diag(eigval))/sum(diag(eigval));
        
        % Orange-Cyan
        tmpvar = [OCDOrespS1{jj}; OCDOrespS2{jj}];
        [~,eigval] = eig(cov(tmpvar'));
        ratioeigOC(jj) = 100*max(diag(eigval))/sum(diag(eigval));
    end
    [aLM,bLM] = ecdf(log10(ratioeigLM));
    [aOC,bOC] = ecdf(log10(ratioeigOC));
    figure(plot_counter),subplot(numsubplots,numsubplots,kk);histogram(ratioeigLM,binsizes2); hold on; histogram(ratioeigOC,binsizes2);hold off;
    figure(plot_counter+1),subplot(numsubplots,numsubplots,kk);plot(bLM,aLM); hold on; plot(bOC,aOC,'r'); hold off;
    Areadiff(kk) = trapz(bLM(~isinf(bLM)),aLM(~isinf(bLM))) - trapz(bOC(~isinf(bOC)),aOC(~isinf(bOC)));
    tmp = max(bOC(~isinf(bOC))) - max(bLM(~isinf(bLM)));
    Areadiff(kk) = Areadiff(kk) + tmp;
    
end
plot_counter = plot_counter + 2;
Areanull = Areadiff - mean(Areadiff)*ones(size(Areadiff));
figure(plot_counter); ecdf(Areanull); hold on; ecdf(Areadiff); set(gcf,'Name','LM and OC cdf area difference');hold off;
plot_counter = plot_counter + 1;

figure(plot_counter);
for ii= 1:trials
    subplot(121);plot(LMDOrespS1{ii},LMDOrespS2{ii},'r'); hold on;
    subplot(122);plot(OCDOrespS2{ii},OCDOrespS1{ii},'b'); hold on;
end
subplot(121); axis equal; axis square; grid on; hold off;
subplot(122); axis equal; axis square; grid on; hold off;
plot_counter = plot_counter + 1;

ratioeigOCtmp = hist(ratioeigOC,binsizes2);
ratioeigLMtmp = hist(ratioeigLM,binsizes1);
figure(plot_counter);bar(binsizes1,[ratioeigOCtmp' ratioeigLMtmp']); legend('OC','LM');
plot_counter = plot_counter + 1;

%% Pretty much the same simulation but here I am using the linear and quadratic fits to the trajectories to decide how linear or non-linear the
% trajectories are. I am doing this for lime-magenta and orange-cyan
% neurons

iters = 1;
trials = 1;%size(munsellsurfs,1);
meandiff = zeros(iters,1);
for aa = 1:iters
    munsellsurfs = randi(totmunsells,[1000 2]);
    mean_cone_exc = mean(munsell(:))*(fundamentals'*dayBasis(:,1));
    tmp = rand(1,2);
    tmp = tmp/norm(tmp);
%         LMDO_S1 = [0.5 -0.5 0.5];
%         OCDO_S1 = [0.5 -0.5 -0.5];
    LMDO_S1 = [tmp(1) -tmp(1) tmp(2)];
    OCDO_S1 = [-tmp(1) tmp(1) tmp(2)];
    LMDO_S2 = -LMDO_S1;
    OCDO_S2 = -OCDO_S1;
    SSE_linearmodelOC = [];
    SSE_linearmodelLM = [];
    SSE_quadmodelOC = [];
    SSE_quadmodelLM = [];
    
    trials = 20;
    LMDOrespS1 = cell(trials,1); LMDOrespS2 = cell(trials,1);
    OCDOrespS1 = cell(trials,1); OCDOrespS2 = cell(trials,1);
    hcsurfaces = zeros(trials,2);
    for jj = 1:trials
        disp(jj);
        % Here is where I will select the munsell edges
        LMSexc1 = []; LMSexc2 = [];
        LMScon1 = []; LMScon2 = [];
        for ii = 1:size(illuminants,1)
            surf1 = munsell(munsellsurfs(jj,1),:); % reflectance spectra of surface 1
            surf2 = munsell(munsellsurfs(jj,2),:); % reflectance spectra of surface 2
            lightfromsurf1 = illuminants(ii,:).*surf1;
            lightfromsurf2 = illuminants(ii,:).*surf2;
            % LMS calculations
            LMSexc1 = [LMSexc1 fundamentals'*lightfromsurf1'];
            LMSexc2 = [LMSexc2 fundamentals'*lightfromsurf2'];
            LMScon1 = [LMScon1 (LMSexc1(:,end) - mean_cone_exc)./mean_cone_exc];
            LMScon2 = [LMScon2 (LMSexc2(:,end) - mean_cone_exc)./mean_cone_exc];
        end
%         keyboard;
        LMDOrespS1{jj} = LMDO_S1 * LMScon1;
        LMDOrespS2{jj} = LMDO_S2 * LMScon2;
        OCDOrespS1{jj} = OCDO_S1 * LMScon1;
        OCDOrespS2{jj} = OCDO_S2 * LMScon2;
        LMDOrespS1{jj}  = (LMDOrespS1{jj} - mean(LMDOrespS1{jj}))+2;
        LMDOrespS2{jj}  = (LMDOrespS2{jj} - mean(LMDOrespS2{jj}))+2;
        OCDOrespS1{jj}  = (OCDOrespS1{jj} - mean(OCDOrespS1{jj}))+2;
        OCDOrespS2{jj}  = (OCDOrespS2{jj} - mean(OCDOrespS2{jj}))+2;
        
         % linear:orange-cyan
        [THETAOC,RHOOC] = cart2pol(OCDOrespS1{jj},OCDOrespS2{jj});
        initguessOC = [OCDOrespS1{jj}' OCDOrespS2{jj}']\ones(numel(OCDOrespS1{jj}),1);
        [linmodelOC,~] = linefit_AD2(RHOOC', THETAOC', logical(zeros(size(OCDOrespS1{jj}))'),logical(ones(size(OCDOrespS1{jj}))'),initguessOC');
        RHOlinOC = 1./(linmodelOC*[cos(THETAOC); sin(THETAOC)]);
        LOOGtmp = RHOlinOC<0;
        [x_linOC,y_linOC] = pol2cart(THETAOC(~LOOGtmp)',RHOlinOC(~LOOGtmp)');
        RHOlinOC(LOOGtmp) = 10000000;
        SSE_linearmodelOC = [SSE_linearmodelOC; sum((log(RHOlinOC)-log(RHOOC)).^2)];
        
        initguess3 = [0 0 0 linmodelOC];
        [quadmodelOC,~] = quadfit_AD2(RHOOC', THETAOC', logical(zeros(size(OCDOrespS1{jj}))'), logical(ones(size(OCDOrespS1{jj}))'),initguess3);
        A = quadmodelOC(1); B = quadmodelOC(2); C = quadmodelOC(3); D = quadmodelOC(4); E = quadmodelOC(5);
        RHOquadOC = [];
        p = [A*cos(THETAOC').^2+B*sin(THETAOC').^2+C*(cos(THETAOC').*sin(THETAOC')) D*cos(THETAOC')+E*sin(THETAOC') -1*ones(numel(THETAOC'),1)];
        for kk = 1:size(p,1)
            rts = roots(p(kk,:));
            if all(rts>0)
                r = min(rts); % if both the roots are positive
            else
                r = max(rts);
            end
            RHOquadOC = [RHOquadOC; r];
        end
        L = RHOquadOC>0 & RHOquadOC==real(RHOquadOC);
        RHOquadOC(~L) = 10000000;
        [x_quadOC,y_quadOC] = pol2cart(THETAOC(L),RHOquadOC(L)');
        SSE_quadmodelOC = [SSE_quadmodelOC; sum((log(RHOquadOC')-log(RHOOC)).^2)];
        
        % linear:lime-magenta
        [THETALM,RHOLM] = cart2pol(LMDOrespS1{jj},LMDOrespS2{jj});
        initguessLM = [LMDOrespS1{jj}' LMDOrespS2{jj}']\ones(numel(LMDOrespS1{jj}),1);
        [linmodelLM,~] = linefit_AD2(RHOLM', THETALM', logical(zeros(size(LMDOrespS1{jj}))'),logical(ones(size(LMDOrespS1{jj}))'),initguessLM');
        RHOlinLM = 1./(linmodelLM*[cos(THETALM); sin(THETALM)]);
        LOOGtmp= RHOlinLM<0;
        [x_linLM,y_linLM] = pol2cart(THETALM(~LOOGtmp)',RHOlinLM(~LOOGtmp)');
        RHOlinLM(LOOGtmp) = 10000000;
        SSE_linearmodelLM = [SSE_linearmodelLM; sum((log(RHOlinLM)-log(RHOLM)).^2)];
        
        initguess3 = [0 0 0 linmodelLM];
        [quadmodelLM,~] = quadfit_AD2(RHOLM', THETALM', logical(zeros(size(LMDOrespS1{jj}))'), logical(ones(size(LMDOrespS1{jj}))'),initguess3);
        A = quadmodelLM(1); B = quadmodelLM(2); C = quadmodelLM(3); D = quadmodelLM(4); E = quadmodelLM(5);
        RHOquadLM = [];
        p = [A*cos(THETALM').^2+B*sin(THETALM').^2+C*(cos(THETALM').*sin(THETALM')) D*cos(THETALM')+E*sin(THETALM') -1*ones(numel(THETALM'),1)];
        for kk = 1:size(p,1)
            rts = roots(p(kk,:));
            if all(rts>0)
                r = min(rts); % if both the roots are positive
            else
                r = max(rts);
            end
            RHOquadLM = [RHOquadLM; r];
        end
        L = RHOquadLM>0 & RHOquadLM==real(RHOquadLM);
        if ~all(L==1)
            disp(jj);
        end
        RHOquadLM(~L) = 10000000;
        [x_quadLM,y_quadLM] = pol2cart(THETALM(L),RHOquadLM(L)');
        SSE_quadmodelLM = [SSE_quadmodelLM; sum((log(RHOquadLM')-log(RHOLM)).^2)];
        
          
        hcsurfaces(jj,1) = ~all(sort(THETALM) == THETALM);
        hcsurfaces(jj,2) = ~all(sort(THETAOC) == THETAOC);
        
        figure(plot_counter); subplot(221);plot(THETALM,RHOLM,'o'); hold on; plot(THETALM,RHOlinLM,'g','Linewidth',2); plot(THETALM,RHOquadLM,'b','Linewidth',2); hold off;
        subplot(222);plot(LMDOrespS1{jj},LMDOrespS2{jj},'o'); hold on; plot(x_linLM,y_linLM,'g','Linewidth',2); plot(x_quadLM,y_quadLM,'b','Linewidth',2); hold off;
        subplot(223);plot(THETAOC,RHOOC,'o'); hold on; plot(THETAOC,RHOlinOC,'g','Linewidth',2); plot(THETAOC,RHOquadOC,'b','Linewidth',2); hold off; drawnow;
        subplot(224);plot(OCDOrespS1{jj},OCDOrespS2{jj},'o'); hold on; plot(x_linOC,y_linOC,'g','Linewidth',2); plot(x_quadOC,y_quadOC,'b','Linewidth',2); hold off;
        plot_counter = plot_counter + 1;
        
    end
    
%     keyboard;
    rratioLM = SSE_linearmodelLM./SSE_quadmodelLM;
    rratioLM = rratioLM(imag(rratioLM)==0);
%     rratioLM(logical(hcsurfaces(:,1)))= 1000;
    rratioLMtmp = log10(rratioLM);
    rratioLMtmp(logical(hcsurfaces(:,1)))= 1000;
    
    rratioOC = SSE_linearmodelOC./SSE_quadmodelOC;
    rratioOC = rratioOC(imag(rratioOC)==0);
    rratioOCtmp = log10(rratioOC);
    rratioOCtmp(logical(hcsurfaces(:,2)))= 1000;
    
    meandiff(aa) = mean(rratioLMtmp) - mean(rratioOCtmp);
    binsizes = 0:0.2:3;
    [cLM,~] = hist(rratioLMtmp,binsizes);
    [cOC,~] = hist(rratioOCtmp,binsizes);
    figure(plot_counter);bar(binsizes,[cOC' cLM']); legend('OC','LM');
    plot_counter = plot_counter + 1;
%     figure(plot_counter),histogram(rratioOCtmp,binsizes,'Normalization','probability'); hold on; histogram(rratioLMtmp,binsizes,'Normalization','probability');
%     plot(mean(rratioLMtmp),0,'gv'); plot(mean(rratioOCtmp),0,'bv');hold off;
    plot_counter = plot_counter + 1;
    
end



figure(plot_counter);
for ii = 1:trials
    subplot(121);plot(LMDOrespS1{ii},LMDOrespS2{ii},'r'); hold on;
    subplot(122);plot(OCDOrespS2{ii},OCDOrespS1{ii},'b'); hold on;
end
subplot(121); axis equal; axis square; grid on; title('Lime-Magenta'); hold off;
subplot(122); axis equal; axis square; grid on; title('Orange-Cyan');hold off;
plot_counter = plot_counter + 1;

%% Pretty much the same simulation but here I am using the linear and quadratic fits in cartesian plane

iters = 1;
munsellsurfs = randi(totmunsells,[10000 2]);
trials = size(munsellsurfs,1);
meandiff = zeros(iters,1);
for aa = 1:iters
    mean_cone_exc = mean(munsell(:))*(fundamentals'*dayBasis(:,1));
    tmp = rand(1,2);
    tmp = tmp/norm(tmp);
%         LMDO_S1 = [0.5 -0.5 0.5];
%         OCDO_S1 = [0.5 -0.5 -0.5];
    LMDO_S1 = [tmp(1) -tmp(1) tmp(2)];
    OCDO_S1 = [-tmp(1) tmp(1) tmp(2)];
    LMDO_S2 = -LMDO_S1;
    OCDO_S2 = -OCDO_S1;
    SSE_linearmodelOC = [];
    SSE_linearmodelLM = [];
    SSE_quadmodelOC = [];
    SSE_quadmodelLM = [];
    
    LMDOrespS1 = cell(trials,1); LMDOrespS2 = cell(trials,1);
    OCDOrespS1 = cell(trials,1); OCDOrespS2 = cell(trials,1);
    hcsurfaces = zeros(trials,2);
    for jj = 1:trials
        disp(jj);
        % Here is where I will select the munsell edges
        LMSexc1 = []; LMSexc2 = [];
        LMScon1 = []; LMScon2 = [];
        for ii = 1:size(illuminants,1)
            surf1 = munsell(munsellsurfs(jj,1),:); % reflectance spectra of surface 1
            surf2 = munsell(munsellsurfs(jj,2),:); % reflectance spectra of surface 2
            lightfromsurf1 = illuminants(ii,:).*surf1;
            lightfromsurf2 = illuminants(ii,:).*surf2;
            % LMS calculations
            LMSexc1 = [LMSexc1 fundamentals'*lightfromsurf1'];
            LMSexc2 = [LMSexc2 fundamentals'*lightfromsurf2'];
            LMScon1 = [LMScon1 (LMSexc1(:,end) - mean_cone_exc)./mean_cone_exc];
            LMScon2 = [LMScon2 (LMSexc2(:,end) - mean_cone_exc)./mean_cone_exc];
        end
%         keyboard;
        LMDOrespS1{jj} = LMDO_S1 * LMScon1;
        LMDOrespS2{jj} = LMDO_S2 * LMScon2;
        OCDOrespS1{jj} = OCDO_S1 * LMScon1;
        OCDOrespS2{jj} = OCDO_S2 * LMScon2;
        LMDOrespS1{jj}  = (LMDOrespS1{jj} - mean(LMDOrespS1{jj}))+2;
        LMDOrespS2{jj}  = (LMDOrespS2{jj} - mean(LMDOrespS2{jj}))+2;
        OCDOrespS1{jj}  = (OCDOrespS1{jj} - mean(OCDOrespS1{jj}))+2;
        OCDOrespS2{jj}  = (OCDOrespS2{jj} - mean(OCDOrespS2{jj}))+2;
        
        % linear:lime-magenta
        [linLM,statslinLM] = fit(LMDOrespS1{jj}',LMDOrespS2{jj}','poly1');
        linmodelLM = [linLM.p1 linLM.p2];
        predlinLM = linmodelLM(1)*LMDOrespS1{jj} + linmodelLM(2);
        SSE_linearmodelLM = [SSE_linearmodelLM; sum((predlinLM-LMDOrespS2{jj}).^2)];
        % non-linear:lime-magenta
        [quadLM,statsquadLM] = fit(LMDOrespS1{jj}',LMDOrespS2{jj}','poly2');
        quadmodelLM = [quadLM.p1 quadLM.p2 quadLM.p3];
        predquadLM = quadmodelLM(1)*LMDOrespS1{jj}.^2 + quadmodelLM(2)*LMDOrespS1{jj} + quadmodelLM(3)*ones(size(LMDOrespS1{jj}));
        SSE_quadmodelLM = [SSE_quadmodelLM; sum((predquadLM-LMDOrespS2{jj}).^2)];
        
        % linear:orange-cyan
        [linOC,statslinOC] = fit(OCDOrespS1{jj}',OCDOrespS2{jj}','poly1');
        linmodelOC = [linOC.p1 linOC.p2];
        predlinOC = linmodelOC(1)*OCDOrespS1{jj} + linmodelOC(2);
        SSE_linearmodelOC = [SSE_linearmodelOC; sum((predlinOC-OCDOrespS2{jj}).^2)];
        % non-linear:lime-magenta
        [quadOC,statsquadOC] = fit(OCDOrespS1{jj}',OCDOrespS2{jj}','poly2');
        quadmodelOC = [quadOC.p1 quadOC.p2 quadOC.p3];
        predquadOC = quadmodelOC(1)*OCDOrespS1{jj}.^2 + quadmodelOC(2)*OCDOrespS1{jj} + quadmodelOC(3)*ones(size(OCDOrespS1{jj}));
        SSE_quadmodelOC = [SSE_quadmodelOC; sum((predquadOC-OCDOrespS2{jj}).^2)];
        
        
        [THETALM,~] = cart2pol(LMDOrespS1{jj},LMDOrespS2{jj});
        [THETAOC,~] = cart2pol(OCDOrespS1{jj},OCDOrespS2{jj});
%         THETALM(THETALM<0) = 2*pi-THETALM(THETALM<0);
%         THETAOC(THETAOC<0) = 2*pi-THETAOC(THETAOC<0);
        
        hcsurfaces(jj,1) = ~(all(sort(THETALM) == THETALM) | all(sort(THETALM) == fliplr(THETALM)));
        hcsurfaces(jj,2) = ~(all(sort(THETAOC) == THETAOC) | all(sort(THETAOC) == fliplr(THETAOC)));
       
%         figure(plot_counter);subplot(121);plot(LMDOrespS1{jj},LMDOrespS2{jj},'o'); hold on; plot(LMDOrespS1{jj},predlinLM,'g','Linewidth',2); plot(LMDOrespS1{jj},predquadLM,'b','Linewidth',2);hold off;
%         subplot(122);plot(OCDOrespS1{jj},OCDOrespS2{jj},'o'); hold on; plot(OCDOrespS1{jj},predlinOC,'g','Linewidth',2);plot(OCDOrespS1{jj},predquadOC,'b','Linewidth',2); hold off;
%         plot_counter = plot_counter + 1;
        
    end
    rratioLM = SSE_linearmodelLM./SSE_quadmodelLM;
    rratioLM = rratioLM(imag(rratioLM)==0);
    rratioLMtmp = log10(rratioLM);
    rratioLMtmp(logical(hcsurfaces(:,1)))= 1000;
    
    rratioOC = SSE_linearmodelOC./SSE_quadmodelOC;
    rratioOC = rratioOC(imag(rratioOC)==0);
    rratioOCtmp = log10(rratioOC);
    rratioOCtmp(logical(hcsurfaces(:,2)))= 1000;
    
    meandiff(aa) = mean(rratioLMtmp) - mean(rratioOCtmp);
    binsizes = 0:0.5:3;
    [cLM,~] = hist(rratioLMtmp,binsizes);
    [cOC,~] = hist(rratioOCtmp,binsizes);
    figure(plot_counter);bar(binsizes,[cOC' cLM']); legend('OC','LM');
    plot_counter = plot_counter + 1;
    
end

figure(plot_counter);
for ii = 1:1000
    subplot(121);plot(LMDOrespS1{ii},LMDOrespS2{ii},'r'); hold on;
    subplot(122);plot(OCDOrespS2{ii},OCDOrespS1{ii},'b'); hold on;
end
subplot(121); axis equal; axis square; grid on; title('Lime-Magenta'); hold off;
subplot(122); axis equal; axis square; grid on; title('Orange-Cyan');hold off;
plot_counter = plot_counter + 1;
