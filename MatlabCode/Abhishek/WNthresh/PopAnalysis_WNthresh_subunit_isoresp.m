% A new script for find any correlation between the signal integration within subunits, between the subunits and the isoresponse curve
% Author - Abhishek De, 10/18
close all; clearvars;
% load all the variables from the subunits data (within subunit interaction from the whitenoise subunit data)
load Fprob1.mat
load Fprob2.mat
load diffmodelpred1.mat
load diffmodelpred2.mat
load LLRprob1.mat
load LLRprob2.mat
load Devprob1.mat
load Devprob2.mat
load sigquadterms1.mat
load sigquadterms2.mat
load AIClin1.mat 
load AICquad1.mat
load AIClin2.mat
load AICquad2.mat
load BIClin1.mat
load BICquad1.mat
load BIClin2.mat
load BICquad2.mat
load mSSElin1.mat
load mSSEquad1.mat
load mSSElin2.mat
load mSSEquad2.mat 

% load all the variables from the subunits data (between subunit interaction from the whitenoise subunit data)
load Fprob.mat
load diffmodelpred.mat
load LLRprob.mat
load Devprob.mat    
load sigquadterms.mat
load AIClin.mat
load AICquad.mat
load BIClin.mat
load BICquad.mat
load mSSElin.mat
load mSSEquad.mat 

% load all the variables from the isoresp fits data
load filename_c.mat % Color cells
load filename_l.mat % Luminance cells
load filenameSTAvsPC1.mat
load filenameSTAsubunits.mat
load S1LMS.mat % cone wts for subunit 1
load S2LMS.mat % cone wts for subunit 2
load linear_modelparams.mat % parameters from the linear fit
load quad_modelparams.mat % parameters from the quadratic fit
load RSSE_linearmodel.mat % errors from robust regression (linear fit)
load RSSE_quadmodel.mat % errors from robust regression (quadratic fit)
load SSE_linearmodel.mat % errors from least square regression (linear fit)
load SSE_quadmodel.mat % errors from least square regression (quadratic fit)
load pixelswithinsubunit1.mat % size of subunit 1 in terms of pixels 
load pixelswithinsubunit2.mat % size of subunit 2 in terms of pixels
load ratioeig.mat 
load RHO_all.mat
load THETA_all.mat
load oog_idx_all.mat
load not_oog_idx_all.mat
absS = abs(S1LMS(3,:)) + abs(S2LMS(3,:));
absSsum = abs(S1LMS(3,:)+S2LMS(3,:));
absSdiff = abs(S1LMS(3,:)-S2LMS(3,:));
absLminusM = abs(S1LMS(1,:)-S1LMS(2,:)) + abs(S2LMS(1,:)-S2LMS(2,:));
absLplusM = abs(S1LMS(1,:)+S1LMS(2,:)) + abs(S2LMS(1,:)+S2LMS(2,:));

% load all the variables from the FR surface fitting 
load linear_FRsurfacefitLL.mat % negative log likelihood from fitting linear model to the FR surface 
load quad_FRsurfacefitLL.mat % negative log likelihood from fitting quadratic model to the FR surface
load diffmodelpredFR.mat % mean square difference between the model predictions

% load data from the Gabor and DOG fits (fitting the RFs with 2 models) 
load DOGBIC.mat
load DOGAIC.mat
load GaborBIC.mat
load GaborAIC.mat
load RFstructure_c.mat
load RFstructure_l.mat
deltaBIC = DOGBIC-GaborBIC;
deltaAIC = DOGAIC-GaborAIC;
RFstructure = [RFstructure_c; RFstructure_l];

% load other variables 
load PCS1RGB.mat
load PCS2RGB.mat
load anglediffWnchecksubunit.mat
load z_scores.mat
load anglebwvectorsRGB.mat
load pruns_linearmodel.mat
load NLindex_r.mat % data driven non-linearity index
load NLindex_p.mat % p-value
load baselineFRstats.mat
load latencySTAcheck.mat
load latencySTAsubunit.mat
load latencybetweensubunits.mat
load latencydiffWNchecksubunit.mat
load numsubunitspikes.mat
load pranksumOFF.mat
load RF_loc.mat
load Pangu.mat
load Maui.mat
RF_loc = sqrt(RF_loc(:,1).^2 + RF_loc(:,2).^2);
Monkey = zeros(numel([Pangu; Maui]),1);
Monkey(Maui) = 1; 
Monkey(Pangu) = 2;

% loading the indices for DO, SO, LUM and hardtoclassifycells
load newOCidx.mat
load newLMidx.mat
load newLUMidx.mat
load newSOidx.mat
load newhardtoclassifyidx.mat
filename = [filename_c; filename_l];
OCidx = newOCidx;
LMidx = newLMidx;
LUMidx = newLUMidx';
SOidx = newSOidx;
hardtoclassifyidx = newhardtoclassifyidx;
idx = [OCidx LMidx LUMidx SOidx hardtoclassifyidx];
% RFstructure = RFstructure(idx);
% RF_loc = RF_loc(idx);
% absS = absS(idx);
% RHO_all = RHO_all(:,idx);
% THETA_all = THETA_all(:,[OCidx LMidx LUMidx SOidx hardtoclassifyidx]);
% oog_idx_all = oog_idx_all(:,[OCidx LMidx LUMidx SOidx hardtoclassifyidx]);
% not_oog_idx_all = not_oog_idx_all(:,[OCidx LMidx LUMidx SOidx hardtoclassifyidx]);
% tmp = [numel(OCidx) numel(LMidx) numel(LUMidx) numel(SOidx) numel(hardtoclassifyidx)];
% tmp = [0 cumsum(tmp)];
% OCidx = tmp(1)+1:tmp(2);
% LMidx = tmp(2)+1:tmp(3);
% LUMidx = tmp(3)+1:tmp(4);
% SOidx = tmp(4)+1:tmp(5);
% hardtoclassifyidx = tmp(5)+1:tmp(6);
% idx = [OCidx LMidx LUMidx SOidx hardtoclassifyidx];

% idx = [OCidx LMidx]; % Only DO cells
% idx = [LUMidx]; % Only luminance cells
% idx =  [OCidx LMidx LUMidx]; % DO & luminance cells
% idx = [hardtoclassifyidx]; % hard to classify cells
% idx = [SOidx]; % SO cells
%
% Fprob1 = Fprob1(idx); Fprob2 = Fprob2(idx);
% diffmodelpred1 = diffmodelpred1(idx); diffmodelpred2 = diffmodelpred2(idx);
% LLRprob1 = LLRprob1(idx); LLRprob2 = LLRprob2(idx);
% Devprob1 = Devprob1(idx); Devprob2 = Devprob2(idx);
% sigquadterms1 = sigquadterms1(idx); 
% sigquadterms2 = sigquadterms2(idx);
% AIClin1 = AIClin1(idx); AICquad1 = AICquad1(idx); AIClin2 = AIClin2(idx); AICquad2 = AICquad2(idx);
% BIClin1 = BIClin1(idx); BICquad1 = BICquad1(idx); BIClin2 = BIClin2(idx); BICquad2 = BICquad2(idx);
% mSSElin1 = mSSElin1(idx); mSSEquad1 = mSSEquad1(idx); mSSElin2 = mSSElin2(idx); mSSEquad2 = mSSEquad2(idx); 
% 
% Fprob = Fprob(idx); diffmodelpred = diffmodelpred(idx); LLRprob = LLRprob(idx); Devprob = Devprob(idx);    
% sigquadterms = sigquadterms(idx); AIClin = AIClin(idx); AICquad = AICquad(idx); 
% BIClin = BIClin(idx); BICquad = BICquad(idx); mSSElin = mSSElin(idx); mSSEquad = mSSEquad(idx);
% 
% linear_modelparams = linear_modelparams(idx); quad_modelparams = quad_modelparams(idx);
% RSSE_linearmodel = RSSE_linearmodel(idx); RSSE_quadmodel = RSSE_quadmodel(idx);
% SSE_linearmodel = SSE_linearmodel(idx); SSE_quadmodel = SSE_quadmodel(idx);
% pixelswithinsubunit1 = pixelswithinsubunit1(idx); pixelswithinsubunit2 = pixelswithinsubunit2(idx); ratioeig = ratioeig(idx);
% 
% linear_FRsurfacefitLL = linear_FRsurfacefitLL(idx); quad_FRsurfacefitLL = quad_FRsurfacefitLL(idx); diffmodelpredFR = diffmodelpredFR(idx);
% DOGBIC = DOGBIC(idx); DOGAIC = DOGAIC(idx); GaborBIC = GaborBIC(idx); GaborAIC = GaborAIC(idx);
% deltaBIC = deltaBIC(idx); deltaAIC = deltaAIC(idx); RFstructure = RFstructure(idx);
% 
% anglediffWNchecksubunit = anglediffWNchecksubunit(idx); z_scores = z_scores(idx); anglebwvectorsRGB = anglebwvectorsRGB(idx);
% pruns_linearmodel = pruns_linearmodel(idx); NLindex_p = NLindex_p(idx); NLindex_r = NLindex_r(idx); baselineFRstats = baselineFRstats(idx);
% latencySTAcheck = latencySTAcheck(idx); latencySTAsubunit = latencySTAsubunit(idx); 
% latencybetweensubunits = latencybetweensubunits(idx); latencydiffWNchecksubunit = latencydiffWNchecksubunit(idx);
% absS = absS(idx); numsubunitspikes = numsubunitspikes(idx); RF_loc = RF_loc(idx); pranksumOFF = pranksumOFF(idx);
% RHO_all = RHO_all(:,idx); THETA_all = THETA_all(:,idx); oog_idx_all = oog_idx_all(:,idx); not_oog_idx_all = not_oog_idx_all(:,idx);

Q_sse = log(SSE_linearmodel./SSE_quadmodel);
[~,ind_sse] = sort(Q_sse);
Q_rob = log(RSSE_linearmodel./RSSE_quadmodel);
[~,ind_rob] = sort(Q_rob);

% Plotting results 
plot_counter = 1;
for ii = 1:numel(filename)
    ind1 = ind_rob(ii);
    ind2 = ind_sse(ii);
    
    % Finding any trend between the isoresponse data and the signal integration within subunits
    figure(plot_counter);
    subplot(321); plot([ii ii],[diffmodelpred1(ind1) diffmodelpred2(ind1)],'r'); hold on;
    subplot(322); plot([ii ii],[diffmodelpred1(ind2) diffmodelpred2(ind2)],'b'); hold on;
    subplot(323); plot([ii ii],[sigquadterms1(ind1) sigquadterms2(ind1)],'r'); hold on;
    subplot(324); plot([ii ii],[sigquadterms1(ind2) sigquadterms2(ind2)],'b'); hold on;
    subplot(325); plot([ii ii],[mSSElin1(ind1)-mSSEquad1(ind1) mSSElin2(ind1)-mSSEquad2(ind1)],'r'); hold on;
    subplot(326); plot([ii ii],[mSSElin1(ind2)-mSSEquad1(ind2) mSSElin2(ind2)-mSSEquad2(ind2)],'b'); hold on;
    
    figure(plot_counter+1);
    subplot(321); plot([ii ii],[AIClin1(ind1)-AICquad1(ind1) AIClin2(ind1)-AICquad2(ind1)],'r'); hold on;
    subplot(322); plot([ii ii],[AIClin1(ind2)-AICquad1(ind2) AIClin2(ind2)-AICquad2(ind2)],'b'); hold on;
    subplot(323); plot([ii ii],[BIClin1(ind1)-BICquad1(ind1) BIClin2(ind1)-BICquad2(ind1)],'r'); hold on;
    subplot(324); plot([ii ii],[BIClin1(ind2)-BICquad1(ind2) BIClin2(ind2)-BICquad2(ind2)],'b'); hold on;
    subplot(325); plot([ii ii],[LLRprob1(ind1) LLRprob2(ind1)],'r'); hold on;
    subplot(326); plot([ii ii],[LLRprob1(ind2) LLRprob2(ind2)],'b'); hold on;
    
    % Finding any trend between the isoresponse data and the fits to the FR surface and the size of the subunits
    figure(plot_counter+2);
    subplot(321); plot([ii ii],[pixelswithinsubunit1(ind1) pixelswithinsubunit2(ind1)],'r'); hold on;
    subplot(322); plot([ii ii],[pixelswithinsubunit1(ind2) pixelswithinsubunit2(ind2)],'b'); hold on;
    subplot(323); plot(ii,linear_FRsurfacefitLL(ind1)-quad_FRsurfacefitLL(ind1),'ro','MarkerFacecolor',[1 0 0],'MarkerEdgecolor',[1 1 1]); hold on;
    subplot(324); plot(ii,linear_FRsurfacefitLL(ind2)-quad_FRsurfacefitLL(ind2),'bo','MarkerFacecolor',[0 0 1],'MarkerEdgecolor',[1 1 1]); hold on;
    subplot(325); plot(ii,diffmodelpredFR(ind1),'ro','MarkerFacecolor',[1 0 0],'MarkerEdgecolor',[1 1 1]); hold on;
    subplot(326); plot(ii,diffmodelpredFR(ind2),'bo','MarkerFacecolor',[0 0 1],'MarkerEdgecolor',[1 1 1]); hold on;
    
    figure(plot_counter+3);
    subplot(321);plot(ii,diffmodelpred(ind1),'ro','MarkerFacecolor',[1 0 0],'MarkerEdgecolor',[1 1 1]); hold on; 
    subplot(322);plot(ii,diffmodelpred(ind2),'bo','MarkerFacecolor',[0 0 1],'MarkerEdgecolor',[1 1 1]); hold on; 
    subplot(323);plot(ii,sigquadterms(ind1),'ro','MarkerFacecolor',[1 0 0],'MarkerEdgecolor',[1 1 1]); hold on; 
    subplot(324);plot(ii,sigquadterms(ind2),'bo','MarkerFacecolor',[0 0 1],'MarkerEdgecolor',[1 1 1]); hold on; 
    subplot(325);plot(ii,mSSElin(ind1)-mSSEquad(ind1),'ro','MarkerFacecolor',[1 0 0],'MarkerEdgecolor',[1 1 1]); hold on; 
    subplot(326);plot(ii,mSSElin(ind2)-mSSEquad(ind2),'bo','MarkerFacecolor',[0 0 1],'MarkerEdgecolor',[1 1 1]); hold on; 
end

[r1,p1] = corr([diffmodelpred1;diffmodelpred2],[Q_rob;Q_rob],'Type','Spearman');
[r2,p2] = corr([sigquadterms1;sigquadterms2],[Q_rob;Q_rob],'Type','Spearman');
[r3,p3] = corr([mSSElin1-mSSEquad1;mSSElin2-mSSEquad2],[Q_rob;Q_rob],'Type','Spearman');
figure(plot_counter); set(gcf,'Name','Subunit non-linearity vs isoresp non-linearity I');
subplot(321); set(gca,'Yscale','log'); xlabel('Q rob ascending order'); ylabel('diff in model pred'); title(strcat(num2str(r1,2),',',' ',num2str(p1,2))); hold off; 
subplot(322); set(gca,'Yscale','log'); xlabel('Q sse ascending order'); ylabel('diff in model pred'); hold off;
subplot(323); xlabel('Q rob ascending order'); ylabel('sig quad terms'); title(strcat(num2str(r2,2),',',' ',num2str(p2,2))); hold off; 
subplot(324); xlabel('Q sse ascending order'); ylabel('sig quad terms'); hold off;
subplot(325); xlabel('Q rob ascending order'); ylabel('mean squared error diff'); title(strcat(num2str(r3,2),',',' ',num2str(p3,2))); hold off; 
subplot(326); xlabel('Q sse ascending order'); ylabel('mean squared error diff'); hold off;

[r1,p1] = corr([AIClin1-AICquad1;AIClin2-AICquad2],[Q_rob;Q_rob],'Type','Spearman');
[r2,p2] = corr([BIClin1-BICquad1;BIClin2-BICquad2],[Q_rob;Q_rob],'Type','Spearman');
figure(plot_counter+1); set(gcf,'Name','Subunit non-linearity vs isoresp non-linearity II');
subplot(321); xlabel('Q rob ascending order'); ylabel('diff in AIC'); title(strcat(num2str(r1,2),',',' ',num2str(p1,2))); hold off; 
subplot(322); xlabel('Q sse ascending order'); ylabel('diff in AIC'); hold off;
subplot(323); xlabel('Q rob ascending order'); ylabel('diff in BIC'); title(strcat(num2str(r2,2),',',' ',num2str(p2,2))); hold off; 
subplot(324); xlabel('Q sse ascending order'); ylabel('diff in BIC'); hold off;
subplot(325); xlabel('Q rob ascending order'); ylabel('LLR prob'); hold off; 
subplot(326); xlabel('Q sse ascending order'); ylabel('LLR prob'); hold off;

[r1,p1] = corr([pixelswithinsubunit1;pixelswithinsubunit1],[Q_rob;Q_rob],'Type','Spearman');
[r2,p2] = corr(linear_FRsurfacefitLL-quad_FRsurfacefitLL,Q_rob,'Type','Spearman');
[r3,p3] = corr(diffmodelpredFR,Q_rob,'Type','Spearman');
figure(plot_counter+2); set(gcf,'Name','Isoresp contour vs FR surface fits');
subplot(321); set(gca,'Yscale','log'); xlabel('Q rob ascending order'); ylabel('subunits size'); hold off; title(strcat(num2str(r1,2),',',' ',num2str(p1,2)));hold off;
subplot(322); set(gca,'Yscale','log'); xlabel('Q sse ascending order'); ylabel('subunits size'); hold off;
subplot(323); set(gca,'Yscale','log'); xlabel('Q rob ascending order'); ylabel('diff in LL fits'); title(strcat(num2str(r2,2),',',' ',num2str(p2,2)));hold off;
subplot(324); set(gca,'Yscale','log'); xlabel('Q sse ascending order'); ylabel('diff in LL fits'); hold off;
subplot(325); set(gca,'Yscale','log'); xlabel('Q rob ascending order'); ylabel('diff in model pred'); title(strcat(num2str(r3,2),',',' ',num2str(p3,2)));hold off;
subplot(326); set(gca,'Yscale','log'); xlabel('Q sse ascending order'); ylabel('diff in model pred'); hold off;

[r1,p1] = corr(diffmodelpred,Q_rob,'Type','Spearman');
[r2,p2] = corr(sigquadterms,Q_rob,'Type','Spearman');
[r3,p3] = corr(mSSElin-mSSEquad,Q_rob,'Type','Spearman');
figure(plot_counter+3); set(gcf,'Name','Isoresp contour vs across subunits');
subplot(321); set(gca,'Yscale','log'); xlabel('Q rob ascending order'); ylabel('diff in model pred'); title(strcat(num2str(r1,2),',',' ',num2str(p1,2))); hold off;
subplot(322); set(gca,'Yscale','log'); xlabel('Q sse ascending order'); ylabel('diff in model pred'); hold off;
subplot(323); xlabel('Q rob ascending order'); ylabel('sig quad terms'); title(strcat(num2str(r2,2),',',' ',num2str(p2,2))); hold off;
subplot(324); xlabel('Q sse ascending order'); ylabel('sig quad terms'); hold off;
subplot(325); xlabel('Q rob ascending order'); ylabel('mean squared error diff'); title(strcat(num2str(r3,2),',',' ',num2str(p3,2))); hold off;
subplot(326); xlabel('Q sse ascending order'); ylabel('mean squared error diff'); hold off;
plot_counter = plot_counter + 4;

% Figure for plotting the relationship between color angle between the checkerboard and subunit stimuli to Q rob and Q sse 
[r1,p1] = corr(anglediffWNchecksubunit,Q_rob,'Type','Spearman');
[r2,p2] = corr(z_scores,Q_rob,'Type','Spearman');
[r3,p3] = corr(anglebwvectorsRGB,Q_rob,'Type','Spearman');
figure(plot_counter); set(gcf,'Name','Isoresp vs measuresments from WN I');
subplot(321); plot(anglediffWNchecksubunit(ind_rob),'ro','MarkerFacecolor',[1 0 0],'MarkerEdgecolor',[1 1 1]); set(gca,'Yscale','log'); xlabel('Q rob ascending order'); ylabel('angle diff check sub'); title(strcat(num2str(r1,2),',',' ',num2str(p1,2)));
subplot(322); plot(anglediffWNchecksubunit(ind_sse),'bo','MarkerFacecolor',[0 0 1],'MarkerEdgecolor',[1 1 1]); set(gca,'Yscale','log'); xlabel('Q sse ascending order'); ylabel('angle diff check sub');
subplot(323); plot(abs(z_scores(ind_rob)),'ro','MarkerFacecolor',[1 0 0],'MarkerEdgecolor',[1 1 1]); xlabel('Q rob ascending order'); ylabel('PC1 zscore');title(strcat(num2str(r2,2),',',' ',num2str(p2,2)));
subplot(324); plot(abs(z_scores(ind_sse)),'bo','MarkerFacecolor',[0 0 1],'MarkerEdgecolor',[1 1 1]); xlabel('Q sse ascending order'); ylabel('PC1 zscore');
subplot(325); plot(anglebwvectorsRGB(ind_rob),'ro','MarkerFacecolor',[1 0 0],'MarkerEdgecolor',[1 1 1]); xlabel('Q rob ascending order'); ylabel('angle diff WN sub'); title(strcat(num2str(r3,2),',',' ',num2str(p3,2)));
subplot(326); plot(anglebwvectorsRGB(ind_sse),'bo','MarkerFacecolor',[0 0 1],'MarkerEdgecolor',[1 1 1]); xlabel('Q sse ascending order'); ylabel('angle diff WN sub'); 
plot_counter = plot_counter + 1;


[r1,p1] = corr(latencySTAcheck,Q_rob,'Type','Spearman');
[r2,p2] = corr(latencySTAsubunit,Q_rob,'Type','Spearman');
[r3,p3] = corr(latencybetweensubunits,Q_rob,'Type','Spearman');
[r4,p4] = corr(latencydiffWNchecksubunit,Q_rob,'Type','Spearman');
figure(plot_counter); set(gcf,'Name','Isoresp vs measuresments from WN II');
subplot(421); plot(latencySTAcheck(ind_rob),'ro','MarkerFacecolor',[1 0 0],'MarkerEdgecolor',[1 1 1]); xlabel('Q rob ascending order'); ylabel('latency check'); title(strcat(num2str(r1,2),',',' ',num2str(p1,2)));
subplot(422); plot(latencySTAcheck(ind_sse),'bo','MarkerFacecolor',[0 0 1],'MarkerEdgecolor',[1 1 1]); xlabel('Q sse ascending order'); ylabel('latency check');
subplot(423); plot(latencySTAsubunit(ind_rob),'ro','MarkerFacecolor',[1 0 0],'MarkerEdgecolor',[1 1 1]); xlabel('Q rob ascending order'); ylabel('latency sub');title(strcat(num2str(r2,2),',',' ',num2str(p2,2)));
subplot(424); plot(latencySTAsubunit(ind_sse),'bo','MarkerFacecolor',[0 0 1],'MarkerEdgecolor',[1 1 1]); xlabel('Q sse ascending order'); ylabel('latency sub');
subplot(425); plot(latencybetweensubunits(ind_rob),'ro','MarkerFacecolor',[1 0 0],'MarkerEdgecolor',[1 1 1]); xlabel('Q rob ascending order'); ylabel('latency bw sub'); title(strcat(num2str(r3,2),',',' ',num2str(p3,2)));
subplot(426); plot(latencybetweensubunits(ind_sse),'bo','MarkerFacecolor',[0 0 1],'MarkerEdgecolor',[1 1 1]); xlabel('Q sse ascending order'); ylabel('latency bw sub'); 
subplot(427); plot(latencydiffWNchecksubunit(ind_rob),'ro','MarkerFacecolor',[1 0 0],'MarkerEdgecolor',[1 1 1]); xlabel('Q rob ascending order'); ylabel('latency bw sub check'); title(strcat(num2str(r4,2),',',' ',num2str(p4,2)));
subplot(428); plot(latencydiffWNchecksubunit(ind_sse),'bo','MarkerFacecolor',[0 0 1],'MarkerEdgecolor',[1 1 1]); xlabel('Q sse ascending order'); ylabel('latency bw sub check');
plot_counter = plot_counter + 1;

% Looking at some other trends in the data 
[r1,p1] = corr(NLindex_r,Q_rob,'Type','Spearman');
[r2,p2] = corr(pruns_linearmodel,Q_rob,'Type','Spearman');
[r3,p3] = corr(ratioeig,Q_rob,'Type','Spearman');
figure(plot_counter); subplot(321); plot(abs(NLindex_r(ind_rob)),'ro','MarkerFacecolor',[1 0 0],'MarkerEdgecolor',[1 1 1]); xlabel('Q rob ascending order'); ylabel('NL index r'); title(strcat(num2str(r1,2),',',' ',num2str(p1,2)));
subplot(322); plot(abs(NLindex_r(ind_sse)),'bo','MarkerFacecolor',[0 0 1],'MarkerEdgecolor',[1 1 1]); xlabel('Q rob ascending order'); ylabel('NL index r'); 
subplot(323); plot(pruns_linearmodel(ind_rob),'ro','MarkerFacecolor',[1 0 0],'MarkerEdgecolor',[1 1 1]); set(gca,'Yscale','log'); xlabel('Q rob ascending order'); ylabel('p runs test lin'); title(strcat(num2str(r2,2),',',' ',num2str(p2,2)));
subplot(324); plot(pruns_linearmodel(ind_sse),'bo','MarkerFacecolor',[0 0 1],'MarkerEdgecolor',[1 1 1]); set(gca,'Yscale','log'); xlabel('Q sse ascending order'); ylabel('p runs test lin'); 
subplot(325); plot(ratioeig(ind_rob),'ro','MarkerFacecolor',[1 0 0],'MarkerEdgecolor',[1 1 1]); set(gca,'Yscale','log'); xlabel('Q rob ascending order'); ylabel('ratio eig'); title(strcat(num2str(r3,2),',',' ',num2str(p3,2)));
subplot(326); plot(ratioeig(ind_sse),'bo','MarkerFacecolor',[0 0 1],'MarkerEdgecolor',[1 1 1]); set(gca,'Yscale','log'); xlabel('Q sse ascending order'); ylabel('ratio eig'); 
plot_counter = plot_counter + 1;

% Looking if how strongly correlated is the model based and the user based approach for determining the RF structure
[h_fisher,p_fisher,stats_fisher] = fishertest([sum(deltaBIC>=0 & RFstructure==1), sum(deltaBIC>=0 & RFstructure==2); sum(deltaBIC<0 & RFstructure==1), sum(deltaBIC<0 & RFstructure==2)]);
mdl = fitglm(deltaBIC,RFstructure-1,'linear','Distribution','binomial','Link','logit');
[h_fisher1,p_fisher1,stats_fisher1] = fishertest([sum(deltaAIC>=0 & RFstructure==1), sum(deltaAIC>=0 & RFstructure==2); sum(deltaAIC<0 & RFstructure==1), sum(deltaAIC<0 & RFstructure==2)]);
mdl1 = fitglm(deltaAIC,RFstructure-1,'linear','Distribution','binomial','Link','logit');
bins = -50:10:300;
figure(plot_counter); set(gcf,'Name','RF structure vs delta BIC');
subplot(231);plot(deltaBIC,RFstructure-1,'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1]); xlabel('delta BIC'); hold on;
plot([-100:1:300],predict(mdl,[-100:1:300]'),'Linewidth',2); axis square; text(100,0.5,strcat('p=',num2str(p_fisher))); hold off;
subplot(232); plot(GaborBIC,DOGBIC,'o','MarkerSize',6,'LineWidth',0.5,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[1 1 1]); hold on; axis square;
set(gca,'Xlim',[-1000 -100],'Ylim',[-1000 -100],'TickDir','out'); line([-1000 -100],[-1000 -100]);xlabel('Gabor BIC'); ylabel('DOG BIC'); title('BIC'); hold off;
subplot(233); histogram(deltaBIC,bins); axis square; xlabel('delta BIC'); ylabel('count'); set(gca,'TickDir','out');hold off;
subplot(234);plot(deltaAIC,RFstructure-1,'o','MarkerSize',5,'LineWidth',0.5,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1]); xlabel('delta AIC'); hold on;
plot([-100:1:300],predict(mdl1,[-100:1:300]'),'Linewidth',2); axis square; text(100,0.5,strcat('p=',num2str(p_fisher1))); hold off;
subplot(235); plot(GaborAIC,DOGAIC,'o','MarkerSize',6,'LineWidth',0.5,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[1 1 1]); hold on; axis square;
set(gca,'Xlim',[-1000 -100],'Ylim',[-1000 -100],'TickDir','out'); line([-1000 -100],[-1000 -100]);xlabel('Gabor AIC'); ylabel('DOG AIC'); title('AIC'); hold off;
subplot(236); histogram(deltaAIC,bins); axis square; xlabel('delta AIC'); ylabel('count'); set(gca,'TickDir','out');hold off;
plot_counter = plot_counter + 1;

% Comparing delta BICs for different cell classes
bins = -50:10:300;
figure(plot_counter); subplot(221); histogram(deltaBIC,bins); hold on; histogram(deltaBIC([OCidx LMidx]),bins); plot(median(deltaBIC),0,'kv','Markerfacecolor',[0 0 1]); plot(median(deltaBIC([OCidx LMidx])),0,'kv','Markerfacecolor',[1 0 0]); axis square; xlabel('delta BIC'); ylabel('count'); title('DO'); set(gca,'TickDir','out');hold off;
subplot(222); histogram(deltaBIC,bins); hold on; histogram(deltaBIC(LUMidx),bins); plot(median(deltaBIC),0,'kv','Markerfacecolor',[0 0 1]); plot(median(deltaBIC([LUMidx])),0,'kv','Markerfacecolor',[1 0 0]); axis square; xlabel('delta BIC'); ylabel('count'); title('LUM'); set(gca,'TickDir','out');hold off;
subplot(223); histogram(deltaBIC,bins); hold on; histogram(deltaBIC(hardtoclassifyidx),bins); plot(median(deltaBIC),0,'kv','Markerfacecolor',[0 0 1]); plot(median(deltaBIC([hardtoclassifyidx])),0,'kv','Markerfacecolor',[1 0 0]); axis square; xlabel('delta BIC'); ylabel('count'); title('htc');  set(gca,'TickDir','out');hold off;
subplot(224); histogram(deltaBIC,bins); hold on; histogram(deltaBIC(SOidx),bins); plot(median(deltaBIC),0,'kv','Markerfacecolor',[0 0 1]); plot(median(deltaBIC(SOidx)),0,'kv','Markerfacecolor',[1 0 0]); axis square; xlabel('delta BIC'); ylabel('count'); title('SO'); set(gca,'TickDir','out');hold off;
plot_counter = plot_counter + 1;

% Plotting the box plot version of the previous analyses
x = [ones(numel([OCidx LMidx]),1); 2*ones(numel([LUMidx]),1); 3*ones(numel([SOidx]),1); 4*ones(numel([hardtoclassifyidx]),1)];
g1 = deltaBIC([OCidx LMidx LUMidx SOidx hardtoclassifyidx]); p1 = anova1(g1,x,'off');
g2 = Q_rob([OCidx LMidx LUMidx SOidx hardtoclassifyidx]); p2 = anova1(g2,x,'off');
g3 = absS([OCidx LMidx LUMidx SOidx hardtoclassifyidx]); p3 = anova1(g3,x,'off');
g4 = z_scores([OCidx LMidx LUMidx SOidx hardtoclassifyidx]); p4 = anova1(g4,x,'off');
g5 = anglebwvectorsRGB([OCidx LMidx LUMidx SOidx hardtoclassifyidx]); p5 = anova1(g5,x,'off');
g6 = anglediffWNchecksubunit([OCidx LMidx LUMidx SOidx hardtoclassifyidx]); p6 = anova1(g6,x,'off');
g7 = NLindex_r([OCidx LMidx LUMidx SOidx hardtoclassifyidx]); p7 = anova1(g7,x,'off');
g8 = ratioeig([OCidx LMidx LUMidx SOidx hardtoclassifyidx]); p8 = anova1(g8,x,'off');
g9 = absSsum([OCidx LMidx LUMidx SOidx hardtoclassifyidx]); p9 = anova1(g9,x,'off');
g10 = absSdiff([OCidx LMidx LUMidx SOidx hardtoclassifyidx]); p10 = anova1(g10,x,'off');
figure(plot_counter); set(gcf,'Name','Boxplots for DO,LUM,SO and htc: ANOVA 1 way');
subplot(251); boxplot(g1,x,'Notch','on','Labels',{'DO','LUM','SO','htc'}); title(strcat('deltaBIC',{' '},num2str(p1,3)));
subplot(252); boxplot(g2,x,'Notch','on','Labels',{'DO','LUM','SO','htc'}); title(strcat('Q rob',{' '},num2str(p2,3)));
subplot(253); boxplot(g3,x,'Notch','on','Labels',{'DO','LUM','SO','htc'}); title(strcat('abs S',{' '},num2str(p3,3)));
subplot(254); boxplot(g4,x,'Notch','on','Labels',{'DO','LUM','SO','htc'}); title(strcat('z score',{' '},num2str(p4,3)));
subplot(255); boxplot(g5,x,'Notch','on','Labels',{'DO','LUM','SO','htc'}); title(strcat('angle bt sub',{' '},num2str(p5,3)));
subplot(256); boxplot(g6,x,'Notch','on','Labels',{'DO','LUM','SO','htc'}); title(strcat('angle bt check sub',{' '},num2str(p6,3)));
subplot(257); boxplot(g7,x,'Notch','on','Labels',{'DO','LUM','SO','htc'}); title(strcat('NL index r',{' '},num2str(p7,3)));
subplot(258); boxplot(g8,x,'Notch','on','Labels',{'DO','LUM','SO','htc'}); title(strcat('% PC1: ratioeig',{' '},num2str(p8,3)));
subplot(259); boxplot(g9,x,'Notch','on','Labels',{'DO','LUM','SO','htc'}); title(strcat('abs S sum',{' '},num2str(p9,3)));
subplot(2,5,10); boxplot(g10,x,'Notch','on','Labels',{'DO','LUM','SO','htc'}); title(strcat('abs S diff',{' '},num2str(p10,3)));
plot_counter = plot_counter + 1;

bins = 0:0.5:10;
figure(plot_counter); subplot(231); boxplot(Q_rob,RFstructure,'Notch','on','Labels',{'CS','Adj'}); title('Q rob');
subplot(232); boxplot(Q_rob,sign(deltaBIC),'Notch','on','Labels',{'DOG','Gabor'}); title('Q rob');
subplot(233); boxplot(Q_rob,Monkey,'Notch','on','Labels',{'Maui','Pangu'}); title('Q rob');
subplot(234); histogram(Q_rob(RFstructure==1),bins); hold on; histogram(Q_rob(RFstructure==2),bins); plot(median(Q_rob(RFstructure==1)),0,'kv','Markerfacecolor',[0 0 1]); plot(median(Q_rob(RFstructure==2)),0,'kv','Markerfacecolor',[1 0 0]); title('Q_rob'); legend('CS','Adj');
subplot(235); histogram(Q_rob(sign(deltaBIC)<0),bins); hold on; histogram(Q_rob(sign(deltaBIC)>0),bins); plot(median(Q_rob(sign(deltaBIC)<0)),0,'kv','Markerfacecolor',[0 0 1]); plot(median(Q_rob(sign(deltaBIC)>0)),0,'kv','Markerfacecolor',[1 0 0]); title('Q_rob'); legend('DOG','Gabor');
subplot(236); histogram(Q_rob(Maui),bins); hold on; histogram(Q_rob(Pangu),bins); plot(median(Q_rob(Maui)),0,'kv','Markerfacecolor',[0 0 1]); plot(median(Q_rob(Pangu)),0,'kv','Markerfacecolor',[1 0 0]);  title('Q_rob'); legend('Maui','Pangu');
plot_counter = plot_counter + 1;

% Looking at the correlation between absolute S cone-signal and the iso-response data
[r1,p1] = corr(absS',Q_rob,'Type','Spearman');
[r2,p2] = corr(numsubunitspikes,Q_rob,'Type','Spearman');
[r3,p3] = corr(RF_loc,Q_rob,'Type','Spearman');
figure(plot_counter); subplot(321), plot(absS(ind_rob),'ro','MarkerFacecolor',[1 0 0],'MarkerEdgecolor',[1 1 1]); xlabel('Q rob ascending order'); ylabel('S'); title(strcat(num2str(r1,2),',',' ',num2str(p1,2)));
subplot(322), plot(absS(ind_sse),'bo','MarkerFacecolor',[0 0 1],'MarkerEdgecolor',[1 1 1]); xlabel('Q rob ascending order'); ylabel('S'); 
subplot(323), plot(numsubunitspikes(ind_rob),'ro','MarkerFacecolor',[1 0 0],'MarkerEdgecolor',[1 1 1]); xlabel('Q rob ascending order'); ylabel('num spikes'); title(strcat(num2str(r2,2),',',' ',num2str(p2,2)));
subplot(324), plot(numsubunitspikes(ind_sse),'bo','MarkerFacecolor',[0 0 1],'MarkerEdgecolor',[1 1 1]); xlabel('Q rob ascending order'); ylabel('num spikes'); 
subplot(325), plot(RF_loc(ind_rob),'ro','MarkerFacecolor',[1 0 0],'MarkerEdgecolor',[1 1 1]); xlabel('Q rob ascending order'); ylabel('ecc'); title(strcat(num2str(r3,2),',',' ',num2str(p3,2)));
subplot(326), plot(RF_loc(ind_sse),'bo','MarkerFacecolor',[0 0 1],'MarkerEdgecolor',[1 1 1]); xlabel('Q rob ascending order'); ylabel('ecc'); 
plot_counter = plot_counter + 1;

% Trying to figure if OFF responses are correlated with the linear or nonlinear isoreponse contour
[r1,p1] = corr(pranksumOFF,Q_rob,'Type','Spearman');
figure(plot_counter); subplot(121);plot(log10(pranksumOFF(ind_rob)),'ro', 'MarkerFacecolor',[1 0 0],'MarkerEdgecolor',[1 1 1]); title(strcat(num2str(r1,2),',',' ',num2str(p1,2))); ('Q rob ascending order'); ylabel('p ranksum');
subplot(122);plot(log10(pranksumOFF(ind_sse)),'bo','MarkerFacecolor',[0 0 1],'MarkerEdgecolor',[1 1 1]); xlabel('Q sse ascending order'); ylabel('p ranksum OFF');
plot_counter = plot_counter + 1;

% Comparing the Q_rob values for DO, LUM and hardtoclassify cells
[p1,h1] = ranksum(Q_rob([OCidx]),Q_rob([LMidx]));
[p2,h2] = ranksum(Q_rob([OCidx LMidx]),Q_rob([LUMidx]));
[p3,h3] = ranksum(Q_rob([OCidx LMidx]),Q_rob([hardtoclassifyidx]));
[p4,h4] = ranksum(Q_rob([LUMidx]),Q_rob([hardtoclassifyidx]));
bins = 0:0.5:10;
figure(plot_counter); subplot(221); histogram(Q_rob,bins); hold on; plot(median(Q_rob),0,'kv','Markerfacecolor',[1 0 0]); xlabel('Q_rob'); title('All cells'); 
subplot(222); histogram(Q_rob([OCidx]),bins); hold on; plot(median(Q_rob([OCidx])),0,'kv','Markerfacecolor',[0 0 1]); histogram(Q_rob([LMidx]),bins); plot(median(Q_rob([LMidx])),0,'kv','Markerfacecolor',[1 0 0]); xlabel('Q_rob'); title('OC vs LM'); 
subplot(223); histogram(Q_rob([OCidx LMidx]),bins); hold on; plot(median(Q_rob([OCidx LMidx])),0,'kv','Markerfacecolor',[0 0 1]); histogram(Q_rob([LUMidx]),bins); plot(median(Q_rob([LUMidx])),0,'kv','Markerfacecolor',[1 0 0]); xlabel('Q_rob'); title('DO vs LUM'); 
subplot(224); histogram(Q_rob([OCidx LMidx]),bins); hold on; plot(median(Q_rob([OCidx LMidx])),0,'kv','Markerfacecolor',[0 0 1]); histogram(Q_rob([hardtoclassifyidx]),bins); plot(median(Q_rob([hardtoclassifyidx])),0,'kv','Markerfacecolor',[1 0 0]); xlabel('Q_rob'); title('DO vs htc'); 
plot_counter = plot_counter + 1;

[r1,p1] = corr(RSSE_linearmodel-RSSE_quadmodel,Q_rob,'Type','Spearman');
figure(plot_counter);
subplot(221); plot(log10(RSSE_linearmodel),log10(RSSE_quadmodel),'o','MarkerFacecolor',[1 0 0],'MarkerEdgecolor',[1 1 1]); xlabel('RSSE lin'); ylabel('RSSE quad');
subplot(222); plot(log10(SSE_linearmodel),log10(SSE_quadmodel),'o','MarkerFacecolor',[0 0 1],'MarkerEdgecolor',[1 1 1]); xlabel('SSE lin'); ylabel('SSE quad');
subplot(223); histogram(RSSE_linearmodel-RSSE_quadmodel,50); xlabel('RSSE lin - quad'); title('Robust regression');
subplot(224); plot(RSSE_linearmodel-RSSE_quadmodel,Q_rob,'o','MarkerFacecolor',[1 0 0],'MarkerEdgecolor',[1 1 1]); xlabel('RSSE lin - quad'); ylabel('Q_rob');
plot_counter = plot_counter + 1;

%% Visualizing the linearity/non-linearity of isoresponse contour in an ascending order (Q_rob) 
% Displaying the 16 most linear cell and most nonlinear cell last
% and some other variables like absS, deltaBIC etc for those cells 

numsubplots = 4;
mostlinear = 1:numsubplots^2;
mostnonlinear = numel(idx):-1:numel(idx)-numsubplots^2+1; 
N = [mostlinear mostnonlinear];
count = 1;
prod_eigvalues = zeros(numel(N),1);
for ii = 1:numel(N)
    ind = ind_rob(N(ii));
    THETA = THETA_all{1,ind};
    THETA = THETA * pi/180; % converting to radians
    if any(THETA>(135*pi/180))
        allthetas = linspace(-pi,pi,100);
        newtheta = linspace(-pi,pi,101);
    else
        allthetas = linspace(-pi/4,3*pi/4,100);
        newtheta = linspace(-pi/4,3*pi/4,101);
    end
    RHO = RHO_all{1,ind};
    oog_idx = oog_idx_all{1,ind};
    not_oog_idx = not_oog_idx_all{1,ind};
    outofgamut = zeros(size(THETA));
    outofgamut(oog_idx) = 1;
    outofgamut = logical(outofgamut);
    [x_orig, y_orig] = pol2cart(THETA,RHO);            
             
    %Plotting the figures
    figure(plot_counter),subplot(numsubplots,numsubplots,count); plot(x_orig(not_oog_idx), y_orig(not_oog_idx),'o','MarkerSize',4,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[1 1 1]); hold on;
    if any(outofgamut)
        hold on; plot(upsample(x_orig(outofgamut),2),upsample(y_orig(outofgamut),2),'k');
    end
    set(gca,'XLim',[-1.5 1.5],'YLim',[-1.5 1.5],'XTick',[],'YTick',[]); title(num2str(Q_rob(ind),3)); drawnow; axis square; hold off;
    
    % calculating product of eigenvalues
    [~,d] = eig([quad_modelparams(ind(1),1) quad_modelparams(ind(1),3)/2; quad_modelparams(ind(1),3)/2 quad_modelparams(ind(1),2)]);
    prod_eigvalues(ii) = prod(diag(d));
%     prod_eigenvalues(ii) = quad_modelparams(ind,3)^2 - 4*quad_modelparams(ind,1)*quad_modelparams(ind,2); 
    
    count = count + 1;
    if count > numsubplots^2
        plot_counter = plot_counter + 1;
        count  = 1;
    end
end
plot_counter = plot_counter + 1;
% Using Wilconxon-rank sum test to test if there are any systematic differences between the least and the most linear cells
x = [ones(numel(mostlinear),1); 2*ones(numel(mostnonlinear),1)];
g1 = deltaBIC(ind_rob([mostlinear mostnonlinear])); [p1,~] = ranksum(deltaBIC(ind_rob(mostlinear)),deltaBIC(ind_rob(mostnonlinear)));
g2 = Q_rob(ind_rob([mostlinear mostnonlinear]));[p2,~] = ranksum(Q_rob(ind_rob(mostlinear)),Q_rob(ind_rob(mostnonlinear)));
g3 = absS(ind_rob([mostlinear mostnonlinear])); [p3,~] = ranksum(absS(ind_rob(mostlinear)),absS(ind_rob(mostnonlinear)));
g4 = z_scores(ind_rob([mostlinear mostnonlinear])); [p4,~] = ranksum(z_scores(ind_rob(mostlinear)),z_scores(ind_rob(mostnonlinear)));
g5 = anglebwvectorsRGB(ind_rob([mostlinear mostnonlinear])); [p5,~] = ranksum(anglebwvectorsRGB(ind_rob(mostlinear)),anglebwvectorsRGB(ind_rob(mostnonlinear)));
g6 = anglediffWNchecksubunit(ind_rob([mostlinear mostnonlinear])); [p6,~] = ranksum(anglediffWNchecksubunit(ind_rob(mostlinear)),anglediffWNchecksubunit(ind_rob(mostnonlinear)));
g7 = NLindex_r(ind_rob([mostlinear mostnonlinear])); [p7,~] = ranksum(NLindex_r(ind_rob(mostlinear)),NLindex_r(ind_rob(mostnonlinear)));
g8 = ratioeig(ind_rob([mostlinear mostnonlinear])); [p8,~] = ranksum(ratioeig(ind_rob(mostlinear)),ratioeig(ind_rob(mostnonlinear)));
g9 = absLminusM(ind_rob([mostlinear mostnonlinear])); [p9,~] = ranksum(absLminusM(ind_rob(mostlinear)),absLminusM(ind_rob(mostnonlinear)));
g10 = absLplusM(ind_rob([mostlinear mostnonlinear])); [p10,~] = ranksum(absLplusM(ind_rob(mostlinear)),absLplusM(ind_rob(mostnonlinear)));
figure(plot_counter); set(gcf,'Name',strcat('Boxplots for ',num2str(numsubplots^2),'lin and non-lin cells: Wilcoxon rank sum test'));
subplot(251); boxplot(g1,x,'Notch','on','Labels',{'lin','nonlin'}); title(strcat('deltaBIC',{' '},num2str(p1)));
subplot(252); boxplot(g2,x,'Notch','on','Labels',{'lin','nonlin'}); title(strcat('Q rob: sanity check',{' '},num2str(p2)));
subplot(253); boxplot(g3,x,'Notch','on','Labels',{'lin','nonlin'}); title(strcat('abs S',{' '},num2str(p3)));
subplot(254); boxplot(g4,x,'Notch','on','Labels',{'lin','nonlin'}); title(strcat('z score',{' '},num2str(p4)));
subplot(255); boxplot(g5,x,'Notch','on','Labels',{'lin','nonlin'}); title(strcat('angle bt sub',{' '},num2str(p5)));
subplot(256); boxplot(g6,x,'Notch','on','Labels',{'lin','nonlin'}); title(strcat('angle bt check sub',{' '},num2str(p6)));
subplot(257); boxplot(g7,x,'Notch','on','Labels',{'lin','nonlin'}); title(strcat('NL index r',{' '},num2str(p7)));
subplot(258); boxplot(g8,x,'Notch','on','Labels',{'lin','nonlin'}); title(strcat('% PC1: ratioeig',{' '},num2str(p8)));
subplot(259); boxplot(g9,x,'Notch','on','Labels',{'lin','nonlin'}); title(strcat('abs L minus M',{' '},num2str(p9)));
subplot(2,5,10); boxplot(g10,x,'Notch','on','Labels',{'lin','nonlin'}); title(strcat('abs L plus M',{' '},num2str(p10)));
plot_counter = plot_counter + 1;

% Want to check how many of the nonlinear cells are narrower-than-linear vs broader-than-linear
broaderthanlinear = prod_eigvalues(numel(N)/2+1:numel(N))>0;
narrowerthanlinear = prod_eigvalues(numel(N)/2+1:numel(N))<0;
figure(plot_counter); set(gcf,'Name','Nature of nonlinearity');
bar([sum(broaderthanlinear) sum(narrowerthanlinear)]); set(gca,'XTickLabel',{'broader','narrower'}); ylabel('No. of cells');
plot_counter = plot_counter + 1;

% Need to look at the signal integration between and across subunits from WN analysis and FR surface fits
pthresh = 0.01;
f1 = diffmodelpred(ind_rob([mostlinear mostnonlinear])); [p11,~] = ranksum(diffmodelpred(ind_rob(mostlinear)),diffmodelpred(ind_rob(mostnonlinear))); % WN - across subunit analysis
f2 = diffmodelpredFR(ind_rob([mostlinear mostnonlinear])); [p12,~] = ranksum(diffmodelpredFR(ind_rob(mostlinear)),diffmodelpredFR(ind_rob(mostnonlinear))); % WN - across subunit analysis
f3 = [diffmodelpred1(ind_rob([mostlinear mostnonlinear])); diffmodelpred2(ind_rob([mostlinear mostnonlinear]))]; % WN - within subunit
[p13,~] = ranksum([diffmodelpred1(ind_rob(mostlinear));diffmodelpred2(ind_rob(mostlinear))],[diffmodelpred1(ind_rob(mostnonlinear)); diffmodelpred2(ind_rob(mostnonlinear))]);
f4 = [sigquadterms1(ind_rob([mostlinear mostnonlinear])); sigquadterms2(ind_rob([mostlinear mostnonlinear]))]; % significant subunits within subunits 
[p14,~] = ranksum([sigquadterms1(ind_rob(mostlinear));sigquadterms2(ind_rob(mostlinear))],[sigquadterms1(ind_rob(mostnonlinear));sigquadterms2(ind_rob(mostnonlinear))]);
f5 = sigquadterms(ind_rob([mostlinear mostnonlinear])); [p15,~] = ranksum(sigquadterms(ind_rob(mostlinear)),sigquadterms(ind_rob(mostnonlinear))); % significant subunits within subunits 
figure(plot_counter); set(gcf,'Name',strcat('Boxplots for ',num2str(numsubplots^2),'lin and non-lin cells'));
subplot(231); boxplot(f1,x,'Notch','on','Labels',{'lin','nonlin'}); title(strcat('Across subunits',{' '},num2str(p11)));
subplot(232); boxplot(f2,x,'Notch','on','Labels',{'lin','nonlin'}); title(strcat('FR surface fits',{' '},num2str(p12)));
subplot(233); boxplot(f3,[x;x],'Notch','on','Labels',{'lin','nonlin'}); title(strcat('Within subunits',{' '},num2str(p13)));
subplot(234); boxplot(f4,[x;x],'Notch','on','Labels',{'lin','nonlin'}); title(strcat('sig quad terms within subunits',{' '},num2str(p14)));
subplot(235); boxplot(f5,x,'Notch','on','Labels',{'lin','nonlin'}); title(strcat('sig quad terms across subunits',{' '},num2str(p15)));
subplot(236); bar([sum([Fprob1(ind_rob(mostlinear)); Fprob2(ind_rob(mostlinear))]<pthresh) sum([Fprob1(ind_rob(mostnonlinear)); Fprob2(ind_rob(mostnonlinear))]<pthresh)]);
set(gca,'xticklabel',{'lin','nonlin'}); ylabel('nonlin subunits'); title(strcat('Fprob<',num2str(pthresh)));
plot_counter = plot_counter + 1;

% Trying to see if a subunit with a higher S cone weight is likely to have non-linear signal integration or not
% I am looking at the absolute S cone weight in individual subunit
newabsS = abs([S1LMS(3,:)'; S2LMS(3,:)']);
newdiffmodelpred = [diffmodelpred1; diffmodelpred2];
newsigquadterms = [sigquadterms1; sigquadterms2];
Fprobsubunit = [Fprob1; Fprob2];
[~,indnewabsS] = sort(newabsS);
numsubunitstopick = 16;
pthresh = 0.01;
leastSconewt = 1:numsubunitstopick;
mostSconewt = numel(newabsS):-1:numel(newabsS)-numsubunitstopick+1;
x = [ones(numel(leastSconewt),1); 2*ones(numel(mostSconewt),1)];
g1 = newdiffmodelpred(indnewabsS([leastSconewt mostSconewt])); p1 =  ranksum(newdiffmodelpred(indnewabsS(leastSconewt)),newdiffmodelpred(indnewabsS(mostSconewt)));
g2 = newsigquadterms(indnewabsS([leastSconewt mostSconewt]));  p2 =  ranksum(newsigquadterms(indnewabsS(leastSconewt)),newsigquadterms(indnewabsS(mostSconewt)));
g3 = Fprobsubunit(indnewabsS([leastSconewt mostSconewt]));  p3 =  ranksum(Fprobsubunit(indnewabsS(leastSconewt)),newsigquadterms(indnewabsS(mostSconewt)));
figure(plot_counter); set(gcf,'Name',strcat('Boxplots for ',num2str(numsubunitstopick),'least and most S cone wt subunits')); 
subplot(231); boxplot(g1,x,'Notch','on','Labels',{'least','most'}); title(strcat('Diff in pred',{' '},num2str(p1)));
subplot(232); boxplot(g2,x,'Notch','on','Labels',{'least','most'}); title(strcat('sig quad terms',{' '},num2str(p2)));
subplot(233); boxplot(g3,x,'Notch','on','Labels',{'least','most'}); title(strcat('Fprob',{' '},num2str(p3)));
subplot(234); histogram(newdiffmodelpred(indnewabsS(leastSconewt))); hold on; histogram(newdiffmodelpred(indnewabsS(mostSconewt))); title('Diff in pred'); hold off; 
subplot(235); histogram(newsigquadterms(indnewabsS(leastSconewt))); hold on; histogram(newsigquadterms(indnewabsS(mostSconewt))); title('sig quad terms'); hold off;
subplot(236); bar([sum(Fprobsubunit(indnewabsS(leastSconewt))<pthresh); sum(Fprobsubunit(indnewabsS(mostSconewt))<pthresh)]); set(gca,'xticklabel',{'least','most'}); ylabel('nonlin subunits');  hold off;
plot_counter = plot_counter + 1;

% Trying to see if the "cumulative" absolute S cone weight can explain some the nonlinearities within and across subunit. 
% This is different from the previous analyses where I was sorting subunits based on S cone weights and not the cell 
[~,indabsS] = sort(absS);
numcellstopick = 16;
leastSconesig = 1:numcellstopick;
mostSconesig = numel(absS):-1:numel(absS)-numcellstopick+1;
x = [ones(numel(leastSconesig),1); 2*ones(numel(mostSconesig),1)];
g1 = Q_rob(indabsS([leastSconesig mostSconesig])); p1 = ranksum(Q_rob(indabsS(leastSconesig)),Q_rob(indabsS(mostSconesig)));
g2 = diffmodelpred(indabsS([leastSconesig mostSconesig])); p2 = ranksum(diffmodelpred(indabsS(leastSconesig)),diffmodelpred(indabsS(mostSconesig)));
g3 = diffmodelpredFR(indabsS([leastSconesig mostSconesig])); p3 = ranksum(diffmodelpredFR(indabsS(leastSconesig)),diffmodelpredFR(indabsS(mostSconesig)));
figure(plot_counter); set(gcf,'Name',strcat('Boxplots for ',num2str(numcellstopick),'least and most S cone signal')); 
subplot(221); boxplot(g1,x,'Notch','on','Labels',{'least','most'}); title(strcat('Q rob',{' '},num2str(p1)));
subplot(222); boxplot(g2,x,'Notch','on','Labels',{'least','most'}); title(strcat('Diff in pred',{' '},num2str(p2)));
subplot(223); boxplot(g3,x,'Notch','on','Labels',{'least','most'}); title(strcat('Diff in pred FR',{' '},num2str(p3)));
plot_counter = plot_counter + 1;



%% Trying some data visualization techniques - tSNE, dendrograms and 
Y_euc = tsne([S1LMS' S2LMS'],'algorithm','exact','Distance','euclidean'); % tSNE - euclidean
Y_cos = tsne([S1LMS' S2LMS'],'algorithm','exact','Distance','cosine'); % tSNE - cosine
Y_mah = tsne([S1LMS' S2LMS'],'algorithm','exact','Distance','mahalanobis'); % tSNE - mahalabobis
Y_che = tsne([S1LMS' S2LMS'],'algorithm','exact','Distance','chebychev'); % tSNE - chebychev
tree = linkage([S1LMS' S2LMS'],'average','euclidean');
Ymds = mdscale(pdist([S1LMS' S2LMS']),2);
colorstepsize = 1/numel(idx);
figure(plot_counter); set(gcf,'Name', 'Data visualization');
for ii = 1:numel(idx)
    subplot(231);plot(Y_euc(ind_rob(ii),1),Y_euc(ind_rob(ii),2),'o','MarkerFaceColor',[ii*colorstepsize 0 1-ii*colorstepsize],'MarkerEdgeColor',[1 1 1]); hold on; 
    subplot(232);plot(Y_cos(ind_rob(ii),1),Y_cos(ind_rob(ii),2),'o','MarkerFaceColor',[ii*colorstepsize 0 1-ii*colorstepsize],'MarkerEdgeColor',[1 1 1]); hold on;
    subplot(233);plot(Y_mah(ind_rob(ii),1),Y_mah(ind_rob(ii),2),'o','MarkerFaceColor',[ii*colorstepsize 0 1-ii*colorstepsize],'MarkerEdgeColor',[1 1 1]); hold on; 
    subplot(234);plot(Y_che(ind_rob(ii),1),Y_che(ind_rob(ii),2),'o','MarkerFaceColor',[ii*colorstepsize 0 1-ii*colorstepsize],'MarkerEdgeColor',[1 1 1]); hold on;
    subplot(235);plot(Ymds(ind_rob(ii),1),Ymds(ind_rob(ii),2),'o','MarkerFaceColor',[ii*colorstepsize 0 1-ii*colorstepsize],'MarkerEdgeColor',[1 1 1]); hold on; 
end
subplot(231); title('tSNE- euclidean'); hold off;
subplot(232); title('tSNE- cosine'); hold off;
subplot(233); title('tSNE- mahalanobis'); hold off;
subplot(234); title('tSNE- chebychev'); hold off;
subplot(235); title('MDS'); hold off;
subplot(236);dendrogram(tree,0); title('Dendrogram');

plot_counter = plot_counter + 1;