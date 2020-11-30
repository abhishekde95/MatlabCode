% Plotting the results from 3 subjects, Propixx Data
% Author- Abhishek De, 08/17

close all; clearvars;
% Foveal Data
load('N:\Abhishek_Psychophysics_data\ProPixx\fitting_data\Abhishek\foveal_fit\fit_data.mat');abhi_fitdata_f = fit_data;
load('N:\Abhishek_Psychophysics_data\ProPixx\fitting_data\Emily\foveal_fit\fit_data.mat'); emily_fitdata_f = fit_data;
load('N:\Abhishek_Psychophysics_data\ProPixx\fitting_data\Greg\foveal_fit\fit_data.mat'); greg_fitdata_f = fit_data;
% Parafoveal Data
load('N:\Abhishek_Psychophysics_data\ProPixx\fitting_data\Abhishek\5_deg_fit\fit_data.mat');abhi_fitdata_p = fit_data;
load('N:\Abhishek_Psychophysics_data\ProPixx\fitting_data\Emily\5_deg_fit\fit_data.mat');emily_fitdata_p = fit_data;
load('N:\Abhishek_Psychophysics_data\ProPixx\fitting_data\Greg\5_deg_fit\fit_data.mat');greg_fitdata_p = fit_data;
clear fit_data
% Obtaining the model Errs
load('N:\Abhishek_Psychophysics_data\ProPixx\fitting_data\Abhishek\foveal_fit\modelErrs.mat'); modelErrs_abhi_f = modelErrs;
load('N:\Abhishek_Psychophysics_data\ProPixx\fitting_data\Abhishek\5_deg_fit\modelErrs.mat'); modelErrs_abhi_p = modelErrs;
load('N:\Abhishek_Psychophysics_data\ProPixx\fitting_data\Emily\foveal_fit\modelErrs.mat'); modelErrs_emily_f = modelErrs;
load('N:\Abhishek_Psychophysics_data\ProPixx\fitting_data\Emily\5_deg_fit\modelErrs.mat'); modelErrs_emily_p = modelErrs;
load('N:\Abhishek_Psychophysics_data\ProPixx\fitting_data\Greg\foveal_fit\modelErrs.mat'); modelErrs_greg_f = modelErrs;
load('N:\Abhishek_Psychophysics_data\ProPixx\fitting_data\Greg\5_deg_fit\modelErrs.mat'); modelErrs_greg_p = modelErrs;
clear modelErrs

abhi_lsf_oc_f = abhi_fitdata_f.model3(1); abhi_lsf_lm_f = abhi_fitdata_f.model3(3);
abhi_hsf_oc_f = abhi_fitdata_f.model4(1); abhi_hsf_lm_f = abhi_fitdata_f.model4(3);
abhi_lsf_oc_p = abhi_fitdata_p.model3(1); abhi_lsf_lm_p = abhi_fitdata_p.model3(3);
abhi_hsf_oc_p = abhi_fitdata_p.model4(1); abhi_hsf_lm_p = abhi_fitdata_p.model4(3);

emily_lsf_oc_f = emily_fitdata_f.model3(1); emily_lsf_lm_f = emily_fitdata_f.model3(3);
emily_hsf_oc_f = emily_fitdata_f.model4(1); emily_hsf_lm_f = emily_fitdata_f.model4(3);
emily_lsf_oc_p = emily_fitdata_p.model3(1); emily_lsf_lm_p = emily_fitdata_p.model3(3);
emily_hsf_oc_p = emily_fitdata_p.model4(1); emily_hsf_lm_p = emily_fitdata_p.model4(3);

greg_lsf_oc_f = greg_fitdata_f.model3(1); greg_lsf_lm_f = greg_fitdata_f.model3(3);
greg_hsf_oc_f = greg_fitdata_f.model4(1); greg_hsf_lm_f = greg_fitdata_f.model4(3);
greg_lsf_oc_p = greg_fitdata_p.model3(1); greg_lsf_lm_p = greg_fitdata_p.model3(3);
greg_hsf_oc_p = greg_fitdata_p.model4(1); greg_hsf_lm_p = greg_fitdata_p.model4(3);

% Plotting the results
% Each figure will have 3 subplots corresponding to 3 subjects
% x - axis is eccentricity: fovea and parafovea
figure(1),subplot(131),errorbar([1; 2],[abhi_lsf_oc_f; abhi_lsf_oc_p],[modelErrs_abhi_f.low_sf; modelErrs_abhi_p.low_sf],'b','Linewidth',1); hold on;
errorbar([1; 2],[abhi_hsf_oc_f; abhi_hsf_oc_p],[modelErrs_abhi_f.high_sf; modelErrs_abhi_p.high_sf],'b','Linewidth',3);hold on;
errorbar([1; 2],[abhi_lsf_lm_f; abhi_lsf_lm_p],[modelErrs_abhi_f.low_sf; modelErrs_abhi_p.low_sf],'m','Linewidth',1);hold on;
errorbar([1; 2],[abhi_hsf_lm_f; abhi_hsf_lm_p],[modelErrs_abhi_f.high_sf; modelErrs_abhi_p.high_sf],'m','Linewidth',3);hold on;
legend('OC-L','OC-H','LM-L','LM-H'); set(gca,'Xlim',[0.5 2.5],'Yscale','log'); title('AD'); ylabel('Contrast');hold off;

subplot(132),errorbar([1; 2],[emily_lsf_oc_f; emily_lsf_oc_p],[modelErrs_emily_f.low_sf; modelErrs_emily_p.low_sf],'b','Linewidth',1); hold on;
errorbar([1; 2],[emily_hsf_oc_f; emily_hsf_oc_p],[modelErrs_emily_f.high_sf; modelErrs_emily_p.high_sf],'b','Linewidth',3);hold on;
errorbar([1; 2],[emily_lsf_lm_f; emily_lsf_lm_p],[modelErrs_emily_f.low_sf; modelErrs_emily_p.low_sf],'m','Linewidth',1);hold on;
errorbar([1; 2],[emily_hsf_lm_f; emily_hsf_lm_p],[modelErrs_emily_f.high_sf; modelErrs_emily_p.high_sf],'m','Linewidth',3);hold on;
legend('OC-L','OC-H','LM-L','LM-H'); set(gca,'Xlim',[0.5 2.5],'Yscale','log'); title('EG'); ylabel('Contrast');hold off;

subplot(133),errorbar([1; 2],[greg_lsf_oc_f; greg_lsf_oc_p],[modelErrs_greg_f.low_sf; modelErrs_greg_p.low_sf],'b','Linewidth',1); hold on;
errorbar([1; 2],[greg_hsf_oc_f; greg_hsf_oc_p],[modelErrs_greg_f.high_sf; modelErrs_greg_p.high_sf],'b','Linewidth',3);hold on;
errorbar([1; 2],[greg_lsf_lm_f; greg_lsf_lm_p],[modelErrs_greg_f.low_sf; modelErrs_greg_p.low_sf],'m','Linewidth',1);hold on;
errorbar([1; 2],[greg_hsf_lm_f; greg_hsf_lm_p],[modelErrs_greg_f.high_sf; modelErrs_greg_p.high_sf],'m','Linewidth',3);hold on;
legend('OC-L','OC-H','LM-L','LM-H'); set(gca,'Xlim',[0.5 2.5],'Yscale','log'); title('GH'); ylabel('Contrast');hold off;

% x - axis is spatial frequency: low and high sf
figure(2),subplot(131),errorbar([1; 2],[abhi_lsf_oc_f; abhi_hsf_oc_f],[modelErrs_abhi_f.low_sf; modelErrs_abhi_f.high_sf],'b','Linewidth',1); hold on;
errorbar([1; 2],[abhi_lsf_oc_p; abhi_hsf_oc_p],[modelErrs_abhi_p.low_sf; modelErrs_abhi_p.high_sf],'b','Linewidth',3); hold on;
errorbar([1; 2],[abhi_lsf_lm_f; abhi_hsf_lm_f],[modelErrs_abhi_f.low_sf; modelErrs_abhi_f.high_sf],'m','Linewidth',1); hold on;
errorbar([1; 2],[abhi_lsf_lm_p; abhi_hsf_lm_p],[modelErrs_abhi_p.low_sf; modelErrs_abhi_p.high_sf],'m','Linewidth',3); hold on;
legend('OC-F','OC-P','LM-F','LM-P'); set(gca,'Xlim',[0.5 2.5],'Yscale','log'); title('AD'); ylabel('Contrast');hold off;

subplot(132),errorbar([1; 2],[emily_lsf_oc_f; emily_hsf_oc_f],[modelErrs_emily_f.low_sf; modelErrs_emily_f.high_sf],'b','Linewidth',1); hold on;
errorbar([1; 2],[emily_lsf_oc_p; emily_hsf_oc_p],[modelErrs_emily_p.low_sf; modelErrs_emily_p.high_sf],'b','Linewidth',3); hold on;
errorbar([1; 2],[emily_lsf_lm_f; emily_hsf_lm_f],[modelErrs_emily_f.low_sf; modelErrs_emily_f.high_sf],'m','Linewidth',1); hold on;
errorbar([1; 2],[emily_lsf_lm_p; emily_hsf_lm_p],[modelErrs_emily_p.low_sf; modelErrs_emily_p.high_sf],'m','Linewidth',3); hold on;
legend('OC-F','OC-P','LM-F','LM-P'); set(gca,'Xlim',[0.5 2.5],'Yscale','log'); title('EG'); ylabel('Contrast');hold off;

subplot(133),errorbar([1; 2],[greg_lsf_oc_f; greg_hsf_oc_f],[modelErrs_greg_f.low_sf; modelErrs_greg_f.high_sf],'b','Linewidth',1); hold on;
errorbar([1; 2],[greg_lsf_oc_p; greg_hsf_oc_p],[modelErrs_greg_p.low_sf; modelErrs_greg_p.high_sf],'b','Linewidth',3); hold on;
errorbar([1; 2],[greg_lsf_lm_f; greg_hsf_lm_f],[modelErrs_greg_f.low_sf; modelErrs_greg_f.high_sf],'m','Linewidth',1); hold on;
errorbar([1; 2],[greg_lsf_lm_p; greg_hsf_lm_p],[modelErrs_greg_p.low_sf; modelErrs_greg_p.high_sf],'m','Linewidth',3); hold on;
legend('OC-F','OC-P','LM-F','LM-P'); set(gca,'Xlim',[0.5 2.5],'Yscale','log'); title('GH'); ylabel('Contrast');hold off;

%% Some additional figures from subject AD for comparing the results of Rig3 and Propixx data
load('C:\Users\Abhishek\Desktop\MatlabCode\Abhishek\for_FVM\foveal\threshold_abhi_f.mat'); abhi_crtfitdata_f = threshold_abhi_f;
load('C:\Users\Abhishek\Desktop\MatlabCode\Abhishek\for_FVM\5_degree\threshold_abhi_p.mat'); abhi_crtfitdata_p = threshold_abhi_p;
load('C:\Users\Abhishek\Desktop\MatlabCode\Abhishek\for_FVM\foveal\modelErrs_abhi_f.mat'); modelErrs_crtabhi_f = modelErrs;
load('C:\Users\Abhishek\Desktop\MatlabCode\Abhishek\for_FVM\5_degree\modelErrs_abhi_p.mat'); modelErrs_crtabhi_p = modelErrs;
clear modelErrs fit_data

lowsf_abhi_par = [1 threshold_abhi_p.lowsf.son; 2 threshold_abhi_p.lowsf.soff];
highsf_abhi_par = [1 threshold_abhi_p.highsf.son; 2 threshold_abhi_p.highsf.soff];
lowsf_abhi_fov = [1 threshold_abhi_f.lowsf.son; 2 threshold_abhi_f.lowsf.soff];
highsf_abhi_fov = [1 threshold_abhi_f.highsf.son; 2 threshold_abhi_f.highsf.soff];

% x - axis is eccentricity in these figures
figure(3), subplot(121),errorbar([1; 2],[abhi_lsf_oc_f; abhi_lsf_oc_p],[modelErrs_abhi_f.low_sf; modelErrs_abhi_p.low_sf],'b','Linewidth',1); hold on;
errorbar([1; 2],[abhi_hsf_oc_f; abhi_hsf_oc_p],[modelErrs_abhi_f.high_sf; modelErrs_abhi_p.high_sf],'b','Linewidth',3);hold on;
errorbar([1; 2],[abhi_lsf_lm_f; abhi_lsf_lm_p],[modelErrs_abhi_f.low_sf; modelErrs_abhi_p.low_sf],'m','Linewidth',1);hold on;
errorbar([1; 2],[abhi_hsf_lm_f; abhi_hsf_lm_p],[modelErrs_abhi_f.high_sf; modelErrs_abhi_p.high_sf],'m','Linewidth',3);hold on;
legend('OC-L','OC-H','LM-L','LM-H'); set(gca,'Xlim',[0.5 2.5],'Yscale','log'); title('AD-Propixx'); ylabel('Contrast');hold off;

subplot(122), errorbar([1; 2],[threshold_abhi_f.lowsf.soff; threshold_abhi_p.lowsf.soff],[modelErrs_crtabhi_f.low_sf(1); modelErrs_crtabhi_p.low_sf(1)],'b','Linewidth',1); hold on;
errorbar([1; 2],[threshold_abhi_f.highsf.soff; threshold_abhi_p.highsf.soff],[modelErrs_crtabhi_f.high_sf(1); modelErrs_crtabhi_p.high_sf(1)],'b','Linewidth',3); hold on;
errorbar([1; 2],[threshold_abhi_f.lowsf.son; threshold_abhi_p.lowsf.son],[modelErrs_crtabhi_f.low_sf(3); modelErrs_crtabhi_p.low_sf(3)],'m','Linewidth',1); hold on;
errorbar([1; 2],[threshold_abhi_f.highsf.son; threshold_abhi_p.highsf.son],[modelErrs_crtabhi_f.high_sf(3); modelErrs_crtabhi_p.high_sf(3)],'m','Linewidth',3); hold on;
legend('OC-L','OC-H','LM-L','LM-H'); set(gca,'Xlim',[0.5 2.5],'Yscale','log'); title('AD-CRT'); ylabel('Contrast');hold off;

% x - axis is spatial frequency in these figures
figure(4),subplot(121),errorbar([1; 2],[abhi_lsf_oc_f; abhi_hsf_oc_f],[modelErrs_abhi_f.low_sf; modelErrs_abhi_f.high_sf],'b','Linewidth',1); hold on;
errorbar([1; 2],[abhi_lsf_oc_p; abhi_hsf_oc_p],[modelErrs_abhi_p.low_sf; modelErrs_abhi_p.high_sf],'b','Linewidth',3); hold on;
errorbar([1; 2],[abhi_lsf_lm_f; abhi_hsf_lm_f],[modelErrs_abhi_f.low_sf; modelErrs_abhi_f.high_sf],'m','Linewidth',1); hold on;
errorbar([1; 2],[abhi_lsf_lm_p; abhi_hsf_lm_p],[modelErrs_abhi_p.low_sf; modelErrs_abhi_p.high_sf],'m','Linewidth',3); hold on;
legend('OC-F','OC-P','LM-F','LM-P'); set(gca,'Xlim',[0.5 2.5],'Yscale','log'); title('AD-Propixx'); ylabel('Contrast');hold off;

subplot(122),errorbar([1; 2],[threshold_abhi_f.lowsf.soff; threshold_abhi_f.highsf.soff],[modelErrs_crtabhi_f.low_sf(1); modelErrs_crtabhi_f.high_sf(1)],'b','Linewidth',1); hold on;
errorbar([1; 2],[threshold_abhi_p.lowsf.soff; threshold_abhi_p.highsf.soff],[modelErrs_crtabhi_p.low_sf(1); modelErrs_crtabhi_p.high_sf(1)],'b','Linewidth',3); hold on;
errorbar([1; 2],[threshold_abhi_f.lowsf.son; threshold_abhi_f.highsf.son],[modelErrs_crtabhi_f.low_sf(3); modelErrs_crtabhi_f.high_sf(3)],'m','Linewidth',1); hold on;
errorbar([1; 2],[threshold_abhi_p.lowsf.son; threshold_abhi_p.highsf.son],[modelErrs_crtabhi_p.low_sf(3); modelErrs_crtabhi_p.high_sf(3)],'m','Linewidth',3); hold on;
legend('OC-F','OC-P','LM-F','LM-P'); set(gca,'Xlim',[0.5 2.5],'Yscale','log'); title('AD-CRT'); ylabel('Contrast');hold off;
