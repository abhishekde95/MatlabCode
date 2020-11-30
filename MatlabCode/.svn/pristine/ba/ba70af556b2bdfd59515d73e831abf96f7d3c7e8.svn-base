% Population_analysis_for_FVM_2016
% Author - Abhishek De

% Foveal data
load('C:\Users\Abhishek\Desktop\MatlabCode\Abhishek\for_FVM\foveal\threshold_abhi_f.mat')
load('C:\Users\Abhishek\Desktop\MatlabCode\Abhishek\for_FVM\foveal\threshold_emily_f.mat')
load('C:\Users\Abhishek\Desktop\MatlabCode\Abhishek\for_FVM\foveal\threshold_greg_f.mat')
% Parafoveal data
load('C:\Users\Abhishek\Desktop\MatlabCode\Abhishek\for_FVM\5_degree\threshold_abhi_p.mat')
load('C:\Users\Abhishek\Desktop\MatlabCode\Abhishek\for_FVM\5_degree\threshold_emily_p.mat')
load('C:\Users\Abhishek\Desktop\MatlabCode\Abhishek\for_FVM\5_degree\threshold_greg_p.mat')

% Standard errors from the fit
load('C:\Users\Abhishek\Desktop\MatlabCode\Abhishek\for_FVM\foveal\modelErrs_abhi_f.mat'); modelErrs_abhi_f = modelErrs;
load('C:\Users\Abhishek\Desktop\MatlabCode\Abhishek\for_FVM\foveal\modelErrs_emily_f.mat'); modelErrs_emily_f = modelErrs;
load('C:\Users\Abhishek\Desktop\MatlabCode\Abhishek\for_FVM\foveal\modelErrs_greg_f.mat'); modelErrs_greg_f = modelErrs;
load('C:\Users\Abhishek\Desktop\MatlabCode\Abhishek\for_FVM\5_degree\modelErrs_abhi_p.mat'); modelErrs_abhi_p = modelErrs;
load('C:\Users\Abhishek\Desktop\MatlabCode\Abhishek\for_FVM\5_degree\modelErrs_emily_p.mat'); modelErrs_emily_p = modelErrs;
load('C:\Users\Abhishek\Desktop\MatlabCode\Abhishek\for_FVM\5_degree\modelErrs_greg_p.mat'); modelErrs_greg_p = modelErrs;
clear modelErrs

lowsf_abhi_par = [1 threshold_abhi_p.lowsf.son; 2 threshold_abhi_p.lowsf.soff];
lowsf_greg_par = [1 threshold_greg_p.lowsf.son; 2 threshold_greg_p.lowsf.soff];
lowsf_emily_par = [1 threshold_emily_p.lowsf.son; 2 threshold_emily_p.lowsf.soff];
highsf_abhi_par = [1 threshold_abhi_p.highsf.son; 2 threshold_abhi_p.highsf.soff];
highsf_greg_par = [1 threshold_greg_p.highsf.son; 2 threshold_greg_p.highsf.soff];
highsf_emily_par = [1 threshold_emily_p.highsf.son; 2 threshold_emily_p.highsf.soff];

lowsf_abhi_fov = [1 threshold_abhi_f.lowsf.son; 2 threshold_abhi_f.lowsf.soff];
lowsf_greg_fov = [1 threshold_greg_f.lowsf.son; 2 threshold_greg_f.lowsf.soff];
lowsf_emily_fov = [1 threshold_emily_f.lowsf.son; 2 threshold_emily_f.lowsf.soff];
highsf_abhi_fov = [1 threshold_abhi_f.highsf.son; 2 threshold_abhi_f.highsf.soff];
highsf_greg_fov = [1 threshold_greg_f.highsf.son; 2 threshold_greg_f.highsf.soff];
highsf_emily_fov = [1 threshold_emily_f.highsf.son; 2 threshold_emily_f.highsf.soff];

figure(1),subplot(221),errorbar(lowsf_abhi_fov(:,1), lowsf_abhi_fov(:,2),[modelErrs_abhi_f.low_sf(1); modelErrs_abhi_f.low_sf(3)],'-','Linewidth',2); hold on;
errorbar(lowsf_greg_fov(:,1), lowsf_greg_fov(:,2),[modelErrs_greg_f.low_sf(1); modelErrs_greg_f.low_sf(3)],'-','Linewidth',2); hold on;
errorbar(lowsf_emily_fov(:,1), lowsf_emily_fov(:,2),[modelErrs_emily_f.low_sf(1); modelErrs_emily_f.low_sf(3)],'-','Linewidth',2); hold on;
legend('AD','GH','EG'); set(gca,'Xlim',[0.5 2.5],'Ylim',[0.02 0.3],'Yscale','log'); title('Low SF'); ylabel('Contrast');hold off;

subplot(222),errorbar(highsf_abhi_fov(:,1), highsf_abhi_fov(:,2),[modelErrs_abhi_f.high_sf(1); modelErrs_abhi_f.high_sf(3)],'-','Linewidth',2); hold on;
errorbar(highsf_greg_fov(:,1), highsf_greg_fov(:,2),[modelErrs_greg_f.high_sf(1); modelErrs_greg_f.high_sf(3)],'-','Linewidth',2); hold on;
errorbar(highsf_emily_fov(:,1), highsf_emily_fov(:,2),[modelErrs_emily_f.high_sf(1); modelErrs_emily_f.high_sf(3)],'-','Linewidth',2); hold on;
legend('AD','GH','EG');set(gca,'Xlim',[0.5 2.5],'Ylim',[0.02 0.3],'Yscale','log'); title('High SF'); ylabel('Contrast');hold off;

subplot(223),errorbar(lowsf_abhi_par(:,1), lowsf_abhi_par(:,2),[modelErrs_abhi_p.low_sf(1); modelErrs_abhi_p.low_sf(3)],'-','Linewidth',2); hold on;
errorbar(lowsf_greg_par(:,1), lowsf_greg_par(:,2),[modelErrs_greg_p.low_sf(1); modelErrs_greg_p.low_sf(3)],'-','Linewidth',2); hold on;
errorbar(lowsf_emily_par(:,1), lowsf_emily_par(:,2),[modelErrs_emily_p.low_sf(1); modelErrs_emily_p.low_sf(3)],'-','Linewidth',2); hold on;
legend('AD','GH','EG'); set(gca,'Xlim',[0.5 2.5],'Ylim',[0.05 0.6],'Yscale','log'); title('Low SF'); ylabel('Contrast');hold off;

subplot(224),errorbar(highsf_abhi_par(:,1), highsf_abhi_par(:,2),[modelErrs_abhi_p.high_sf(1); modelErrs_abhi_p.high_sf(3)],'-','Linewidth',2); hold on;
errorbar(highsf_greg_par(:,1), highsf_greg_par(:,2),[modelErrs_greg_p.high_sf(1); modelErrs_greg_p.high_sf(3)],'-','Linewidth',2); hold on;
errorbar(highsf_emily_par(:,1), highsf_emily_par(:,2),[modelErrs_emily_p.high_sf(1); modelErrs_emily_p.high_sf(3)],'-','Linewidth',2); hold on;
legend('AD','GH','EG');set(gca,'Xlim',[0.5 2.5],'Ylim',[0.05 0.6],'Yscale','log'); title('High SF'); ylabel('Contrast');hold off;

