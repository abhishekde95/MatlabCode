% Data Analysis for an individual subject for all conditions
% Author: Abhishek De, 10/17

close all; clearvars;
% This is super hardcoded, someday I will come up with a better way to pull
% up the data. As of now I am loading all the files where I have used
% grating stimulus.
load('N:\Abhishek_Psychophysics_data\Rig_1_data\fitting_data\SLMTF\Abhishek\foveal_2deg\dataSLMTF.mat'); data1 = dataSLMTF; % foveal 2 deg cone fundamental data
load('N:\Abhishek_Psychophysics_data\Rig_1_data\fitting_data\SLMTF\Abhishek\foveal\dataSLMTF.mat'); data2 = dataSLMTF; % foveal 2 deg cone fundamental data
load('N:\Abhishek_Psychophysics_data\Rig_1_data\fitting_data\SLMTF\Abhishek\parafoveal\dataSLMTF.mat'); data3 = dataSLMTF; % 5 deg data with gratings 
load('N:\Abhishek_Psychophysics_data\Rig_1_data\fitting_data\SLMTF_KMwsf\Abhishek\20degdata\0.5cyclesperdeg\dataSLMTF_KMwsf.mat'); data4 = dataSLMTF_KMwsf; % 20 deg data with 0.5 cycles/deg grating
load('N:\Abhishek_Psychophysics_data\Rig_1_data\fitting_data\SLMTF_KMwsf\Abhishek\20degdata\1.0cyclesperdeg\dataSLMTF_KMwsf.mat'); data5 = dataSLMTF_KMwsf; % 20 deg data with 1.0 cycles/deg grating

% Now am loading the data where I have used gaussian blob stimulus. Here I am using lime, magenta, orange and cyan as my stimulus.
load('N:\Abhishek_Psychophysics_data\Rig_1_data\fitting_data\SLMTF_KM\Abhishek\foveal_2deg\dataSLMTF_KM.mat'); data6 = dataSLMTF_KM; % foveal data with 2 deg fundamentals
load('N:\Abhishek_Psychophysics_data\Rig_1_data\fitting_data\SLMTF_KM\Abhishek\20degdata\dataSLMTF_KM.mat');  data7 = dataSLMTF_KM; % 20 deg ecc data 

% First plotting the asymmetries of the poles along the intermediate
% direction as a function of eccentricities, plot1 - lime-magenta, plot2 - orange-cyan
figure(1),subplot(121),errorbar([data6.fit_data.modellm_ind(1) data7.fit_data.modellm_ind(1)],[data6.modelErrs.lm(1) data7.modelErrs.lm(1)],'-','Linewidth',2); hold on;
errorbar([data6.fit_data.modellm_ind(3) data7.fit_data.modellm_ind(3)],[data6.modelErrs.lm(3) data7.modelErrs.lm(3)],'-','Linewidth',2);
text(0.8,0.8, strcat('p=',num2str(1-chi2cdf(data6.lm_LLratio,1)))); text(1.8,0.8, strcat('p=',num2str(1-chi2cdf(data7.lm_LLratio,1))));
legend('lime','magenta','location','SouthEast');set(gca,'Xlim',[0.5 2.5],'XTick',[1 2],'XTickLabel',{'fovea','20 deg'},'Yscale','log'); title('lime-magenta'); ylabel('Contrast');hold off;
subplot(122),errorbar([data6.fit_data.modeloc_ind(1) data7.fit_data.modeloc_ind(1)],[data6.modelErrs.oc(1) data7.modelErrs.oc(1)],'-','Linewidth',2); hold on;
errorbar([data6.fit_data.modeloc_ind(3) data7.fit_data.modeloc_ind(3)],[data6.modelErrs.oc(3) data7.modelErrs.oc(3)],'-','Linewidth',2);
text(0.8,0.8, strcat('p=',num2str(1-chi2cdf(data6.oc_LLratio,1)))); text(1.8,0.8, strcat('p=',num2str(1-chi2cdf(data7.oc_LLratio,1))));
legend('orange','cyan','location','SouthEast');set(gca,'Xlim',[0.5 2.5],'XTick',[1 2],'XTickLabel',{'fovea','20 deg'},'Yscale','log'); title('orange-cyan'); ylabel('Contrast');hold off;

% Here I will plotting the 0.5 cycles per deg data for all the eccentricities 
figure(2), subplot(221), errorbar([data1.fit_data.model3(3) data3.fit_data.model3(3) data4.fit_data.model_ind(3)],[data1.modelErrs.low_sf(3) data3.modelErrs.low_sf(3) data4.modelErrs(3) ],'-','Linewidth',2); hold on;
errorbar([ data1.fit_data.model3(1) data3.fit_data.model3(1) data4.fit_data.model_ind(1)],[data1.modelErrs.low_sf(1) data3.modelErrs.low_sf(1) data4.modelErrs(1)],'-','Linewidth',2);
text(0.8,0.8, strcat('p=',num2str(1-chi2cdf(data1.low_sf_LLratio,1)))); text(1.8,0.8, strcat('p=',num2str(1-chi2cdf(data3.low_sf_LLratio,1)))); text(2.8,0.8, strcat('p=',num2str(1-chi2cdf(data4.LLratio,1))));
legend('lime-magenta','orange-cyan','location','SouthEast');set(gca,'Xlim',[0.5 3.5],'XTick',[1 2 3],'XTickLabel',{'fovea','5 deg','20 deg'},'Yscale','log'); title('SF = 0.5'); ylabel('Contrast');hold off;
subplot(222), errorbar([data2.fit_data.model3(3) data3.fit_data.model3(3) data4.fit_data.model_ind(3)],[data2.modelErrs.low_sf(3) data3.modelErrs.low_sf(3) data4.modelErrs(3) ],'-','Linewidth',2); hold on;
errorbar([ data2.fit_data.model3(1) data3.fit_data.model3(1) data4.fit_data.model_ind(1)],[data2.modelErrs.low_sf(1) data3.modelErrs.low_sf(1) data4.modelErrs(1)],'-','Linewidth',2);
text(0.8,0.8, strcat('p=',num2str(1-chi2cdf(data2.low_sf_LLratio,1)))); text(1.8,0.8, strcat('p=',num2str(1-chi2cdf(data3.low_sf_LLratio,1)))); text(2.8,0.8, strcat('p=',num2str(1-chi2cdf(data4.LLratio,1))));
legend('lime-magenta','orange-cyan','location','SouthEast');set(gca,'Xlim',[0.5 3.5],'XTick',[1 2 3],'XTickLabel',{'fovea','5 deg','20 deg'},'Yscale','log'); title('SF = 0.5 all 10 deg fund'); ylabel('Contrast');hold off;
% Here I am plotting the data for high spatial frequency which is 3 cycles/deg for fovea and 5 deg and 1 cycles/ deg for 20 deg ecc
subplot(223), errorbar([data1.fit_data.model4(3) data3.fit_data.model4(3) data5.fit_data.model_ind(3)],[data1.modelErrs.high_sf(3) data3.modelErrs.high_sf(3) data5.modelErrs(3) ],'-','Linewidth',2); hold on;
errorbar([ data1.fit_data.model4(1) data3.fit_data.model4(1) data5.fit_data.model_ind(1)],[data1.modelErrs.high_sf(1) data3.modelErrs.high_sf(1) data5.modelErrs(1)],'-','Linewidth',2);
text(0.8,0.8, strcat('p=',num2str(1-chi2cdf(data1.high_sf_LLratio,1)))); text(1.8,0.8, strcat('p=',num2str(1-chi2cdf(data3.high_sf_LLratio,1)))); text(2.8,0.8, strcat('p=',num2str(1-chi2cdf(data5.LLratio,1))));
legend('lime-magenta','orange-cyan','location','SouthEast');set(gca,'Xlim',[0.5 3.5],'XTick',[1 2 3],'XTickLabel',{'fovea','5 deg','20 deg'},'Yscale','log'); title('high SF'); ylabel('Contrast');hold off;
subplot(224), errorbar([data2.fit_data.model4(3) data3.fit_data.model4(3) data5.fit_data.model_ind(3)],[data1.modelErrs.high_sf(3) data3.modelErrs.high_sf(3) data5.modelErrs(3) ],'-','Linewidth',2); hold on;
errorbar([ data2.fit_data.model4(1) data3.fit_data.model4(1) data5.fit_data.model_ind(1)],[data1.modelErrs.high_sf(1) data3.modelErrs.high_sf(1) data5.modelErrs(1)],'-','Linewidth',2);
text(0.8,0.8, strcat('p=',num2str(1-chi2cdf(data2.high_sf_LLratio,1)))); text(1.8,0.8, strcat('p=',num2str(1-chi2cdf(data3.high_sf_LLratio,1)))); text(2.8,0.8, strcat('p=',num2str(1-chi2cdf(data5.LLratio,1))));
legend('lime-magenta','orange-cyan','location','SouthEast');set(gca,'Xlim',[0.5 3.5],'XTick',[1 2 3],'XTickLabel',{'fovea','5 deg','20 deg'},'Yscale','log'); title('high SF all 10 deg fund'); ylabel('Contrast');hold off;

% Till now I have been plotting detection thresholds as a function of
% eccentricity. Now I am gonna plot detection thresholds as a function of
% spatial frequency

figure(3),subplot(141),errorbar([data6.fit_data.model_ind(3) data1.fit_data.model3(3) data1.fit_data.model4(3)],[data6.modelErrs.model_ind(3) data1.modelErrs.low_sf(3) data1.modelErrs.high_sf(3)],'-','Linewidth',2); hold on;
errorbar([data6.fit_data.model_ind(1) data1.fit_data.model3(1) data1.fit_data.model4(1)],[data6.modelErrs.model_ind(1) data1.modelErrs.low_sf(1) data1.modelErrs.high_sf(1)],'-','Linewidth',2);
text(0.8,0.8, strcat('p=',num2str(1-chi2cdf(data6.LLratio,1)))); text(1.8,0.8, strcat('p=',num2str(1-chi2cdf(data1.low_sf_LLratio,1)))); text(2.8,0.8, strcat('p=',num2str(1-chi2cdf(data1.high_sf_LLratio,1))));
legend('lime-magenta','orange-cyan','location','SouthEast');set(gca,'Xlim',[0.5 3.5],'Ylim',[0.01 1],'XTick',[1 2 3],'XTickLabel',{'0','0.5','3'},'Yscale','log'); title('Foveal data 2 deg fund'); 
xlabel('SF'); ylabel('Contrast');hold off;
subplot(142),errorbar([data2.fit_data.model3(3) data2.fit_data.model4(3)],[data2.modelErrs.low_sf(3) data2.modelErrs.high_sf(3)],'-','Linewidth',2); hold on;
errorbar([data2.fit_data.model3(1) data2.fit_data.model4(1)],[data2.modelErrs.low_sf(1) data2.modelErrs.high_sf(1)],'-','Linewidth',2);
text(0.8,0.8, strcat('p=',num2str(1-chi2cdf(data2.low_sf_LLratio,1)))); text(1.8,0.8, strcat('p=',num2str(1-chi2cdf(data2.high_sf_LLratio,1))));
legend('lime-magenta','orange-cyan','location','SouthEast');set(gca,'Xlim',[0.5 2.5],'Ylim',[0.01 1],'XTick',[1 2],'XTickLabel',{'0.5','3'},'Yscale','log'); title('Foveal data'); 
xlabel('SF'); ylabel('Contrast');hold off;
% 5 deg data
subplot(143),errorbar([data3.fit_data.model3(3) data3.fit_data.model4(3)],[data3.modelErrs.low_sf(3) data3.modelErrs.high_sf(3)],'-','Linewidth',2); hold on;
errorbar([data3.fit_data.model3(1) data3.fit_data.model4(1)],[data3.modelErrs.low_sf(1) data3.modelErrs.high_sf(1)],'-','Linewidth',2);
text(0.8,0.8, strcat('p=',num2str(1-chi2cdf(data3.low_sf_LLratio,1)))); text(1.8,0.8, strcat('p=',num2str(1-chi2cdf(data3.high_sf_LLratio,1))));
legend('lime-magenta','orange-cyan','location','SouthEast');set(gca,'Xlim',[0.5 2.5],'XTick',[1 2 3],'XTickLabel',{'0.5','3'},'Yscale','log'); title('5 deg data'); 
xlabel('SF'); ylabel('Contrast');hold off;
% 20 deg data
subplot(144),errorbar([data7.fit_data.model_ind(3) data4.fit_data.model_ind(3) data5.fit_data.model_ind(3)],[data7.modelErrs.model_ind(3) data4.modelErrs(3) data5.modelErrs(3)],'-','Linewidth',2); hold on;
errorbar([data7.fit_data.model_ind(1) data4.fit_data.model_ind(1) data5.fit_data.model_ind(1)],[data7.modelErrs.model_ind(1) data5.modelErrs(1) data5.modelErrs(1)],'-','Linewidth',2);
text(0.8,0.8, strcat('p=',num2str(1-chi2cdf(data7.LLratio,1)))); text(1.8,0.8, strcat('p=',num2str(1-chi2cdf(data4.LLratio,1)))); text(2.8,0.8, strcat('p=',num2str(1-chi2cdf(data5.LLratio,1))));
legend('lime-magenta','orange-cyan','location','SouthEast');set(gca,'Xlim',[0.5 3.5],'Ylim',[0.01 1],'XTick',[1 2 3],'XTickLabel',{'0','0.5','1.0'},'Yscale','log'); title('20 deg data'); 
xlabel('SF'); ylabel('Contrast');hold off;
